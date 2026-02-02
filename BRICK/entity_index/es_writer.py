from __future__ import annotations

import os
import json
from dataclasses import dataclass
import logging
from typing import Iterable, List, Optional, Sequence, Union

from elasticsearch import Elasticsearch, helpers, exceptions  # 添加 exceptions 导入

from .pipeline import IndexArtifacts


@dataclass
class ESConfig:
    host: Union[str, Sequence[str]]
    username: Optional[str] = None
    password: Optional[str] = None
    verify_certs: bool = True

    @classmethod
    def from_env(cls) -> "ESConfig":
        """
        功能：从环境变量读取 Elasticsearch 连接配置，返回 `ESConfig` 实例。
        构建逻辑：优先解析 `ES_CONFIG`（JSON 字典），其次读取 `ES_HOSTS` 或 `ES_HOST`，以及 `ES_USERNAME`、`ES_PASSWORD`、`ES_VERIFY_CERTS` 等变量，设置默认值后封装为配置对象。
        出参入参：无入参；返回 `ESConfig`。
        数据流：上游在初始化 ES 客户端前调用此方法；下游由 `create_client` 创建真正的客户端。
        模块内调用链：`create_client` 可依赖本方法获取配置。
        模块外调用链：索引写入脚本或服务。
        """
        config_dict: dict[str, object] = {}
        raw_config = os.environ.get("ES_CONFIG")
        if raw_config:
            try:
                config_dict = json.loads(raw_config)
            except json.JSONDecodeError:
                logging.getLogger(__name__).warning(
                    "Invalid ES_CONFIG JSON, falling back to individual environment variables.",
                )

        def _coerce_hosts(value: Optional[object]) -> List[str]:
            if not value:
                return []
            if isinstance(value, str):
                items = [value]
            else:
                try:
                    items = list(value)  # type: ignore[arg-type]
                except TypeError:
                    items = []
            hosts: List[str] = []
            for item in items:
                if isinstance(item, str) and item.strip():
                    hosts.append(item.strip())
            return hosts

        hosts: List[str] = []
        hosts = _coerce_hosts(config_dict.get("hosts"))
        if not hosts:
            hosts_env = os.environ.get("ES_HOSTS")
            if hosts_env:
                hosts = [item.strip() for item in hosts_env.split(",") if item.strip()]
        if not hosts:
            single_host = config_dict.get("host") or os.environ.get("ES_HOST")
            if isinstance(single_host, str) and single_host.strip():
                hosts = [single_host.strip()]
        if not hosts:
            hosts = ["http://localhost:9200"]

        username = config_dict.get("username", os.environ.get("ES_USERNAME"))
        username = username if isinstance(username, str) and username.strip() else None

        password = config_dict.get("password", os.environ.get("ES_PASSWORD"))
        password = password if isinstance(password, str) and password.strip() else None

        verify_value = config_dict.get("verify_certs", os.environ.get("ES_VERIFY_CERTS"))
        if isinstance(verify_value, bool):
            verify = verify_value
        else:
            verify = str(verify_value).lower() != "false" if verify_value is not None else True

        host_value: Union[str, Sequence[str]] = hosts
        if len(hosts) == 1:
            host_value = hosts[0]

        return cls(host=host_value, username=username, password=password, verify_certs=verify)

    def create_client(self) -> Elasticsearch:
        """
        功能：基于配置创建 Elasticsearch 客户端实例。
        构建逻辑：设置证书校验参数；如配置了用户名/密码，则注入基础认证。
        出参入参：无额外入参；返回 `Elasticsearch` 客户端。
        数据流：写入函数使用此客户端执行索引操作。
        模块内调用链：`index_string_documents`、`index_vector_documents` 使用该客户端。
        模块外调用链：索引写入脚本。
        """
        kwargs = {
            "verify_certs": self.verify_certs,
            "headers": {
                "accept": "application/vnd.elasticsearch+json; compatible-with=8",
                "content-type": "application/vnd.elasticsearch+json; compatible-with=8",
            },
        }
        if self.username and self.password:
            kwargs.update({"basic_auth": (self.username, self.password)})
        return Elasticsearch(
            self.host,
            request_timeout=180,
            retry_on_timeout=True,
            max_retries=3,
            **kwargs,
        )


def ensure_index(es: Elasticsearch, index_name: str, settings: dict, recreate: bool = False) -> None:
    """
    功能：确保目标索引存在；在需要重建时删除旧索引再按 settings 创建。
    构建逻辑：调用 ES indices API 检查索引是否存在，结合 `recreate` 标志执行 delete/create 操作。
    出参入参：入参为 ES 客户端、索引名、settings、重建标志；无返回值。
    数据流：字符串/向量写入前的准备步骤。
    模块内调用链：`index_string_documents`、`index_vector_documents`。
    模块外调用链：无。
    """
    exists_kwargs = {"request_timeout": 180}
    logging.getLogger(__name__).info(f"Ensuring index: {index_name}, recreate: {recreate}, exists_kwargs: {exists_kwargs}")
    logging.getLogger(__name__).info(f"Index settings: {json.dumps(settings, indent=2)}")  # 打印索引设置
    def _index_exists() -> bool:
        try:
            return es.indices.exists(index=index_name, **exists_kwargs)  # type: ignore[attr-defined]
        except exceptions.BadRequestError as exc:
            status = getattr(exc.meta, "status", None) if hasattr(exc, "meta") else None
            if status == 400:
                logger = logging.getLogger(__name__)
                logger.warning(
                    "indices.exists returned 400 for %s; falling back to cat.indices for existence check",
                    index_name,
                )
                compat_headers = {
                    "accept": "application/vnd.elasticsearch+json; compatible-with=8",
                    "content-type": "application/vnd.elasticsearch+json; compatible-with=8",
                }
                try:
                    cat_response = es.cat.indices(  # type: ignore[attr-defined]
                        index=index_name,
                        format="json",
                        headers=compat_headers,
                    )
                except exceptions.BadRequestError as cat_exc:
                    logger.error("cat.indices also failed for %s: %s", index_name, cat_exc)
                    raise
                return bool(cat_response)
            raise

    try:
        exists = _index_exists()
        if recreate and exists:
            es.indices.delete(index=index_name, **exists_kwargs)  # type: ignore[attr-defined]
            exists = False
        if not exists:
            es.indices.create(  # type: ignore[attr-defined]
                index=index_name,
                body=settings,
                timeout="180s",
                master_timeout="180s"
            )
    except exceptions.BadRequestError as e:  # 捕获详细错误信息
        logger = logging.getLogger(__name__)
        logger.error("Failed to ensure index %s: %s", index_name, e)
        try:
            logger.error("Error body: %s", json.dumps(e.body, indent=2) if e.body else "No body")
        except TypeError:
            logger.error("Error body (non-JSON): %r", e.body)
        if e.meta:
            meta_payload = getattr(e.meta, "raw", None) or repr(e.meta)
            logger.error("Error meta: %s", meta_payload)
        else:
            logger.error("Error meta: No meta")
        raise


def _index_name(node: str, kind: str, prefix: Optional[str] = None) -> str:
    """
    功能：根据节点名称与索引类型生成统一的索引名称，例如 `entity_gene_string`。
    构建逻辑：读取可选的索引前缀（函数参数 > `ES_INDEX_PREFIX` 环境变量 > 默认 `entity`），将节点名转换为小写并与类型字符串拼接。
    出参入参：入参为节点名、索引类型；出参为字符串索引名。
    数据流：字符串/向量写入函数使用该名称定位目标索引。
    模块内调用链：`index_string_documents`、`index_vector_documents`。
    模块外调用链：无。
    """
    base_prefix = prefix or os.environ.get("ES_INDEX_PREFIX") or "entity"
    base_prefix = base_prefix.strip() or "entity"
    index_name = f"{base_prefix}_{node.lower()}_{kind}"
    logging.getLogger(__name__).info(f"Generated index name: {index_name}")  # 添加日志
    return index_name


def _maybe_limit(documents: List[dict], limit: Optional[int]) -> List[dict]:
    """
    功能：根据 `limit` 参数对文档列表进行截断，便于抽样或调试。
    构建逻辑：当未指定限制或限制大于等于文档数时直接返回原列表，否则切取前 `limit` 条。
    出参入参：入参为文档列表和可选的限制值；出参为处理后的列表。
    数据流：写入函数在批量入库前调用以控制写入数量。
    模块内调用链：`index_string_documents`、`index_vector_documents`。
    模块外调用链：无。
    """
    if limit is None or limit >= len(documents):
        logging.getLogger(__name__).info(f"Loaded all documents: {len(documents)}")  # 添加日志
        return documents
    logging.getLogger(__name__).info(f"Loaded limited documents: {limit}")  # 添加日志
    return documents[:limit]


def _actions(index_name: str, documents: Iterable[dict]):
    """
    功能：生成 Elasticsearch bulk API 所需的动作/文档结构。
    构建逻辑：遍历文档集合，逐条产出 `{\"_index\": index_name, \"_source\": doc}`。
    出参入参：入参为索引名和文档迭代器；出参为可迭代的 bulk 动作。
    数据流：供 `helpers.bulk` 调用执行批量写入。
    模块内调用链：`index_string_documents`、`index_vector_documents`。
    模块外调用链：无。
    """
    for doc in documents:
        yield {"_index": index_name, "_source": doc}


def _chunked(documents: List[dict], chunk_size: int) -> Iterable[List[dict]]:
    """
    功能：按批大小切分文档列表，便于控制 bulk 写入规模。
    构建逻辑：使用 range 步进切片列表，逐批返回子列表。
    出参入参：入参为文档列表与批大小；出参为批次迭代器。
    数据流：写入函数调用本工具逐批处理文档。
    模块内调用链：`index_string_documents`、`index_vector_documents`。
    模块外调用链：无。
    """
    for start in range(0, len(documents), chunk_size):
        yield documents[start : start + chunk_size]


def index_string_documents(
    es: Elasticsearch,
    artifact: IndexArtifacts,
    recreate: bool = False,
    doc_limit: Optional[int] = None,
    chunk_size: int = 500,
    index_prefix: Optional[str] = None,
) -> str:
    """
    功能：批量写入字符串索引文档，可控制索引重建、写入数量与批大小。
    构建逻辑：确定索引名并调用 `ensure_index`，按需截断文档列表，利用 `_chunked` 和 `_actions` 搭配 `helpers.bulk` 写入，过程中记录进度日志。
    出参入参：入参包含 ES 客户端、索引工件、重建标志、文档限制、批大小；出参为索引名。
    数据流：输出索引用于字符串检索通道。
    模块内调用链：`_index_name`、`ensure_index`、`_maybe_limit`、`_chunked`、`_actions`。
    模块外调用链：`index_artifacts`。
    """
    index_name = _index_name(artifact.node, "string", prefix=index_prefix)
    ensure_index(es, index_name, artifact.string_settings, recreate=recreate)
    docs = _maybe_limit(artifact.string_documents, doc_limit)
    logger = logging.getLogger(__name__)
    total = len(docs)
    if total == 0:
        logger.info("No string documents to index for node %s", artifact.node)
        return index_name
    processed = 0
    for chunk in _chunked(docs, chunk_size):
        helpers.bulk(  # type: ignore[attr-defined]
            es,
            _actions(index_name, chunk),
            chunk_size=chunk_size,
        )
        processed += len(chunk)
        logger.info(
            "Indexed %d/%d string documents for node %s",
            processed,
            total,
            artifact.node,
        )
    return index_name


def index_vector_documents(
    es: Elasticsearch,
    artifact: IndexArtifacts,
    recreate: bool = False,
    doc_limit: Optional[int] = None,
    chunk_size: int = 200,
    index_prefix: Optional[str] = None,
) -> Optional[str]:
    """
    功能：批量写入向量索引文档，若工件未包含向量数据则返回 None。
    构建逻辑：当存在向量文档与 settings 时，确保索引可用，按需截断文档后以批量方式写入，并记录进度日志。
    出参入参：入参为 ES 客户端、索引工件、重建标志、写入限制、批大小；出参为索引名或 None。
    数据流：写入完成后的索引用于向量/语义检索。
    模块内调用链：复用 `_index_name`、`ensure_index`、`_maybe_limit`、`_chunked`、`_actions`。
    模块外调用链：`index_artifacts`。
    """
    if not artifact.vector_documents or not artifact.vector_settings:
        return None
    index_name = _index_name(artifact.node, "vector", prefix=index_prefix)
    ensure_index(es, index_name, artifact.vector_settings, recreate=recreate)
    docs = _maybe_limit(artifact.vector_documents, doc_limit)
    logger = logging.getLogger(__name__)
    total = len(docs)
    if total == 0:
        logger.info("No vector documents to index for node %s", artifact.node)
        return index_name
    processed = 0
    for chunk in _chunked(docs, chunk_size):
        try:
            helpers.bulk(  # type: ignore[attr-defined]
                es,
                _actions(index_name, chunk),
                chunk_size=chunk_size,
                request_timeout=180,
            )
        except helpers.BulkIndexError as exc:  # type: ignore[attr-defined]
            samples = exc.errors[:5]
            for sample in samples:
                logger.error("Bulk index error sample for node %s: %s", artifact.node, sample)
            logger.error(
                "Bulk indexing failed for node %s; %d document(s) failed.",
                artifact.node,
                len(exc.errors),
            )
            raise
        processed += len(chunk)
        logger.info(
            "Indexed vector batch %d/%d for node %s",
            processed,
            total,
            artifact.node,
        )
    return index_name


def index_artifacts(
    es: Elasticsearch,
    artifacts: Iterable[IndexArtifacts],
    recreate: bool = False,
    sample_limit: Optional[int] = None,
    string_chunk_size: int = 500,
    vector_chunk_size: int = 200,
    index_prefix: Optional[str] = None,
) -> None:
    """
    功能：遍历索引工件，依次写入字符串与向量索引，实现完整的索引部署。
    构建逻辑：对每个工件先调用 `index_string_documents` 再调用 `index_vector_documents`，支持重建、抽样与批大小配置。
    出参入参：入参为 ES 客户端、工件迭代器、重建标志、采样上限、字符串/向量批大小；无返回值。
    数据流：索引构建链路的最后一步，写入结果供在线检索使用。
    模块内调用链：调用 `index_string_documents`、`index_vector_documents`。
    模块外调用链：索引部署脚本或自动化流程。
    """
    for artifact in artifacts:
        index_string_documents(
            es,
            artifact,
            recreate=recreate,
            doc_limit=sample_limit,
            chunk_size=string_chunk_size,
            index_prefix=index_prefix,
        )
        index_vector_documents(
            es,
            artifact,
            recreate=recreate,
            doc_limit=sample_limit,
            chunk_size=vector_chunk_size,
            index_prefix=index_prefix,
        )


__all__ = [
    "ESConfig",
    "ensure_index",
    "index_string_documents",
    "index_vector_documents",
    "index_artifacts",
]
