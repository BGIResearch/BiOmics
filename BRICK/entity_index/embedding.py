"""
功能：封装嵌入模型客户端与批量调用逻辑，为实体索引构建或检索提供标准化的向量生成接口。
构建逻辑：定义配置数据类、OpenAI 兼容客户端以及批处理函数，并通过环境变量读取鉴权及模型信息。
数据流：上游多为索引构建流程或混合检索模块需要向量时调用；下游是 Elasticsearch 等存储或检索模块消费生成的 embedding。
调用链：`pipeline.build_index_artifacts` 调用 `build_batch_embedder` 组合 `EmbeddingClient` 完成向量生成；混合检索亦可复用这些接口。
"""

from __future__ import annotations

import os
from dataclasses import dataclass
import logging
import time
from typing import Callable, List, Sequence

try:  # optional dependency
    from openai import OpenAI  # type: ignore
except ImportError as exc:  # pragma: no cover
    OpenAI = None  # type: ignore
    _IMPORT_ERROR = exc
else:
    _IMPORT_ERROR = None

from .entity_types import EntityRecord


@dataclass
class EmbeddingConfig:
    api_key: str
    base_url: str
    model: str = "text-embedding-v3"
    timeout: float = 300.0

    @classmethod
    def from_env(cls) -> "EmbeddingConfig":
        """
        功能：从环境变量读取嵌入服务配置，生成 `EmbeddingConfig` 实例。
        构建逻辑：依次读取 `EMBEDDING_API_KEY`、`EMBEDDING_BASE_URL`、`EMBEDDING_MODEL`、`EMBEDDING_TIMEOUT`，校验必需项并设置默认模型与超时时间。
        出参入参：无入参；返回 `EmbeddingConfig` 对象。
        数据流：上游在初始化嵌入客户端或批处理器时调用；下游为 `EmbeddingClient`、`build_batch_embedder` 等模块。
        模块内调用链：`EmbeddingClient.__init__`（可调用本方法）。
        模块外调用链：索引构建/检索模块在准备嵌入配置时调用。
        """
        api_key = os.environ.get("EMBEDDING_API_KEY")
        base_url = os.environ.get("EMBEDDING_BASE_URL")
        model = os.environ.get("EMBEDDING_MODEL", "text-embedding-v3")
        timeout_raw = os.environ.get("EMBEDDING_TIMEOUT")
        timeout = float(timeout_raw) if timeout_raw else 300.0

        missing = [name for name, value in [
            ("EMBEDDING_API_KEY", api_key),
            ("EMBEDDING_BASE_URL", base_url),
        ] if not value]
        if missing:
            raise RuntimeError(f"Missing embedding environment variables: {', '.join(missing)}")
        return cls(api_key=api_key, base_url=base_url, model=model, timeout=timeout)


class EmbeddingClient:
    def __init__(self, config: EmbeddingConfig):
        """
        功能：根据配置初始化 OpenAI 兼容的嵌入客户端，供文本向量化调用。
        构建逻辑：确保 OpenAI SDK 可用，使用配置中的 key、base_url、timeout 创建客户端，并记录模型名。
        出参入参：入参为 `EmbeddingConfig`；无返回值。
        数据流：上游根据配置生成客户端；下游调用 `embed_texts`/`embed_text` 获取 embedding。
        模块内调用链：无。
        模块外调用链：`build_batch_embedder`、其他嵌入需求逻辑。
        """
        if OpenAI is None:  # pragma: no cover
            raise RuntimeError(
                "openai package is required for embeddings but is not installed"
            ) from _IMPORT_ERROR
        self._client = OpenAI(
            api_key=config.api_key,
            base_url=config.base_url,
            timeout=config.timeout,
        )
        self._model = config.model

    def embed_texts(self, texts: Sequence[str]) -> List[List[float]]:
        """
        功能：批量请求嵌入向量，与输入文本一一对应。
        构建逻辑：调用 OpenAI embeddings API，将文本列表转换为 embedding 列表。
        出参入参：入参为字符串序列；出参为二维浮点列表。
        数据流：`build_batch_embedder` 等批处理逻辑调用；输出用于索引或检索。
        模块内调用链：`embed_text` 调用此方法。
        模块外调用链：`build_batch_embedder`、实时检索等使用。
        """
        response = self._client.embeddings.create(model=self._model, input=list(texts))
        return [item.embedding for item in response.data]

    def embed_text(self, text: str) -> List[float]:
        """
        功能：生成单个文本的 embedding。
        构建逻辑：复用 `embed_texts`，传入单元素列表并返回第一项。
        出参入参：入参为字符串；出参为浮点向量。
        数据流：供单条向量需求使用。
        模块内调用链：调用 `embed_texts`。
        模块外调用链：上游模块直接调用获取单条 embedding。
        """
        return self.embed_texts([text])[0]


def record_to_text(record: EntityRecord, max_aliases: int = 10) -> str:
    """
    功能：将 `EntityRecord` 转换为嵌入模型输入文本，统一模板。
    构建逻辑：拼接主名称、定义及限定数量的别名，形成格式化描述。
    出参入参：入参为实体记录与最大别名数；出参为描述字符串。
    数据流：`build_batch_embedder` 调用后将文本传给嵌入客户端。
    模块内调用链：`build_batch_embedder` -> `record_to_text`。
    模块外调用链：无。
    """
    aliases = record.aliases[:max_aliases]
    alias_str = ", ".join(aliases) if aliases else "None"
    definition = record.definition or ""
    return f"{record.primary_name}. {definition}. Aliases: {alias_str}."
def build_batch_embedder(
    client: EmbeddingClient,
    batch_size: int = 10,
    max_retries: int = 3,
    retry_backoff: float = 2.0,
) -> Callable[[List[EntityRecord], str], List[List[float]]]:
    logger = logging.getLogger(__name__)
    logger.info("Using batch_size=%d for embedding requests", batch_size)

    def embed_records(records: List[EntityRecord], node: str) -> List[List[float]]:
        """
        功能：批量为实体记录生成 embedding，包含重试与日志记录。
        构建逻辑：按 `batch_size` 切分记录，调用 `record_to_text` 构造文本并向客户端请求 embedding；失败时指数退避重试，超出次数抛错。
        出参入参：入参为实体记录列表及节点名称；出参为 embedding 列表。
        数据流：`pipeline.build_index_artifacts` 和混合检索在生成向量时调用；输出供索引写入。
        模块内调用链：`embed_records` -> `record_to_text`、`EmbeddingClient.embed_texts`。
        模块外调用链：`pipeline`、检索模块使用返回的函数进行批量嵌入。
        """
        embeddings: List[List[float]] = []
        total = len(records)
        for start in range(0, total, batch_size):
            batch = records[start : start + batch_size]
            texts = [record_to_text(record) for record in batch]
            attempt = 0
            while True:
                try:
                    vectors = client.embed_texts(texts)
                    embeddings.extend(vectors)
                    logger.info(
                        "Embedded batch %d-%d/%d for node %s",
                        start + 1,
                        start + len(batch),
                        total,
                        node,
                    )
                    break
                except Exception as exc:  # noqa: BLE001
                    is_timeout = isinstance(exc, TimeoutError)
                    attempt += 1
                    logger.error(
                        "Embedding batch %d-%d for node %s failed (attempt %d/%d): %s",
                        start + 1,
                        start + len(batch),
                        node,
                        attempt,
                        max_retries,
                        exc,
                    )
                    if attempt >= max_retries:
                        logger.error(
                            "Embedding aborted after %d attempts for node %s batch %d-%d",
                            attempt,
                            node,
                            start + 1,
                            start + len(batch),
                        )
                        raise
                    time.sleep(retry_backoff * attempt)
        return embeddings

    return embed_records


__all__ = [
    "EmbeddingConfig",
    "EmbeddingClient",
    "record_to_text",
    "build_batch_embedder",
]
