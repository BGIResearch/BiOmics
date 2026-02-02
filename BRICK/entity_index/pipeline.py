"""
功能：将实体 catalog、字符串索引设置与向量嵌入生成流程整合为统一工件，供外部写入 ES 或其他存储。
构建逻辑：遍历节点配置，生成字符串索引文档与 settings，按需调用嵌入函数补充向量文档，并封装为 `IndexArtifacts`。
数据流：上游依赖于 `catalog.build_catalog` 和可选的外部嵌入函数；下游则是写入模块（如 `es_writer`）或检索模块使用生成的索引 artificats。
调用链：典型流程为 `build_index_artifacts` -> （外部写入器），混合检索模块也可重用其输出构建索引。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Optional

from .catalog import build_catalog
from .config import NODE_CONFIGS, NodeConfig
from .index_builders import (
    dictionary_document,
    dictionary_index_settings,
    string_document,
    string_index_settings,
    vector_document,
    vector_index_settings,
)
from .entity_types import EntityRecord
from .loader import load_dictionary_records


EmbeddingFn = Callable[[List[EntityRecord], str], List[List[float]]]


@dataclass
class IndexArtifacts:
    node: str
    config: NodeConfig
    string_settings: Dict[str, Any]
    string_documents: List[Dict[str, Any]]
    vector_settings: Optional[Dict[str, Any]] = None
    vector_documents: Optional[List[Dict[str, Any]]] = None


def build_index_artifacts(
    node_names: Optional[Iterable[str]] = None,
    embedding_fn: Optional[EmbeddingFn] = None,
    embedding_dimension: int = 1024,
    include_vectors: bool = True,
    catalog_sample_limit: Optional[int] = None,
) -> List[IndexArtifacts]:
    """
    功能：汇总实体 catalog、字符串/向量索引配置，生成 `IndexArtifacts` 列表供写入模块使用。
    构建逻辑：
    1. 调用 `build_catalog` 获取节点 -> 实体记录；
    2. 为每个节点构造字符串索引文档与 settings；
    3. 若节点启用向量且允许，使用 `embedding_fn` 生成 embedding，并构造向量文档与 settings；
    4. 将结果封装进 `IndexArtifacts`。
    出参入参：入参可指定节点集合、嵌入函数、向量维度、是否包含向量、可选的 catalog 限制；出参为 `IndexArtifacts` 列表。
    数据流：上游索引构建脚本/服务调用本方法；下游 `es_writer` 或其他写入器消费工件写入 ES。
    模块内调用链：依赖 `build_catalog`、`string_document`、`vector_document`、`vector_index_settings` 等。
    模块外调用链：`es_writer`、混合检索初始化流程等使用输出结果。
    """
    selected_nodes = list(node_names) if node_names else list(NODE_CONFIGS.keys())
    entity_nodes = [name for name in selected_nodes if not NODE_CONFIGS[name].is_dictionary]
    catalog = build_catalog(entity_nodes)
    artifacts: List[IndexArtifacts] = []

    for node in selected_nodes:
        config = NODE_CONFIGS[node]

        if config.is_dictionary:
            dict_records = load_dictionary_records(config)
            if catalog_sample_limit and catalog_sample_limit > 0:
                dict_records = dict_records[:catalog_sample_limit]
            string_docs = [dictionary_document(record) for record in dict_records]
            artifact = IndexArtifacts(
                node=node,
                config=config,
                string_settings=dictionary_index_settings(),
                string_documents=string_docs,
            )
            artifacts.append(artifact)
            continue

        records = catalog.get(node, [])
        if catalog_sample_limit and catalog_sample_limit > 0:
            # TEMP: restrict catalog size during embedding experiments; revert after testing.
            records = records[:catalog_sample_limit]
        string_docs = [string_document(record) for record in records]

        artifact = IndexArtifacts(
            node=node,
            config=config,
            string_settings=string_index_settings(),
            string_documents=string_docs,
        )

        if config.enable_vector and include_vectors:
            if embedding_fn is None:
                raise ValueError(
                    f"Node {node} requires vector embeddings but no embedding_fn was provided"
                )
            vectors = embedding_fn(records, node)
            if len(vectors) != len(records):
                raise ValueError(
                    f"Embedding function returned {len(vectors)} vectors for {len(records)} records"
                )
            vector_docs = [
                vector_document(record, list(vector))
                for record, vector in zip(records, vectors)
            ]
            artifact.vector_settings = vector_index_settings(embedding_dimension)
            artifact.vector_documents = vector_docs

        artifacts.append(artifact)

    return artifacts


__all__ = ["IndexArtifacts", "build_index_artifacts"]
