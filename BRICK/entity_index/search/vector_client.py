"""
功能：封装基于 embedding 的语义检索流程，包括向量生成、ES knn/script_score 查询，以及结果归一化，供混合检索的向量通道使用。
构建逻辑：读取 `settings` 提供的 embedding 配置，调用通用嵌入客户端生成向量，再构造 ES 查询并将命中转换为内部结构。
数据流：上游为 `NormalizedQuery` 与 `HybridSearchConfig`；下游与字符串通道一起进入 `fusion`。
调用链：`HybridEntitySearcher` 在语义检索阶段调用本模块函数获取向量候选。
"""

from __future__ import annotations

import logging
from typing import Dict, Iterable, List, Optional, Sequence

from elasticsearch import Elasticsearch  # type: ignore

from .schema import ChannelScores, HybridHit, NormalizedQuery
from .settings import HybridSearchConfig
from ..embedding import EmbeddingClient

TYPE_KEY_TO_CTYPES: Dict[str, Sequence[str]] = {
    "Gene|Protein": ("Gene", "Protein"),
    "Mutation": ("Mutation",),
    "Chemical": ("Chemical",),
    "Disease|Phenotype": ("Disease", "Phenotype"),
    "Process|Function|Pathway|Cell_Component": (
        "Process",
        "Function",
        "Pathway",
        "Cell_Component",
    ),
    "Species": ("Species",),
    "Cell|Tissue": ("Cell", "Tissue"),
}

logger = logging.getLogger(__name__)


def _build_query_embedding(
    client: EmbeddingClient,
    texts: Iterable[str],
) -> List[float]:
    """
    功能：调用嵌入客户端生成查询向量，通常针对单条合并文本。
    构建逻辑：将候选词拼接为文本后调用客户端生成 embedding，返回单个向量。
    出参入参：入参为嵌入客户端与待转换文本序列；出参为向量列表（单条）。
    数据流：`search_vector_channel` 在发起 ES 检索前调用本函数。
    模块内调用链：内部辅助使用。
    模块外调用链：无。
    """
    parts = [text.strip() for text in texts if isinstance(text, str) and text.strip()]
    if not parts:
        return []
    query_text = " ; ".join(parts)
    return client.embed_text(query_text)


def _search_type(
    es: Elasticsearch,
    index_name: str,
    type_key: str,
    embedding: List[float],
    size: int,
) -> List[HybridHit]:
    """
    功能：在向量索引上执行语义检索，返回同一类型的向量命中列表。
    构建逻辑：优先发起 knn 查询，若 ES 版本不支持则回退到 script_score，并将命中转换为 `HybridHit`，仅填充向量得分。
    出参入参：入参为 ES 客户端、索引名、类型键、查询向量、返回数量；出参为 `HybridHit` 列表。
    数据流：`search_vector_channel` 对每个类型调用本函数。
    模块内调用链：无（独立使用 embedding）。
    模块外调用链：无。
    """
    if not embedding:
        return []

    filters: List[Dict[str, object]] = []
    ctype_filters = TYPE_KEY_TO_CTYPES.get(type_key)
    if ctype_filters:
        filters.append({"terms": {"ctype": list(ctype_filters)}})

    knn_body: Dict[str, object] = {
        "size": size,
        "track_total_hits": False,
        "knn": {
            "field": "embedding",
            "query_vector": embedding,
            "k": size,
            "num_candidates": max(size * 4, size),
        },
    }
    if filters:
        knn_body["knn"]["filter"] = {"bool": {"filter": filters}}

    try:
        response = es.search(index=index_name, body=knn_body)
    except Exception as exc:  # noqa: BLE001
        logger.warning(
            "Vector knn search failed for type %s on index %s, fallback to script_score: %s",
            type_key,
            index_name,
            exc,
        )
        query: Dict[str, object] = {"match_all": {}}
        if filters:
            query = {"bool": {"filter": filters}}
        body = {
            "size": size,
            "track_total_hits": False,
            "query": {
                "script_score": {
                    "query": query,
                    "script": {
                        "source": "cosineSimilarity(params.query_vector, 'embedding') + 1.0",
                        "params": {"query_vector": embedding},
                    },
                }
            },
        }
        response = es.search(index=index_name, body=body)

    hits: List[HybridHit] = []
    for item in response.get("hits", {}).get("hits", []):  # type: ignore[assignment]
        source = item.get("_source", {}) or {}
        entity_id = source.get("entity_id")
        primary_name = source.get("primary_name")
        if not entity_id or not primary_name:
            continue
        score = float(item.get("_score") or 0.0)
        hits.append(
            HybridHit(
                entity_id=entity_id,
                primary_name=primary_name,
                type_key=type_key,
                node_type=source.get("ctype", type_key),
                scores=ChannelScores(string_score=0.0, vector_score=score),
                final_score=score,
                matched_alias=None,
                metadata=source.get("metadata") or {},
            )
        )
    return hits


def search_vector_channel(
    es: Elasticsearch,
    normalized_query: NormalizedQuery,
    config: HybridSearchConfig,
) -> List[HybridHit]:
    """
    功能：执行向量检索流程，返回所有类型的语义命中集合。
    构建逻辑：
        1. 若配置禁用向量索引直接返回空列表；
        2. 构造嵌入客户端并生成查询向量；
        3. 遍历类型调用 `_search_type`，汇总 `HybridHit` 列表。
    出参入参：入参为 ES 客户端、标准化查询、混合检索配置；出参为向量通道命中集合。
    数据流：输出与字符串通道结果一起进入 `fusion`。
    模块内调用链：依赖 `_build_query_embedding`、`_search_type`。
    模块外调用链：`hybrid_searcher`。
    """
    if not config.vector_index_name or config.beta == 0.0:
        return []

    try:
        embedding_client = EmbeddingClient(config.embedding)
    except Exception as exc:  # noqa: BLE001
        logger.error("Failed to initialize embedding client: %s", exc)
        return []

    results: List[HybridHit] = []
    size = max(1, normalized_query.options.top_k)
    for type_key, candidates in normalized_query.type_candidates.items():
        embedding = _build_query_embedding(embedding_client, candidates)
        if not embedding:
            continue
        try:
            hits = _search_type(
                es,
                index_name=config.vector_index_name,
                type_key=type_key,
                embedding=embedding,
                size=size,
            )
        except Exception as exc:  # noqa: BLE001
            logger.error(
                "Vector search failed for type %s on index %s: %s",
                type_key,
                config.vector_index_name,
                exc,
            )
            continue
        results.extend(hits)
    return results


__all__ = [
    "search_vector_channel",
]
