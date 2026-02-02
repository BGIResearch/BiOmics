"""
功能：封装与 Elasticsearch 字符串检索相关的查询模板与调用流程，支持精确、模糊、拼音等多种匹配方式，为混合检索提供字符串候选。
构建逻辑：组合 `adapters` 标准化后的查询词与 `text_utils` 输出的拼音 token 构造 ES 查询；对返回结果做归一化处理并转换为内部命中结构。
数据流：上游为 `adapters` 输出的 `NormalizedQuery`；下游是 `fusion` 模块与 `hybrid_searcher`，共同处理字符串候选与向量候选。
调用链：`HybridEntitySearcher` 将在字符串通道阶段调用本模块函数，然后再与向量检索结果合并。
"""

from __future__ import annotations

import logging
from typing import Dict, Iterable, List, Sequence

from elasticsearch import Elasticsearch  # type: ignore

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

from .schema import ChannelScores, HybridHit, NormalizedQuery
from ..text_utils import to_pinyin_tokens


def _build_query_tokens(candidates: Iterable[str]) -> Dict[str, List[str]]:
    """
    功能：根据候选词生成检索所需的查询 token，包括原词、拼音等。
    构建逻辑：遍历候选词，收集原词列表，并调用 `to_pinyin_tokens` 生成拼音 token 去重。
    出参入参：入参为候选词迭代器；出参为包含 `terms` 与 `pinyin` 的字典。
    数据流：`search_string_channel` 在发起 ES 查询前调用本函数准备查询参数。
    模块内调用链：当前模块内部使用。
    模块外调用链：无。
    """
    terms: List[str] = []
    pinyin_tokens: List[str] = []
    seen_terms = set()
    seen_pinyin = set()
    for candidate in candidates:
        if candidate is None:
            continue
        value = candidate.strip()
        if not value:
            continue
        marker = value.lower()
        if marker not in seen_terms:
            seen_terms.add(marker)
            terms.append(value)
        for token in to_pinyin_tokens(value):
            lower_token = token.lower()
            if lower_token in seen_pinyin:
                continue
            seen_pinyin.add(lower_token)
            pinyin_tokens.append(token)
    return {"terms": terms, "pinyin": pinyin_tokens}


def _search_type(  # noqa: PLR0913
    es: Elasticsearch,
    index_name: str,
    type_key: str,
    candidates: Iterable[str],
    size: int,
) -> List[HybridHit]:
    """
    功能：在指定索引上执行字符串检索，返回同一类型的命中列表。
    构建逻辑：
        1. 调用 `_build_query_tokens` 构造多字段查询；
        2. 向 ES 发起检索并按得分排序；
        3. 将命中转换为 `HybridHit`，仅填充字符串得分。
    出参入参：入参为 ES 客户端、索引名、类型键、候选词、返回数量；出参为 `HybridHit` 列表。
    数据流：`search_string_channel` 对每个类型调用本函数并合并结果。
    模块内调用链：依赖 `_build_query_tokens`。
    模块外调用链：无。
    """
    tokens = _build_query_tokens(candidates)
    if not tokens["terms"] and not tokens["pinyin"]:
        return []

    should_clauses: List[Dict[str, object]] = []
    for term in tokens["terms"]:
        term_lower = term.lower()
        should_clauses.extend(
            [
                # 精确短语匹配（权重最高，解决 "T cell" vs "T cell domain"）
                {"match_phrase": {"primary_name": {"query": term, "boost": 20.0}}},
                {"match_phrase": {"aliases": {"query": term, "boost": 15.0}}},
                # term 精确匹配（大小写不敏感）
                {"term": {"primary_name": {"value": term, "boost": 10.0}}},
                {"term": {"primary_name": {"value": term_lower, "boost": 10.0}}},
                {"term": {"primary_name": {"value": term.upper(), "boost": 10.0}}},
                {"term": {"aliases": {"value": term, "boost": 8.0}}},
                {"term": {"aliases": {"value": term_lower, "boost": 8.0}}},
                {"term": {"aliases": {"value": term.upper(), "boost": 8.0}}},
                # 分词匹配
                {"match": {"search_terms": {"query": term, "boost": 3.0}}},
                {"match": {"search_terms.ngram": {"query": term, "boost": 1.5}}},
                {"match": {"search_terms.prefix": {"query": term, "boost": 1.0}}},
            ]
        )
    for token in tokens["pinyin"]:
        should_clauses.append(
            {"term": {"pinyin_terms": {"value": token, "boost": 1.0}}}
        )

    if not should_clauses:
        return []

    filters: List[Dict[str, object]] = []
    ctype_filters = TYPE_KEY_TO_CTYPES.get(type_key)
    if ctype_filters:
        filters.append({"terms": {"ctype": list(ctype_filters)}})

    body: Dict[str, object] = {
        "size": size,
        "track_total_hits": False,
        "query": {
            "bool": {
                "should": should_clauses,
                "minimum_should_match": 1,
                "filter": filters,
            }
        },
    }

    try:
        response = es.search(index=index_name, body=body)
    except Exception as exc:  # noqa: BLE001
        logger.error(
            "String search failed for type %s on index %s: %s",
            type_key,
            index_name,
            exc,
        )
        return []

    hits: List[HybridHit] = []
    candidate_markers = {term.lower(): term for term in tokens["terms"]}
    for item in response.get("hits", {}).get("hits", []):  # type: ignore[assignment]
        source = item.get("_source", {}) or {}
        entity_id = source.get("entity_id")
        primary_name = source.get("primary_name")
        if not entity_id or not primary_name:
            continue
        aliases = source.get("aliases") or []
        matched_alias = None
        for alias in [primary_name, *aliases]:
            if alias and alias.lower() in candidate_markers:
                matched_alias = alias
                break
        score = float(item.get("_score") or 0.0)
        hits.append(
            HybridHit(
                entity_id=entity_id,
                primary_name=primary_name,
                type_key=type_key,
                node_type=source.get("ctype", type_key),
                scores=ChannelScores(string_score=score, vector_score=0.0),
                final_score=score,
                matched_alias=matched_alias,
                metadata=source.get("metadata") or {},
            )
        )
    return hits


def search_string_channel(
    es: Elasticsearch,
    normalized_query: NormalizedQuery,
    index_name: str,
) -> List[HybridHit]:
    """
    功能：对所有类型执行字符串检索并汇总命中结果。
    构建逻辑：
        1. 根据 `NormalizedQuery` 遍历类型候选；
        2. 调用 `_search_type` 获取单类型命中；
        3. 合并所有类型的 `HybridHit` 列表并返回。
    出参入参：入参为 ES 客户端、标准化查询、字符串索引名；出参为字符串通道命中集合。
    数据流：字符串检索结果将与向量通道结果在 `fusion` 中合并。
    模块内调用链：依赖 `_search_type`。
    模块外调用链：`hybrid_searcher`。
    """
    results: List[HybridHit] = []
    size = max(1, normalized_query.options.top_k)
    for type_key, candidates in normalized_query.type_candidates.items():
        size_per_type = max(size, len(candidates))
        hits = _search_type(
            es,
            index_name=index_name,
            type_key=type_key,
            candidates=candidates,
            size=size_per_type,
        )
        results.extend(hits)
    return results


__all__ = [
    "search_string_channel",
]
