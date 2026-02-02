"""
功能：整合字符串与向量检索命中结果，按照类型权重与通道权重计算最终得分，并输出 TopK 候选，供混合检索模块返回。
构建逻辑：提供分数归一化、去重、权重融合等函数，读取 `settings` 中的默认参数（α/β、type_mix），同时保留命中渠道信息以便诊断输出。
数据流：上游为 `string_client` 与 `vector_client` 提供的候选列表；下游是 `hybrid_searcher` 或其他调用者，将融合后的结果转化为标准化输出。
调用链：`HybridEntitySearcher` 在两个通道检索完成后调用本模块，生成最终排序结果，并交给结果标准化模块。
"""

from __future__ import annotations

from typing import Dict, List

from .schema import ChannelScores, HybridHit, NormalizedQuery


def _merge_channel_scores(existing: ChannelScores, incoming: ChannelScores) -> ChannelScores:
    """
    功能：合并同一实体的渠道得分，字符串/向量分值取最大或叠加策略（待实现）。
    构建逻辑：比较已有得分与新得分，按规则合并后返回新的 `ChannelScores`。
    出参入参：入参为已有得分与新得分；出参为合并后的得分对象。
    数据流：`merge_hits` 在聚合重复实体时调用。
    模块内调用链：本模块内部使用。
    模块外调用链：无。
    """
    return ChannelScores(
        string_score=max(existing.string_score, incoming.string_score),
        vector_score=max(existing.vector_score, incoming.vector_score),
    )


def merge_hits(
    string_hits: List[HybridHit],
    vector_hits: List[HybridHit],
) -> Dict[str, HybridHit]:
    """
    功能：按实体 ID 聚合字符串与向量命中，保留渠道得分与最高得分实体信息。
    构建逻辑：以实体 ID 为 key 构建字典，遇到重复时调用 `_merge_channel_scores` 更新得分，并根据得分策略更新 `HybridHit`。
    出参入参：入参为字符串通道与向量通道命中列表；出参为实体 ID -> `HybridHit` 的字典。
    数据流：`apply_fusion` 调用本函数作为融合前置步骤。
    模块内调用链：依赖 `_merge_channel_scores`。
    模块外调用链：无。
    """
    merged: Dict[str, HybridHit] = {}

    def _clone_hit(hit: HybridHit) -> HybridHit:
        return HybridHit(
            entity_id=hit.entity_id,
            primary_name=hit.primary_name,
            type_key=hit.type_key,
            node_type=hit.node_type,
            scores=ChannelScores(
                string_score=hit.scores.string_score,
                vector_score=hit.scores.vector_score,
            ),
            final_score=hit.final_score,
            matched_alias=hit.matched_alias,
            metadata=dict(hit.metadata),
        )

    for hit in string_hits:
        merged[hit.entity_id] = _clone_hit(hit)

    for hit in vector_hits:
        existing = merged.get(hit.entity_id)
        if existing is None:
            merged[hit.entity_id] = _clone_hit(hit)
            continue
        existing.scores = _merge_channel_scores(existing.scores, hit.scores)
        existing.final_score = max(existing.final_score, hit.final_score)
        if not existing.matched_alias and hit.matched_alias:
            existing.matched_alias = hit.matched_alias
        if hit.metadata:
            existing.metadata.update(hit.metadata)
        if not existing.node_type:
            existing.node_type = hit.node_type
        if not existing.primary_name:
            existing.primary_name = hit.primary_name
    return merged


def apply_type_weights(
    merged_hits: Dict[str, HybridHit],
    type_mix: Dict[str, float],
    alpha: float,
    beta: float,
) -> List[HybridHit]:
    """
    功能：根据类型权重与通道权重计算最终得分，并输出排序后的命中列表。
    构建逻辑：
        1. 遍历合并后的命中，按 `alpha`/`beta` 组合字符串与向量得分；
        2. 乘以 `type_mix` 中对应类型的权重；
        3. 更新 `final_score` 并按分值排序。
    出参入参：入参为合并命字典、类型权重、通道权重；出参为排序后的命中列表。
    数据流：`apply_fusion` 使用该函数得到最终排名。
    模块内调用链：无。
    模块外调用链：`hybrid_searcher` 直接消费结果。
    """
    ranked: List[HybridHit] = []
    for hit in merged_hits.values():
        weight = type_mix.get(hit.type_key, 1.0)
        combined = alpha * hit.scores.string_score + beta * hit.scores.vector_score
        hit.final_score = combined * weight
        ranked.append(hit)
    ranked.sort(key=lambda item: item.final_score, reverse=True)
    return ranked


def apply_fusion(
    normalized_query: NormalizedQuery,
    string_hits: List[HybridHit],
    vector_hits: List[HybridHit],
    type_mix: Dict[str, float],
    alpha: float,
    beta: float,
) -> List[HybridHit]:
    """
    功能：执行完整的融合流程，输出截断后的 TopK 命中列表。
    构建逻辑：
        1. 调用 `merge_hits` 聚合字符串/向量命中；
        2. 调用 `apply_type_weights` 计算最终得分并排序；
        3. 根据 `normalized_query.options.top_k` 截断列表。
    出参入参：入参为标准化查询、两通道命中、类型权重、通道权重；出参为融合后的 `HybridHit` 列表。
    数据流：融合结果直接传递给 `hybrid_searcher` 进行响应构建。
    模块内调用链：依赖 `merge_hits`、`apply_type_weights`。
    模块外调用链：`hybrid_searcher`。
    """
    merged = merge_hits(string_hits, vector_hits)
    weighted_hits = apply_type_weights(merged, type_mix, alpha, beta)
    top_k = max(1, normalized_query.options.top_k)
    return weighted_hits[:top_k]


__all__ = [
    "apply_fusion",
]
