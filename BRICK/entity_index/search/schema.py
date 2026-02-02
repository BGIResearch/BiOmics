"""
功能：定义混合检索流程内部使用的数据结构，例如规范化后的查询对象、检索命中、诊断信息与日志记录。
构建逻辑：通过 dataclass/TypedDict 描述请求参数、命中结果与响应格式，使 `adapters`、`string_client`、`vector_client`、`fusion`、`hybrid_searcher` 等模块之间的数据传递清晰可信。
数据流：上游为 `adapters` 解析后的标准化查询；下游模块使用这些结构读取/写入检索命中与响应内容。
调用链：混合检索入口与各子模块会 import 本文件暴露的类型定义，用于创建和返回统一的数据对象。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

# 支持的实体类型列表，保持与设计文档及上游 payload 一致
HYBRID_TYPE_KEYS: List[str] = [
    "Gene|Protein",
    "Mutation",
    "Chemical",
    "Disease|Phenotype",
    "Process|Function|Pathway|Cell_Component",
    "Species",
    "Cell|Tissue",
]
"""
功能：列出混合检索支持的全部类型键，供输入校验与结果组织使用。
构建逻辑：与上游约定的 JSON 结构保持同步，后续如需新增类型需从此列表扩展。
数据流：`adapters` 在校验/填充候选时引用该列表。
"""


@dataclass
class QueryOptions:
    """
    功能：封装检索入口的可选参数，例如 `top_k`、诊断输出开关、调试模式等。
    构建逻辑：从上游 payload 的 `options` 字段提取或填充默认值，统一传递给后续模块。
    数据流：`adapters` 负责构建该对象，`hybrid_searcher`、`fusion` 等模块读取其中的配置。
    """

    top_k: int
    return_diagnostics: bool = False
    debug: bool = False
    type_mix_override: Optional[Dict[str, float]] = None


@dataclass
class NormalizedQuery:
    """
    功能：表示经过 `adapters` 标准化后的检索请求，包含 query_id、类型候选与统一选项。
    构建逻辑：在原始 JSON 中补齐缺失类型、清洗候选词并附带 `QueryOptions`。
    数据流：作为输入传递给字符串/向量检索通道与融合模块。
    """

    query_id: Optional[str]
    type_candidates: Dict[str, List[str]]
    options: QueryOptions


@dataclass
class ChannelScores:
    """
    功能：记录单个命中在字符串/向量通道的得分，方便融合与诊断输出。
    构建逻辑：字符串检索和向量检索分别填充对应得分，未命中默认为 0。
    数据流：`fusion` 在合并命中时使用并更新该结构；诊断响应引用这些数值。
    """

    string_score: float = 0.0
    vector_score: float = 0.0


@dataclass
class HybridHit:
    """
    功能：表示融合阶段的单个实体命中，包含基础信息、所属类型及渠道得分。
    构建逻辑：由 `fusion` 组合字符串/向量命中时生成，保留 `ChannelScores` 及匹配到的别名。
    数据流：`fusion` 输出该结构列表；`hybrid_searcher` 根据它生成标准化结果与诊断信息。
    """

    entity_id: str
    primary_name: str
    type_key: str
    node_type: str
    scores: ChannelScores
    final_score: float
    matched_alias: Optional[str] = None
    metadata: Dict[str, object] = field(default_factory=dict)


@dataclass
class DiagnosticEntry:
    """
    功能：记录命中的详细诊断信息，便于审查排序与渠道贡献。
    构建逻辑：`hybrid_searcher` 在生成响应时，将 `HybridHit` 结构转换为诊断条目。
    数据流：当 `return_diagnostics=True` 或 debug 模式开启时返回给调用方。
    """

    type_key: str
    entity_id: str
    primary_name: str
    node_type: str
    final_score: float
    channel_scores: ChannelScores
    matched_alias: Optional[str] = None
    llm_adjusted: bool = False
    extra: Dict[str, object] = field(default_factory=dict)


@dataclass
class HybridResponse:
    """
    功能：混合检索最终返回的结构，包含标准化结果、可选诊断信息与日志。
    构建逻辑：`hybrid_searcher` 根据融合结果组成 standardized 字典，并在需要时附加 diagnostics 与 logs。
    数据流：作为对外响应返回给上游服务或调用方。
    """

    query_id: Optional[str]
    standardized: Dict[str, List[str]]
    diagnostics: Optional[List[DiagnosticEntry]] = None
    logs: Optional[Dict[str, object]] = None


__all__ = [
    "HYBRID_TYPE_KEYS",
    "QueryOptions",
    "NormalizedQuery",
    "ChannelScores",
    "HybridHit",
    "DiagnosticEntry",
    "HybridResponse",
]
