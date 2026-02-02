"""
混合检索配置入口。

- 约定从 `.env` 中读取 HYBRID_* 变量，如未设置则回退至设计文档中的默认值。
- 将嵌入模型配置与 Elasticsearch 索引名集中在此，让后续模块统一引用。
- 优先复用 `_settings.py` 已有的嵌入/ES 配置，必要时在此覆盖。
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Dict, Optional

from ..embedding import EmbeddingConfig
from ..es_writer import ESConfig

# 默认类型权重配置
DEFAULT_TYPE_MIX: Dict[str, float] = {
    "Gene|Protein": 0.25,
    "Disease|Phenotype": 0.25,
    "Process|Function|Pathway|Cell_Component": 0.20,
    "Chemical": 0.10,
    "Species": 0.10,
    "Cell|Tissue": 0.10,
    "Mutation": 0.0,
}
"""
功能：定义混合检索类型权重的默认配置，作为融合阶段 `type_mix` 的回退值。
构建逻辑：根据设计文档中约定的比例设定权重字典。
数据流：`fusion` 或入口模块在未提供自定义权重时使用该常量。
"""

DEFAULT_TOP_K = 10
"""
功能：定义整合后返回结果的默认 TopK 数量。
构建逻辑：当上游 payload 未指定 `options.top_k` 时使用该值。
数据流：`adapters` 在标准化请求时读取该常量。
"""


@dataclass
class HybridSearchConfig:
    """
    功能：封装混合检索运行所需的核心配置，包括 ES、嵌入与检索策略参数。
    构建逻辑：组合 ES 配置、嵌入配置、索引名称、默认权重等信息，供入口模块统一持有。
    数据流：`hybrid_searcher` 初始化时接收该对象，后续各子模块通过属性读取所需参数。
    """

    es: ESConfig
    embedding: EmbeddingConfig
    string_index_name: str
    vector_index_name: Optional[str]
    alias_index_name: str
    type_mix: Dict[str, float]
    top_k: int = DEFAULT_TOP_K
    alpha: float = 0.6
    beta: float = 0.4


def get_default_index_name(suffix: str, fallback: str) -> str:
    """
    功能：根据 `HYBRID_ES_INDEX` 环境变量或 fallback，生成字符串/向量索引名。
    构建逻辑：读取 `HYBRID_ES_INDEX`，若设置则与 suffix 拼接；否则返回 fallback 名称。
    出参入参：入参为 suffix（如 "string"）和 fallback；出参为索引名字符串。
    数据流：`get_search_config` 在装配配置时调用。
    模块内调用链：仅在本模块内部使用。
    模块外调用链：无。
    """
    base = os.environ.get("HYBRID_ES_INDEX")
    if base:
        return f"{base}_{suffix}" if suffix else base
    return fallback


def _sanitize_prefix(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    sanitized = value.strip()
    if not sanitized:
        return None
    sanitized = sanitized.rstrip("*")
    for suffix in ("_string", "_vector"):
        if sanitized.endswith(suffix):
            sanitized = sanitized[: -len(suffix)]
    sanitized = sanitized.rstrip("_-")
    return sanitized or None


def _resolve_alias_index(default_prefix: str = "entity") -> str:
    explicit = os.environ.get("HYBRID_ALIAS_INDEX")
    if explicit:
        return explicit
    for candidate in (
        os.environ.get("ES_INDEX_PREFIX"),
        os.environ.get("HYBRID_ES_INDEX"),
    ):
        prefix = _sanitize_prefix(candidate)
        if prefix:
            return f"{prefix}_all_name2id_string"
    return f"{default_prefix}_all_name2id_string"


def _load_type_mix() -> Dict[str, float]:
    """
    功能：从环境变量解析混合检索类型权重，若未设置则返回默认值。
    构建逻辑：读取 `HYBRID_TYPE_MIX`，期待 JSON 字典格式；解析失败时抛出异常。
    出参入参：无入参；出参为类型权重字典。
    数据流：`get_search_config` 调用以获取权重配置。
    模块内调用链：仅在本模块内使用。
    模块外调用链：无。
    """
    raw = os.environ.get("HYBRID_TYPE_MIX")
    if not raw:
        return DEFAULT_TYPE_MIX
    raw = raw.strip()
    candidates = [raw]
    if raw and raw[0] in {"'", '"'} and raw[-1] == raw[0]:
        candidates.append(raw[1:-1].strip())
    parsed = None
    for candidate in candidates:
        if not candidate:
            continue
        obj: Optional[Dict[str, float]] = None
        try:
            obj = json.loads(candidate)
        except json.JSONDecodeError:
            from ast import literal_eval

            try:
                obj = literal_eval(candidate)
            except (ValueError, SyntaxError):
                continue
        if isinstance(obj, dict):
            parsed = obj
            break
    if parsed is None:
        raise ValueError("HYBRID_TYPE_MIX 必须是 JSON 字典字符串")
    if not isinstance(parsed, dict):
        raise ValueError("HYBRID_TYPE_MIX 必须是 JSON 字典字符串")
    return {k: float(v) for k, v in parsed.items()}


def get_search_config() -> HybridSearchConfig:
    """
    功能：从环境变量与默认配置组装混合检索所需的 `HybridSearchConfig`。
    构建逻辑：
        1. 通过 `ESConfig.from_env()`、`EmbeddingConfig.from_env()` 获取基础连接配置；
        2. 解析类型权重、索引名称、alpha/beta/top_k 等参数；
        3. 若禁用向量检索，自动清空向量索引名称并将 beta 置 0；
        4. 返回封装后的配置对象。
    出参入参：无入参；出参为 `HybridSearchConfig`。
    数据流：`hybrid_searcher` 初始化或 DI 容器装配时调用；下游模块通过 config 读取参数。
    模块内调用链：依赖 `_load_type_mix`、`get_default_index_name`。
    模块外调用链：混合检索入口调用该函数以获取配置。
    """
    es_cfg = ESConfig.from_env()
    embedding_cfg = EmbeddingConfig.from_env()
    type_mix = _load_type_mix()

    string_index = get_default_index_name("string", "entity_hybrid_string")
    vector_index = get_default_index_name("vector", "entity_hybrid_vector")
    alias_index = _resolve_alias_index()

    alpha = float(os.environ.get("HYBRID_ALPHA", "0.6"))
    beta = float(os.environ.get("HYBRID_BETA", "0.4"))
    top_k = int(os.environ.get("HYBRID_TOP_K", str(DEFAULT_TOP_K)))

    if os.environ.get("HYBRID_DISABLE_VECTOR", "false").lower() == "true":
        vector_index = None
        beta = 0.0

    return HybridSearchConfig(
        es=es_cfg,
        embedding=embedding_cfg,
        string_index_name=string_index,
        vector_index_name=vector_index,
        alias_index_name=alias_index,
        type_mix=type_mix,
        top_k=top_k,
        alpha=alpha,
        beta=beta,
    )


__all__ = [
    "HybridSearchConfig",
    "DEFAULT_TYPE_MIX",
    "DEFAULT_TOP_K",
    "get_search_config",
]
