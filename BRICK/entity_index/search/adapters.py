"""
功能：解析与校验上游传入的检索 JSON 负载，统一填充默认选项，并将各类型候选整理为内部结构，供后续字符串/向量检索使用。
构建逻辑：在正式实现中将提供数据类或函数，用于读取 `settings.py` 的默认配置、校验七大类型字段是否存在、剔除非法值，并生成标准化的查询对象。
数据流：上游来自 BiomicsAgent 或其他抽取模块生成的半结构化 JSON；下游为 `string_client`、`vector_client` 等检索通道以及最终的 `hybrid_searcher`。
调用链：混合检索入口（`HybridEntitySearcher`）会首先调用本模块的适配函数获取规范化查询，再依次进入各检索通道。
"""

from __future__ import annotations

import json
from typing import Dict, Iterable, List, Mapping, MutableMapping

from .schema import HYBRID_TYPE_KEYS, NormalizedQuery, QueryOptions
from .settings import DEFAULT_TOP_K, HybridSearchConfig


def _sanitize_candidates(candidates: Iterable[str]) -> List[str]:
    """
    功能：清洗候选字符串，去除空值、首尾空格与重复项。
    构建逻辑：遍历输入序列，对非空字符串进行 strip 并使用集合去重，保持原始顺序。
    出参入参：入参为候选字符串迭代器；出参为清洗后的列表。
    数据流：`normalize_payload` 在处理各类型候选时调用本函数。
    模块内调用链：仅在当前模块中使用。
    模块外调用链：无。
    """
    seen = set()
    cleaned: List[str] = []
    for candidate in candidates:
        if candidate is None:
            continue
        if not isinstance(candidate, str):
            candidate = str(candidate)
        value = candidate.strip()
        if not value:
            continue
        marker = value.lower()
        if marker in seen:
            continue
        seen.add(marker)
        cleaned.append(value)
    return cleaned


def _prepare_type_candidates(raw_payload: Mapping[str, object]) -> Dict[str, List[str]]:
    """
    功能：从原始 JSON 中提取各类型候选列表，确保缺失类型补空并执行清洗。
    构建逻辑：按 `HYBRID_TYPE_KEYS` 迭代读取列表值，若类型不存在则赋空列表，再调用 `_sanitize_candidates` 处理。
    出参入参：入参为原始 payload；出参为字典，key 为类型、value 为清洗后的候选列表。
    数据流：`normalize_payload` 组合标准化查询时调用。
    模块内调用链：依赖 `_sanitize_candidates`。
    模块外调用链：无。
    """
    prepared: Dict[str, List[str]] = {}
    for type_key in HYBRID_TYPE_KEYS:
        raw_value = raw_payload.get(type_key, [])
        if isinstance(raw_value, str):
            candidates: Iterable[str] = [raw_value]
        elif isinstance(raw_value, Iterable):
            candidates = raw_value  # type: ignore[assignment]
        else:
            candidates = []
        prepared[type_key] = _sanitize_candidates(candidates)
    return prepared


def _prepare_options(raw_options: Mapping[str, object], config: HybridSearchConfig) -> QueryOptions:
    """
    功能：解析 payload 中的 `options` 字段，填充默认的 `top_k`、诊断和调试开关，并支持类型权重覆盖。
    构建逻辑：读取 `top_k`、`return_diagnostics`、`debug` 等键，若缺失则采用配置中的默认值，
              同时允许传入自定义 `type_mix` 用于覆盖默认权重。
    出参入参：入参为原始 options 与混合检索配置；出参为 `QueryOptions`。
    数据流：`normalize_payload` 在构建 `NormalizedQuery` 时调用。
    模块内调用链：无。
    模块外调用链：无。
    """
    options_dict: MutableMapping[str, object]
    if isinstance(raw_options, MutableMapping):
        options_dict = raw_options
    else:
        options_dict = {}

    def _to_bool(value: object, default: bool) -> bool:
        if value is None:
            return default
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, float)):
            return bool(value)
        if isinstance(value, str):
            lowered = value.strip().lower()
            if lowered in {"true", "1", "yes", "y"}:
                return True
            if lowered in {"false", "0", "no", "n", ""}:
                return False
        return default

    raw_top_k = options_dict.get("top_k", config.top_k or DEFAULT_TOP_K)
    try:
        top_k = int(raw_top_k)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        top_k = config.top_k or DEFAULT_TOP_K
    if top_k <= 0:
        top_k = config.top_k or DEFAULT_TOP_K

    return_diagnostics = _to_bool(options_dict.get("return_diagnostics"), False)
    debug = _to_bool(options_dict.get("debug"), False)

    override_raw = options_dict.get("type_mix_override", options_dict.get("type_mix"))
    override_dict: Dict[str, float] | None = None
    if override_raw:
        parsed: Dict[str, object] | None = None
        if isinstance(override_raw, Mapping):
            parsed = dict(override_raw.items())
        elif isinstance(override_raw, str):
            try:
                parsed_obj = json.loads(override_raw)
            except json.JSONDecodeError:
                parsed = None
            else:
                if isinstance(parsed_obj, Mapping):
                    parsed = dict(parsed_obj.items())
        if parsed:
            override_dict = {}
            for key, value in parsed.items():
                if key not in HYBRID_TYPE_KEYS:
                    continue
                try:
                    override_dict[key] = float(value)  # type: ignore[arg-type]
                except (TypeError, ValueError):
                    continue
            if not override_dict:
                override_dict = None

    return QueryOptions(
        top_k=top_k,
        return_diagnostics=return_diagnostics,
        debug=debug,
        type_mix_override=override_dict,
    )


def normalize_payload(payload: Mapping[str, object], config: HybridSearchConfig) -> NormalizedQuery:
    """
    功能：校验并标准化上游 JSON 请求，输出 `NormalizedQuery` 对象。
    构建逻辑：
        1. 读取 `query_id` 与 `options` 字段；
        2. 调用 `_prepare_type_candidates` 获得清洗后的候选；
        3. 调用 `_prepare_options` 生成统一配置；
        4. 返回封装后的 `NormalizedQuery`。
    出参入参：入参为原始 payload 与混合检索配置；出参为 `NormalizedQuery`。
    数据流：`HybridEntitySearcher` 在执行检索前调用本函数。
    模块内调用链：依赖 `_prepare_type_candidates`、`_prepare_options`。
    模块外调用链：`hybrid_searcher`、测试脚本。
    """
    query_id = payload.get("query_id")
    if query_id is not None and not isinstance(query_id, str):
        query_id = str(query_id)

    raw_options: Mapping[str, object]
    options_payload = payload.get("options")
    if isinstance(options_payload, Mapping):
        raw_options = options_payload
    else:
        raw_options = {}

    type_candidates = _prepare_type_candidates(payload)
    options = _prepare_options(raw_options, config)
    return NormalizedQuery(
        query_id=query_id,
        type_candidates=type_candidates,
        options=options,
    )


__all__ = [
    "normalize_payload",
]
