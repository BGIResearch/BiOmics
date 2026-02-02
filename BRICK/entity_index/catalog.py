"""
功能：依据节点配置汇总实体记录，并在需要时合并 overlay 类型的别名，生成索引构建所需的实体目录。
构建逻辑：遍历配置表，调用 `loader.build_entity_records` 读取各节点的基础信息，同时收集 overlay-only 节点的同义词，以统一的 `EntityRecord` 列表形式输出。
数据流：上游为索引构建或检索模块，需要获取标准化实体集合；下游如 `pipeline.build_index_artifacts`、`index_builders` 根据返回结果构造字符串或向量索引。
调用链：例如 `pipeline` 调用 `build_catalog`，随后对每个节点生成索引文档；混合检索模块也可直接使用 catalog 进行本地匹配。
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, List, Mapping, Optional

from .config import NODE_CONFIGS
from .loader import build_entity_records, load_synonym_map
from .entity_types import EntityRecord, NodeConfig


def _collect_overlay_synonyms(configs: Mapping[str, NodeConfig]) -> Dict[str, List[str]]:
    """
    功能：收集 overlay-only 节点的别名映射，供主节点在构建 catalog 时补充使用。
    构建逻辑：遍历配置字典，筛选 `overlay_only=True` 的节点，使用 `load_synonym_map` 读取其 JSON 别名，并按实体 ID 聚合到字典中。
    出参入参：入参为节点配置映射；出参为 `entity_id -> 别名列表` 的字典。
    数据流：在 `build_catalog` 中作为前置步骤，输出结果用于 `_merge_aliases` 合并别名。
    模块内调用链：`build_catalog` -> `_collect_overlay_synonyms`。
    模块外调用链：无。
    """
    overlay_synonyms: Dict[str, List[str]] = defaultdict(list)
    for config in configs.values():
        if not config.overlay_only:
            continue
        for entity_id, aliases in load_synonym_map(config).items():
            overlay_synonyms[entity_id].extend(aliases)
    return overlay_synonyms


def _merge_aliases(existing: List[str], extras: List[str], primary_name: str, entity_id: str) -> List[str]:
    """
    功能：将主节点已有别名与 overlay 别名合并去重，避免与主名称或 ID 冲突。
    构建逻辑：维护一个去重集合，遍历两个列表的拼接结果，过滤空字符串及与主名称/ID 相同的别名，并保留首个出现顺序。
    出参入参：入参为原始别名、额外别名、主名称、实体 ID；出参为处理后的别名列表。
    数据流：在 `build_catalog` 中更新 `EntityRecord.aliases`，供后续索引与检索流程使用。
    模块内调用链：`build_catalog` -> `_merge_aliases`。
    模块外调用链：无。
    """
    merged: List[str] = []
    seen = {primary_name.lower(), entity_id.lower()}
    for alias in existing + extras:
        cleaned = alias.strip()
        if not cleaned:
            continue
        lowered = cleaned.lower()
        if lowered in seen:
            continue
        seen.add(lowered)
        merged.append(cleaned)
    return merged


def build_catalog(
    node_names: Optional[Iterable[str]] = None,
    configs: Mapping[str, NodeConfig] = NODE_CONFIGS,
) -> Dict[str, List[EntityRecord]]:
    """
    功能：加载指定节点类型的实体记录，并整合 overlay 别名，生成“节点名→实体列表”的目录。
    构建逻辑：
    1. 确定需要构建的节点集合（排除 overlay-only 类型）；
    2. 预先收集 overlay 节点别名；
    3. 遍历实体记录，根据需要附加额外别名后写入结果字典。
    出参入参：入参为可选节点名列表与配置映射；出参为 dict，value 为 `EntityRecord` 列表。
    数据流：上游为索引构建或检索初始化；下游如 `pipeline.build_index_artifacts` 读取 catalog 生成索引，混合检索可直接查询此目录。
    模块内调用链：依次调用 `_collect_overlay_synonyms`、`build_entity_records`、`_merge_aliases`。
    模块外调用链：`pipeline`, `search` 等模块使用此函数。
    """
    if node_names is None:
        selected = [
            name
            for name, cfg in configs.items()
            if not cfg.overlay_only and not cfg.is_dictionary
        ]
    else:
        selected = [name for name in node_names if not configs[name].is_dictionary]
    overlay_synonyms = _collect_overlay_synonyms(configs)

    catalog: Dict[str, List[EntityRecord]] = {}
    for node_name in selected:
        config = configs[node_name]
        if config.overlay_only or config.is_dictionary:
            # overlay-only configs have no catalog of their own
            continue
        records: List[EntityRecord] = []
        for record in build_entity_records(config):
            extra_aliases = overlay_synonyms.get(record.entity_id, [])
            if extra_aliases:
                record.aliases = _merge_aliases(
                    record.aliases,
                    extra_aliases,
                    record.primary_name,
                    record.entity_id,
                )
            records.append(record)
        catalog[node_name] = records
    return catalog


__all__ = ["build_catalog"]
