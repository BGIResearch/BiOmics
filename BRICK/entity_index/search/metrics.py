"""
功能：提供混合检索相关的统计与评估工具（例如语义索引数量对比、命中率统计），方便与既有基线进行对比。
构建逻辑：计划封装读取 ES 索引计数或本地 catalog 的统计函数，输出结构化报表，供 mentor 验证索引构建成果。
数据流：上游接入混合检索或索引构建完成后的数据；下游可将统计结果写入日志或导出给分析脚本。
调用链：在构建完索引后，部署脚本或测试脚本调用本模块函数，生成统计报表与对比数据。
"""

# TODO: 后续实现具体的统计函数



from __future__ import annotations

from typing import Dict, Iterable, List, Optional

from elasticsearch import Elasticsearch  # type: ignore

from ..entity_types import EntityRecord
from ..pipeline import IndexArtifacts


def summarize_catalog_counts(catalog: Dict[str, List[EntityRecord]]) -> Dict[str, int]:
    """
    功能：统计本地 catalog 中各节点类型的实体数量，供与 ES 索引对比。
    构建逻辑：遍历 catalog 的节点→实体列表映射，计算列表长度并形成字典。
    出参入参：入参为 catalog（节点名 -> 实体列表）；出参为节点名 -> 数量的字典。
    数据流：索引构建完成后调用，输出作为离线基线。
    模块内调用链：其他统计函数可复用该结果。
    模块外调用链：索引验证脚本、Notebook。
    """
    counts: Dict[str, int] = {}
    for node, records in catalog.items():
        counts[node] = len(records)
    return counts


def fetch_es_counts(es: Elasticsearch, index_names: Iterable[str]) -> Dict[str, int]:
    """
    功能：查询 Elasticsearch 中指定索引的文档数量，返回索引 -> 数量的映射。
    构建逻辑：遍历索引名称，调用 ES `count` API，取返回的 `count` 值。
    出参入参：入参为 ES 客户端与索引名称序列；出参为索引名 -> 文档数的字典。
    数据流：与 catalog 统计对比，检查写入是否完整。
    模块内调用链：可结合 `compare_counts` 输出报表。
    模块外调用链：部署验证脚本。
    """
    counts: Dict[str, int] = {}
    for name in index_names:
        counts[name] = es.count(index=name)["count"]  # type: ignore[index]
    return counts


def compare_counts(
    catalog_counts: Dict[str, int],
    string_index_counts: Dict[str, int],
    vector_index_counts: Optional[Dict[str, int]] = None,
) -> List[Dict[str, object]]:
    """
    功能：比较 catalog 数量与 ES 字符串/向量索引数量，生成差异报表。
    构建逻辑：遍历 catalog 统计，读取对应字符串/向量索引 count，计算差值得到汇总行。
    出参入参：入参为 catalog 统计、字符串索引统计、可选向量索引统计；出参为字典列表，含节点、各类数量与差值。
    数据流：在日志或 Notebook 中展示，辅助 QA。
    模块内调用链：可结合 `summarize_catalog_counts`、`fetch_es_counts` 使用。
    模块外调用链：索引验收流程。
    """
    report: List[Dict[str, object]] = []
    for node, cat_count in catalog_counts.items():
        string_count = string_index_counts.get(f"entity_{node.lower()}_string", 0)
        vector_count = (
            vector_index_counts.get(f"entity_{node.lower()}_vector", 0)
            if vector_index_counts
            else None
        )
        entry: Dict[str, object] = {
            "node": node,
            "catalog_count": cat_count,
            "string_count": string_count,
            "string_delta": string_count - cat_count,
        }
        if vector_index_counts is not None:
            entry["vector_count"] = vector_count
            entry["vector_delta"] = (vector_count or 0) - cat_count
        report.append(entry)
    return report


def summarize_artifacts_counts(artifacts: Iterable[IndexArtifacts]) -> Dict[str, Dict[str, int]]:
    """
    功能：从 `IndexArtifacts` 列表统计各节点的字符串/向量文档数量，用于写入前自检。
    构建逻辑：遍历工件，计算字符串文档数与向量文档数，组织成嵌套字典。
    出参入参：入参为工件迭代器；出参为节点 -> {"string": 数量, "vector": 数量}。
    数据流：写入 ES 之前或之后调用，快速检查写入量是否符合预期。
    模块内调用链：可与 `compare_counts` 搭配使用。
    模块外调用链：测试脚本、监控任务。
    """
    stats: Dict[str, Dict[str, int]] = {}
    for artifact in artifacts:
        stats[artifact.node] = {
            "string": len(artifact.string_documents),
            "vector": len(artifact.vector_documents or []),
        }
    return stats
