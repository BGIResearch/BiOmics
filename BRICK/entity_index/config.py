"""
功能：集中维护实体节点的基础配置，包括词表路径、是否生成语义向量、别名来源等，作为索引构建的单一事实来源。
构建逻辑：定义 `NODE_CONFIGS` 映射，将每种节点类型与其数据文件、字段映射、额外元信息绑定，供管线统一读取。
数据流：上游通常为需要查询节点配置的索引构建流程（`pipeline.build_index_artifacts`）或检索模块；下游则是词表加载器（`loader`）和 catalog 组装器（`catalog`），按照配置定位数据文件。
调用链：索引构建或检索逻辑首先读取 `NODE_CONFIGS`，再据此调用 `loader.build_entity_records` 等函数完成后续处理。
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

from .entity_types import NodeConfig

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "entity_vocab"
TABLE_DIR = DATA_DIR / "nodes_table"
NAME_DIR = DATA_DIR / "nodes_name_json"

DEFAULT_COLUMNS = {
    "id": "id",
    "name": "name",
    "type": "type",
    "definition": "def",
    "synonyms": "synonym",
}

NODE_CONFIGS: Dict[str, NodeConfig] = {
    "GENE": NodeConfig(
        name="GENE",
        table_path=TABLE_DIR / "GENE.node.csv",
        enable_vector=False,
        columns=DEFAULT_COLUMNS,
        synonyms_path=NAME_DIR / "Gene.json",
        id_prefix="NCBI",
    ),
    "GENE_mm": NodeConfig(
        name="GENE_mm",
        table_path=TABLE_DIR / "GENE_mm.node.csv",
        enable_vector=False,
        columns=DEFAULT_COLUMNS,
        synonyms_path=None,
        id_prefix="NCBI",
        extra_metadata={"species": "Mus musculus"},
    ),
    "GO": NodeConfig(
        name="GO",
        table_path=TABLE_DIR / "GO-PLUS.node.csv",
        enable_vector=True,
        columns=DEFAULT_COLUMNS,
        synonyms_path=NAME_DIR / "GO.json",
        id_prefix="GO",
    ),
    "CL": NodeConfig(
        name="CL",
        table_path=TABLE_DIR / "CL.node.csv",
        enable_vector=True,
        columns=DEFAULT_COLUMNS,
        synonyms_path=NAME_DIR / "CL.json",
        id_prefix="CL",
    ),
    "DOID": NodeConfig(
        name="DOID",
        table_path=TABLE_DIR / "DOID.node.csv",
        enable_vector=True,
        columns=DEFAULT_COLUMNS,
        synonyms_path=NAME_DIR / "DOID.json",
        id_prefix="DOID",
    ),
    "HPO": NodeConfig(
        name="HPO",
        table_path=TABLE_DIR / "HPO.node.csv",
        enable_vector=True,
        columns=DEFAULT_COLUMNS,
        synonyms_path=NAME_DIR / "HPO.json",
        id_prefix="HP",
    ),
    "HSAPDV": NodeConfig(
        name="HSAPDV",
        table_path=TABLE_DIR / "HSAPDV.node.csv",
        enable_vector=False,
        columns=DEFAULT_COLUMNS,
        synonyms_path=None,
        id_prefix="HSAPDV",
    ),
    "KEGG": NodeConfig(
        name="KEGG",
        table_path=TABLE_DIR / "KEGG.node.csv",
        enable_vector=False,
        columns=DEFAULT_COLUMNS,
        synonyms_path=NAME_DIR / "KEGG.json",
        id_prefix="KEGG",
    ),
    "MESH": NodeConfig(
        name="MESH",
        table_path=TABLE_DIR / "MESH.node.csv",
        enable_vector=True,
        columns=DEFAULT_COLUMNS,
        synonyms_path=NAME_DIR / "MESH.json",
        id_prefix="MESH",
    ),
    "MMUSDV": NodeConfig(
        name="MMUSDV",
        table_path=TABLE_DIR / "MMUSDV.node.csv",
        enable_vector=False,
        columns=DEFAULT_COLUMNS,
        synonyms_path=None,
        id_prefix="MMUSDV",
    ),
    "NCBITAXON": NodeConfig(
        name="NCBITAXON",
        table_path=TABLE_DIR / "NCBITAXON.node.csv",
        enable_vector=False,
        columns=DEFAULT_COLUMNS,
        synonyms_path=None,
        id_prefix="NCBITAXON",
    ),
    "Stage": NodeConfig(
        name="Stage",
        table_path=TABLE_DIR / "HSAPDV.node.csv",
        enable_vector=False,
        columns=DEFAULT_COLUMNS,
        synonyms_path=NAME_DIR / "Stage.json",
        id_prefix="HSAPDV",
    ),
    "UBERON": NodeConfig(
        name="UBERON",
        table_path=TABLE_DIR / "UBERON.node.csv",
        enable_vector=True,
        columns=DEFAULT_COLUMNS,
        synonyms_path=NAME_DIR / "UBERON.json",
        id_prefix="UBERON",
    ),
    "ALL_NAME2ID": NodeConfig(
        name="ALL_NAME2ID",
        table_path=None,
        enable_vector=False,
        columns=DEFAULT_COLUMNS,
        synonyms_path=None,
        dictionary_path=NAME_DIR / "all_name2id.json",
        is_dictionary=True,
    ),
}
"""Per-node configuration matching the decisions captured in target_plan.md."""
