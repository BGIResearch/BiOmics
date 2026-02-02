"""
功能：定义实体索引流程中的核心数据结构（节点配置与实体记录），为后续加载、索引构建和检索提供统一格式。
构建逻辑：使用 dataclass 描述节点配置、实体记录及附属元信息，确保在传递过程中字段含义清晰。
数据流：上游模块（`config`、`loader`）创建这些数据对象，下游模块（`catalog`、`pipeline`、检索模块）读取并加工。
调用链：例如 `loader.build_entity_records` 返回 `EntityRecord`，`pipeline.build_index_artifacts` 接收并生成索引文档。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional


@dataclass(frozen=True)
class NodeConfig:
    """Configuration for a single ontology/node type index."""

    name: str  # human-friendly label, e.g. "GO"
    table_path: Optional[Path]
    enable_vector: bool
    columns: Dict[str, str]
    synonyms_path: Optional[Path] = None
    id_prefix: Optional[str] = None
    description_fields: Optional[List[str]] = None
    extra_metadata: Dict[str, str] = field(default_factory=dict)
    overlay_only: bool = False
    dictionary_path: Optional[Path] = None
    is_dictionary: bool = False

    def resolve_table(self) -> Optional[Path]:
        if self.table_path is None:
            return None
        path = self.table_path
        if not path.exists():
            raise FileNotFoundError(f"Table file not found for {self.name}: {path}")
        return path

    def resolve_dictionary(self) -> Optional[Path]:
        if self.dictionary_path is None:
            return None
        path = self.dictionary_path
        if not path.exists():
            raise FileNotFoundError(f"Dictionary file not found for {self.name}: {path}")
        return path

    def resolve_synonyms(self) -> Optional[Path]:
        if self.synonyms_path is None:
            return None
        path = self.synonyms_path
        if not path.exists():
            raise FileNotFoundError(f"Synonym file not found for {self.name}: {path}")
        return path


@dataclass
class EntityRecord:
    entity_id: str
    primary_name: str
    ctype: str
    definition: Optional[str]
    aliases: List[str]
    source: str
    metadata: Dict[str, str] = field(default_factory=dict)

    def canonical_aliases(self) -> Iterable[str]:
        for alias in self.aliases:
            cleaned = alias.strip()
            if cleaned:
                yield cleaned


@dataclass
class DictionaryRecord:
    alias: str
    entity_id: str
    ctype: str
    source: str
    metadata: Dict[str, str] = field(default_factory=dict)


__all__ = ["NodeConfig", "EntityRecord", "DictionaryRecord"]
