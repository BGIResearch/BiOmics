"""
功能：读取不同节点类型的词表文件（CSV/JSON/压缩格式），并标准化为 `EntityRecord`，供 catalog 与索引构建使用。
构建逻辑：结合节点配置解析字段映射，提取主键、标准名称、定义、别名及额外元信息，按迭代器逐条产出实体记录。
数据流：上游为 `catalog.build_catalog` 或其他需要实体基础信息的模块；下游则包括索引生成（`pipeline`、`index_builders`）以及检索模块。
调用链：典型路径是 `pipeline -> catalog -> loader`，在执行 `build_entity_records` 时读取配置并遍历词表。
"""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional

from .entity_types import DictionaryRecord, EntityRecord, NodeConfig


def _open_text(path: Path):
    """
    功能：根据文件后缀选择合适的文本读取方式，兼容普通文本与 gzip 压缩文件。
    构建逻辑：判断路径后缀是否为 `.gz`，若是则使用 `gzip.open` 以文本模式读取，否则直接调用 `Path.open`。
    出参入参：入参为文件路径 `Path`；出参为可作为上下文管理器的文件句柄。
    数据流：被 `iter_table_rows` 调用，作为读取实体表格数据的底层 I/O。
    模块内调用链：`iter_table_rows` -> `_open_text`。
    模块外调用链：`pipeline`/`catalog` 间接通过 `build_entity_records` 使用。
    """
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open(mode="r", encoding="utf-8")


def iter_table_rows(config: NodeConfig) -> Iterator[Dict[str, str]]:
    """
    功能：按行遍历节点配置对应的实体表格文件，输出每行的字典表示。
    构建逻辑：通过配置解析出表格路径，使用 `_open_text` 打开文件，并交给 `csv.DictReader` 逐行读取。
    出参入参：入参为节点配置 `NodeConfig`；出参为生成器，每次返回一行数据的字典。
    数据流：上游为 `build_entity_records` 传入的节点配置；下游为实体构建逻辑，将每行转换为 `EntityRecord`。
    模块内调用链：`build_entity_records` -> `iter_table_rows`。
    模块外调用链：`catalog.build_catalog` / `pipeline` 间接调用。
    """
    table_path = config.resolve_table()
    if table_path is None:
        raise ValueError(f"Node {config.name} does not have a table file configured")
    with _open_text(table_path) as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            yield row


def load_synonym_map(config: NodeConfig) -> Dict[str, List[str]]:
    """
    功能：加载与节点对应的别名映射表，将 JSON 文件中的 id→别名列表标准化为纯字符串列表。
    构建逻辑：读取 `NodeConfig` 中声明的 `synonyms_path`，解析 JSON，并过滤掉非字符串的条目。
    出参入参：入参为节点配置；出参为字典，key 为实体 ID，value 为别名列表。
    数据流：上游由 `build_entity_records` 在构建实体时调用；下游用于补充 `EntityRecord.aliases`。
    模块内调用链：`build_entity_records` -> `load_synonym_map`。
    模块外调用链：当 catalog 合并 overlay 别名时，也可能直接调用该函数。
    """
    path = config.resolve_synonyms()
    if path is None:
        return {}
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    normalized = {}
    for key, values in data.items():
        if isinstance(values, list):
            normalized[key] = [v for v in values if isinstance(v, str)]
    return normalized


def parse_aliases(raw: Optional[str]) -> List[str]:
    """
    功能：解析词表中管道符分隔的别名字符串，生成标准化列表。
    构建逻辑：处理空值或特殊占位符（如 "Undef"），否则按 `|` 拆分并过滤空项。
    出参入参：入参为原始别名字串；出参为别名单元列表。
    数据流：在 `build_entity_records` 中用于解析 CSV 的 synonym 列；下游继续合并外部别名并去重。
    模块内调用链：`build_entity_records` -> `parse_aliases`。
    模块外调用链：无。
    """
    if not raw or raw.strip() in {"Undef", "undef"}:
        return []
    return [alias for alias in raw.split("|") if alias]


def build_entity_records(config: NodeConfig) -> Iterable[EntityRecord]:
    """
    功能：根据节点配置加载词表数据与别名映射，构建标准化的 `EntityRecord` 序列。
    构建逻辑：
    1. 先加载 JSON 别名映射与描述字段；
    2. 遍历 CSV 表格行，提取实体 ID、名称、类型、定义、别名；
    3. 合并外部别名并去重；
    4. 组装 metadata（包含额外字段与配置补充信息）；
    5. 产出 `EntityRecord`。
    出参入参：入参为 `NodeConfig`；出参为可迭代的 `EntityRecord` 序列。
    数据流：上游由 `catalog.build_catalog` 调用，最终供索引构建与检索模块消费。
    模块内调用链：内部依次调用 `load_synonym_map`、`iter_table_rows`、`parse_aliases`。
    模块外调用链：`catalog`、`pipeline`、`search` 等使用此函数获取实体数据。
    """
    synonyms_map = load_synonym_map(config)
    description_fields = config.description_fields or []

    for row in iter_table_rows(config):
        entity_id = row.get(config.columns["id"], "").strip()
        name = row.get(config.columns["name"], "").strip()
        ctype = row.get(config.columns.get("type", "type"), config.name).strip() or config.name
        definition = row.get(config.columns.get("definition", "def"))
        definition = definition.strip() if isinstance(definition, str) else None
        if definition in {"", "Undef"}:
            definition = None
        raw_aliases = row.get(config.columns.get("synonyms", "synonym"))
        aliases = parse_aliases(raw_aliases)

        aliases.extend(synonyms_map.get(entity_id, []))

        deduped_aliases: List[str] = []
        seen = set()
        for alias in aliases:
            cleaned = alias.strip()
            if not cleaned:
                continue
            if cleaned.lower() == name.lower() or cleaned == entity_id:
                continue
            lowered = cleaned.lower()
            if lowered in seen:
                continue
            seen.add(lowered)
            deduped_aliases.append(cleaned)

        metadata: Dict[str, str] = {}
        metadata.update(config.extra_metadata)
        for field_name in description_fields:
            value = row.get(field_name)
            if value:
                metadata[field_name] = value

        yield EntityRecord(
            entity_id=entity_id,
            primary_name=name,
            ctype=ctype,
            definition=definition,
            aliases=deduped_aliases,
            source=config.name,
            metadata=metadata,
        )


def load_dictionary_records(config: NodeConfig) -> List[DictionaryRecord]:
    """
    功能：读取聚合字典（alias -> entity_id）映射，生成标准化的 `DictionaryRecord` 列表。
    构建逻辑：
    1. 解析配置中的 `dictionary_path`；
    2. 遍历 JSON 键值对，过滤空 alias 或 target；
    3. 根据目标 ID 推断类型（默认使用冒号前缀），并记录来源信息。
    出参入参：入参为 `NodeConfig`；出参为聚合字典记录列表。
    数据流：聚合索引构建或检索模块可调用本函数，将别名映射写入索引或本地缓存。
    模块内调用链：无。
    模块外调用链：`pipeline` 在构建 all_name2id 索引时使用。
    """
    path = config.resolve_dictionary()
    if path is None:
        raise ValueError(f"Node {config.name} does not declare a dictionary_path")

    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    records: List[DictionaryRecord] = []
    for alias, target in data.items():
        if not isinstance(alias, str) or not alias.strip():
            continue
        if not isinstance(target, str) or not target.strip():
            continue
        alias_clean = alias.strip()
        target_clean = target.strip()
        if ":" in target_clean:
            ctype = target_clean.split(":", 1)[0]
        else:
            ctype = config.extra_metadata.get("target_type", config.name)
        records.append(
            DictionaryRecord(
                alias=alias_clean,
                entity_id=target_clean,
                ctype=ctype,
                source=config.name,
                metadata=config.extra_metadata.copy(),
            )
        )
    return records
