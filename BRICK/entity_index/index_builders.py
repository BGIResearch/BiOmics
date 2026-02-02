"""
功能：提供实体索引的字符串/向量文档构造函数以及对应的 Elasticsearch settings/ mappings，使得 pipeline 能快速拼装索引工件。
构建逻辑：定义字符串文档的字段结构、向量文档的 embedding 存储格式，并封装生成 settings 的辅助函数。
数据流：上游为 `pipeline.build_index_artifacts`、`catalog` 等模块提供的实体记录；下游则是写入器或检索模块使用这些文档来创建索引。
调用链：索引构建时 pipeline 调用本模块的 `string_document`、`vector_document`、`string_index_settings`、`vector_index_settings` 等函数。
"""

from __future__ import annotations

import math
from typing import Any, Dict, Iterable, Iterator, List, Optional

from .text_utils import to_pinyin_tokens
from .entity_types import DictionaryRecord, EntityRecord


EDGE_NGRAM_MIN = 1
EDGE_NGRAM_MAX = 20


def string_index_settings() -> Dict[str, Any]:
    """
    功能：生成字符串索引所需的 Elasticsearch settings 与 mappings，覆盖精确、模糊、前缀检索及拼音字段。
    构建逻辑：设置 edge_ngram 分析器、search_as_you_type 字段，并定义实体文档的属性类型。
    出参入参：无入参；返回包含 settings/mappings 的 dict。
    数据流：被 `pipeline.build_index_artifacts` 调用，用于创建字符串索引结构；写入器将据此调用 ES 创建索引。
    模块内调用链：无。
    模块外调用链：`pipeline` -> `string_index_settings` -> `es_writer`。
    """
    return {
        "settings": {
            "analysis": {
                "filter": {
                    "edge_ngram_filter": {
                        "type": "edge_ngram",
                        "min_gram": EDGE_NGRAM_MIN,
                        "max_gram": EDGE_NGRAM_MAX,
                    }
                },
                "analyzer": {
                    "edge_ngram_analyzer": {
                        "type": "custom",
                        "tokenizer": "standard",
                        "filter": ["lowercase", "asciifolding", "edge_ngram_filter"],
                    }
                },
            }
        },
        "mappings": {
            "dynamic": "false",
            "properties": {
                "entity_id": {"type": "keyword"},
                "primary_name": {"type": "keyword"},
                "ctype": {"type": "keyword"},
                "source": {"type": "keyword"},
                "definition": {"type": "text"},
                "aliases": {"type": "keyword"},
                "search_terms": {
                    "type": "text",
                    "analyzer": "standard",
                    "fields": {
                        "ngram": {"type": "text", "analyzer": "edge_ngram_analyzer"},
                        "prefix": {"type": "search_as_you_type"},
                    },
                },
                "pinyin_terms": {"type": "keyword"},
                "metadata": {"type": "object", "enabled": False},
            },
        },
    }


def dictionary_index_settings() -> Dict[str, Any]:
    """
    功能：为聚合字典构建字符串索引 settings/mappings，仅包含别名检索所需字段。
    构建逻辑：复用字符串索引的 edge_ngram/prefix 设置，字段精简为 alias/target 信息。
    出参入参：无入参；返回 settings/mappings dict。
    数据流：聚合索引构建前由 `pipeline` 调用。
    模块内调用链：无。
    模块外调用链：`pipeline` -> `dictionary_index_settings` -> `es_writer`。
    """
    base = string_index_settings()
    base["mappings"]["properties"] = {
        "alias": {"type": "keyword"},
        "entity_id": {"type": "keyword"},
        "ctype": {"type": "keyword"},
        "source": {"type": "keyword"},
        "search_terms": {
            "type": "text",
            "analyzer": "standard",
            "fields": {
                "ngram": {"type": "text", "analyzer": "edge_ngram_analyzer"},
                "prefix": {"type": "search_as_you_type"},
            },
        },
        "pinyin_terms": {"type": "keyword"},
        "metadata": {"type": "object", "enabled": False},
    }
    return base


def string_document(record: EntityRecord) -> Dict[str, Any]:
    """
    功能：将 `EntityRecord` 转换为字符串索引的文档，包含主名称、别名、拼音 token 等字段。
    构建逻辑：组装搜索词列表，调用 `to_pinyin_tokens` 生成拼音并去重，写入文档。
    出参入参：入参为实体记录；出参为可写入 ES 的字典。
    数据流：`pipeline` 在生成字符串索引文档列表时调用，结果用于 bulk 写入或其他缓存。
    模块内调用链：依赖 `to_pinyin_tokens`。
    模块外调用链：`pipeline` -> `string_document` -> `es_writer`。
    """
    search_terms = [record.primary_name] + list(record.aliases)
    pinyin_terms: List[str] = []
    for term in search_terms:
        pinyin_terms.extend(to_pinyin_tokens(term))
    # deduplicate while preserving order
    seen = set()
    dedup_pinyin: List[str] = []
    for token in pinyin_terms:
        lowered = token.lower()
        if lowered in seen:
            continue
        seen.add(lowered)
        dedup_pinyin.append(token)

    return {
        "entity_id": record.entity_id,
        "primary_name": record.primary_name,
        "ctype": record.ctype,
        "source": record.source,
        "definition": record.definition,
        "aliases": record.aliases,
        "search_terms": search_terms,
        "pinyin_terms": dedup_pinyin,
        "metadata": record.metadata,
    }


def dictionary_document(record: DictionaryRecord) -> Dict[str, Any]:
    """
    功能：将 `DictionaryRecord` 转换为聚合字典索引文档，保留别名与目标实体信息。
    构建逻辑：以别名作为 search_terms，生成拼音 token，合并元数据字段。
    出参入参：入参为字典记录；出参为写入 ES 的字典。
    数据流：聚合字典索引构建阶段由 `pipeline` 调用。
    模块内调用链：依赖 `to_pinyin_tokens`。
    模块外调用链：`pipeline` -> `dictionary_document` -> `es_writer`。
    """
    search_terms = [record.alias]
    pinyin_terms: List[str] = []
    for term in search_terms:
        pinyin_terms.extend(to_pinyin_tokens(term))
    seen = set()
    dedup_pinyin: List[str] = []
    for token in pinyin_terms:
        lowered = token.lower()
        if lowered in seen:
            continue
        seen.add(lowered)
        dedup_pinyin.append(token)

    return {
        "alias": record.alias,
        "entity_id": record.entity_id,
        "ctype": record.ctype,
        "source": record.source,
        "search_terms": search_terms,
        "pinyin_terms": dedup_pinyin,
        "metadata": record.metadata,
    }


def vector_index_settings(dimension: int) -> Dict[str, Any]:
    """
    功能：生成向量索引的 Elasticsearch mappings，定义 dense_vector 字段及实体元数据字段。
    构建逻辑：指定向量维度、相似度、附加元数据字段类型。
    出参入参：入参为 embedding 维度；出参为包含 mappings 的 dict。
    数据流：`pipeline` 调用后传递给写入器创建向量索引。
    模块内调用链：无。
    模块外调用链：`pipeline` -> `vector_index_settings` -> `es_writer`。
    """
    return {
        "mappings": {
            "dynamic": "false",
            "properties": {
                "entity_id": {"type": "keyword"},
                "primary_name": {"type": "keyword"},
                "ctype": {"type": "keyword"},
                "source": {"type": "keyword"},
                "definition": {"type": "text"},
                "aliases": {"type": "keyword"},
                "metadata": {"type": "object", "enabled": False},
                "embedding": {
                    "type": "dense_vector",
                    "dims": dimension,
                    "index": True,
                    "similarity": "cosine",
                },
                "embedding_norm": {"type": "float"},
            },
        },
    }


def compute_embedding_norm(vector: Iterable[float]) -> float:
    """
    功能：计算 embedding 的 L2 范数，供余弦相似度或向量归一化使用。
    构建逻辑：遍历向量元素求平方和并开方。
    出参入参：入参为向量迭代器；出参为 float 型范数。
    数据流：`vector_document` 在缺省范数时调用；结果随文档写入索引。
    模块内调用链：`vector_document` -> `compute_embedding_norm`。
    模块外调用链：无。
    """
    total = 0.0
    for value in vector:
        total += value * value
    return math.sqrt(total)


def vector_document(
    record: EntityRecord,
    embedding: List[float],
    embedding_norm: Optional[float] = None,
) -> Dict[str, Any]:
    """
    功能：构造单条向量索引文档，将实体元数据与 embedding 组合写入。
    构建逻辑：若未传入范数则调用 `compute_embedding_norm`，随后输出包含 embedding/metadata 的字典。
    出参入参：入参为实体记录、向量及可选范数；出参为向量文档。
    数据流：`pipeline` 为每个实体调用，结果供写入器 bulk 导入。
    模块内调用链：调用 `compute_embedding_norm`（当范数缺失时）。
    模块外调用链：`pipeline` -> `vector_document` -> `es_writer`。
    """
    if embedding_norm is None:
        embedding_norm = compute_embedding_norm(embedding)
    return {
        "entity_id": record.entity_id,
        "primary_name": record.primary_name,
        "ctype": record.ctype,
        "source": record.source,
        "definition": record.definition,
        "aliases": record.aliases,
        "metadata": record.metadata,
        "embedding": embedding,
        "embedding_norm": embedding_norm,
    }


def iter_bulk_actions(
    index_name: str,
    documents: Iterable[Dict[str, Any]],
) -> Iterator[Dict[str, Any]]:
    """
    功能：生成 Elasticsearch bulk API 需要的动作/文档序列。
    构建逻辑：为每个文档前置一条 `{"index": {"_index": index_name}}` 动作，再输出文档本身。
    出参入参：入参为索引名和文档迭代器；出参为迭代器（动作与文档交替）。
    数据流：写入器调用该函数准备 bulk 请求，最终由 ES 客户端发送。
    模块内调用链：无。
    模块外调用链：`es_writer`。
    """
    for doc in documents:
        yield {"index": {"_index": index_name}}
        yield doc


__all__ = [
    "compute_embedding_norm",
    "dictionary_document",
    "dictionary_index_settings",
    "iter_bulk_actions",
    "string_document",
    "string_index_settings",
    "vector_document",
    "vector_index_settings",
]
