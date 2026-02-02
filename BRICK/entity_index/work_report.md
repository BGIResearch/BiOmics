# 索引搭建工作报告

## 全局
- 依据 `target_plan.md` 更新节点检索策略记录（`target_plan.md:30-41`），确保实际实现与规划一致。

## 新增模块 `entity_index`
- `entity_index/__init__.py`
  - 目的：提供包级导出，外部可直接引用配置与实体类型。
- `entity_index/types.py`
  - 目的：定义 `NodeConfig`、`EntityRecord` 数据结构，描述每个节点的索引配置与实体基础字段。
- `entity_index/config.py`
  - 目的：集中列出全部节点的配置（表路径、语义开关、overlay-only 设定等），作为索引构建的单一事实来源。
- `entity_index/loader.py`
  - 目的：读取 CSV/CSV.GZ/JSON 词表并清洗别名，构建 `EntityRecord` 列表，自动剔除空值与重复。
- `entity_index/catalog.py`
  - 目的：合并 overlay 别名（如 Stage）并按节点输出实体清单，供后续索引写入使用。
- `entity_index/index_builders.py`
  - 目的：提供字符串/向量索引的 ES mappings 与文档序列化逻辑（edge-ngram、search-as-you-type、拼音可选、dense_vector + norm 等）。
- `entity_index/text_utils.py`
  - 目的：封装拼音转换为可选依赖（缺省时自动返回空列表）。
- `entity_index/pipeline.py`
  - 目的：将上述组件整合成 `IndexArtifacts`，统一输出字符串/向量索引所需的 settings 与文档集合，并在启用语义检索时要求外部提供 embedding 函数。

## 验证
- `python -m compileall entity_index`
  - 目的：快速校验新模块语法正确，确保引入仓库后可直接使用。

> 所有改动均发生在新建的 `entity_index` 目录内，未触碰既有 RAG/LLM 问答与知识图谱模块，符合 `Inhibition.md` 约束。

## 2025-10-27 更新
- `entity_index/entity_types.py`: 为 `NodeConfig` 扩展 `dictionary_path`/`is_dictionary`，新增 `DictionaryRecord` 以支撑聚合字典。
- `entity_index/config.py`: 新增 `ALL_NAME2ID` 节点配置，并将 Stage 切换为常规 CSV 节点，便于独立写入字符串索引。
- `entity_index/loader.py`: 增补 `load_dictionary_records`，从 `all_name2id.json` 构造聚合记录。
- `entity_index/index_builders.py`: 新增 `dictionary_index_settings` 与 `dictionary_document`，输出 alias→ID 索引文档。
- `entity_index/catalog.py` / `entity_index/pipeline.py`: 兼容 `is_dictionary` 节点，允许聚合字典与常规实体并行生成索引工件。
- `entity_index/search/settings.py`: 暴露 `alias_index_name`，支持通过 `HYBRID_ALIAS_INDEX`/`ES_INDEX_PREFIX` 自动定位聚合字典索引。
- `SearchClient.py`: 新增 `search_alias_dictionary` 方法，对外提供别名归并查询；初始化时默认指向远端 `brick_index_all_name2id_string`。
- `entity_index/search.executed.ipynb`: 添加 “聚合字典检索（all_name2id）” 示例单元，演示 `EGFR` 别名在聚合索引中的归并效果。
- `.env`: 将 `ES_INDEX_PREFIX` 固定为 `brick_index` 并补充 `HYBRID_ALIAS_INDEX`，避免再次触发通配符写入错误。
