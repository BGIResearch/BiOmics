"""
功能：向调用方统一暴露实体索引构建相关的核心对象，避免直接依赖旧有的 RAG/LLM 代码。
构建逻辑：通过在包级导入节点配置、实体类型定义、嵌入客户端与索引管线等组件，实现开箱即用的入口。
数据流：上游为需要加载实体 catalog 或执行索引写入/检索的业务模块；下游分别由 `config`、`types`、`embedding`、`pipeline` 等子模块提供数据读取、嵌入生成与索引组装能力。
调用链：在业务侧执行 `import BRICK.entity_index` 时，先进入本模块，再按需引用各子模块以完成实体索引的构建或读写。
"""

from .config import NODE_CONFIGS, NodeConfig  # noqa: F401
from .entity_types import EntityRecord  # noqa: F401
from .embedding import (  # noqa: F401
    EmbeddingClient,
    EmbeddingConfig,
    build_batch_embedder,
)
from .pipeline import build_index_artifacts  # noqa: F401
