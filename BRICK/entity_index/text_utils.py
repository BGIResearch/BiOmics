"""
功能：提供字符串检索相关的文本处理辅助函数（目前包含拼音转换），用于提升模糊匹配能力。
构建逻辑：尝试按需导入第三方库 `pypinyin`，若依赖缺失则自动降级为空实现，避免阻塞主流程。
数据流：上游为输入规范化或字符串检索模块；下游提供拼音 token 等结果，供 `string_client` 或其他检索组件使用。
调用链：混合检索时 `QueryNormalizer`、`StringRetriever` 等会调用这里的工具以生成拼音查询词。
"""

from __future__ import annotations

from typing import List

try:  # optional dependency
    from pypinyin import lazy_pinyin  # type: ignore
except ImportError:  # pragma: no cover - optional
    lazy_pinyin = None


def to_pinyin_tokens(text: str) -> List[str]:
    if not text:
        return []
    if lazy_pinyin is None:
        return []
    return lazy_pinyin(text)
