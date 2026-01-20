from token import OP
from pydantic import BaseModel, Field, model_validator
from langchain_core.messages import AnyMessage, BaseMessage
from langgraph.graph.message import add_messages
from datetime import datetime
from typing import Annotated, Literal, Union, List, Dict, Any, Optional, Tuple

from pynndescent.pynndescent_ import process_candidates

class QueryState(BaseModel):
    # 用户输入的原始问题
    source_question: Optional[str] = None
    # 翻译后的用户问题
    translated_question: Optional[str] = None
    # 选择脚本
    chosen_script: Optional[str|int] = None
    # 提取到的术语
    terminology: Optional[List[Tuple[str, str]]] = None
    # 任务3对应的目标类型
    target_type: Optional[str] = None
    # 混合检索结果
    searched_res:Optional[List[List]] = None
    # 经过rerank后大模型给出的实体参数
    entity_param: Optional[list] = None
    # 经过rerank后大模型给出的类型参数
    type_param: Optional[str] = None
    # 表格存储文件夹
    table_dir: Optional[str] = None
    # 关系查询结果表格
    relation_frame: Optional[Any] = None
    # 查询最终结果
    final_res:Optional[str] = None

    class Config:
        arbitrary_types_allowed = True  # 允许 DataFrame 等任意类型

