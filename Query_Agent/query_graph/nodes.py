import os
import sys
import subprocess
import time
from scanpy import read_10x_h5

# 添加项目根目录到 sys.path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from .create_agent import create_agent
from .state import QueryState
from langchain_core.messages import SystemMessage, HumanMessage
import json
import tempfile
from scripts import mrq, rq, tq
from datetime import datetime
from pathlib import Path
from langchain_core.output_parsers import JsonOutputParser, StrOutputParser
from .llm import basic_llm as model
from search_hybrid import search_entity




output_parser = JsonOutputParser()

def translator(state: QueryState) -> QueryState:
    print("translator started")
    agent = create_agent(
        "translator",
        "translator",
        [],
        "translator",
        state.model_dump()
    )
    
    result = agent.invoke({"messages": [HumanMessage(content="Please translate the query into English.")]})
    if "messages" in result:
        result = result["messages"][-1].content
        print("translator result:", result)
    
    try:
        result = json.loads(result)# type: ignore
        print("Translator json success.")
        print("translated question:",result["translated_question"])
        print("translator output:",result["output"])
    except json.JSONDecodeError:
        print("Translator json fail.")
        result = {}
    updated_state = state.model_copy(
        update={
            "translated_question": result.get("translated_question", []),
            "output": result.get("output", [])
        }
    )
    
    print("translator finished")
    return updated_state


def optant(state: QueryState) -> QueryState:
    print("optant started")
    agent = create_agent(
        "optant",
        "optant",
        [],
        "optant",
        state.model_dump()
    )
    result = agent.invoke({"messages": [HumanMessage(content="Please choose the correct script for the query.")]})
    if "messages" in result:
        result = result["messages"][-1].content
        print("optant result:", result)
    try:
        result = json.loads(result)# type: ignore
        print("optant json success.")
        print("chosen script:",result["chosen_script"])
        print("optant output:",result["output"])
    except json.JSONDecodeError:
        print("optant json fail.")
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    updated_state = state.model_copy(
        update={
            "chosen_script": result.get("chosen_script", []),
            "output": result.get("output", [])
        }
    )
    
    print("optant finished")
    return updated_state

def extractor(state: QueryState) -> QueryState:
    print("extractor started")
    agent = create_agent(
        "extractor",
        "extractor",
        [],
        "extractor",
        state.model_dump()
    )
    result = agent.invoke({"messages": [HumanMessage(content="Please extract the key information from the query.")]})
    if "messages" in result:
        result = result["messages"][-1].content
        print("extractor result:", result)
    
    try:
        result = json.loads(result)# type: ignore
        print("extractor json success.")
        # 将 JSON 数组转换为元组列表
        if "terminology" in result and isinstance(result["terminology"], list):
            result["terminology"] = [tuple(item) if isinstance(item, list) else item for item in result["terminology"]]
        print("terminology:",result["terminology"])
        print("target type:",result["target_type"])
        print("extractor output:",result["output"])
    except json.JSONDecodeError:
        print("extractor json fail.")
        result = {}

    updated_state = state.model_copy(
        update={
            "terminology": result.get("terminology", None),
            "target_type": result.get("target_type", None),
            "output": result.get("output", None)
        }
    )
    
    print("extractor finished")
    return updated_state

def hybrid_searcher(state: QueryState) -> QueryState:
    print("hybrid_searcher started")
    terms = state.terminology if state.terminology else []
    searched_res = search_entity(terms)
    out_content = "hybrid_searcher searched_res:" + str(searched_res)
    print(out_content)
    updated_state = state.model_copy(
        update={
            "searched_res": searched_res,
            "output": out_content
        }
    )
    
    print("hybrid_searcher finished")
    return updated_state

def reranker(state:QueryState) -> QueryState:
    print("reranker started")
    type_param = ""
    terminology = state.terminology 
    if terminology:
        for term in terminology:
            e, t = term
            type_param += t + "|"
    else:
        print("No terminology provided")
    # 去掉末尾的 "|"
    type_param = type_param.rstrip("|")
    print("type_param:", type_param)
    agent = create_agent(
        "reranker",
        "reranker",
        [],
        "reranker",
        state.model_dump()
    )
    result = agent.invoke({"messages": [HumanMessage(content="Please rerank the search results based on relevance to the query.")]})
    if "messages" in result:
        result = result["messages"][-1].content
        print("reranker result:", result)
    
    try:
        result = json.loads(result)# type: ignore
        print("reranker json success.")
        print("entity_param:",result["entity_param"])
        print("type_param:", type_param)
        print("reranker output:",result["output"])
    except json.JSONDecodeError:
        print("reranker json fail.")
        result = {}
    

    updated_state = state.model_copy(
        update={
            "entity_param": result.get("entity_param", []),
            "type_param": type_param,
            "output": result.get("output", [])
        }
    )
    
    print("reranker finished")
    return updated_state

def executor(state: QueryState) -> QueryState:
    print("executor started")
    chosen_script = state.chosen_script
    
    # 使用 table_dir 作为保存目录，默认回退到 output_files
    save_dir = state.table_dir if state.table_dir else "output_files"
    # 确保是绝对路径
    save_dir = os.path.abspath(save_dir)
    os.makedirs(save_dir, exist_ok=True)
    file_stamp = os.path.join(save_dir, str(int(time.time())))
    print(f"[DEBUG] CSV 将保存到: {save_dir}")
    
    ques = state.source_question if state.source_question else ""
    entity_param = state.entity_param if state.entity_param else []
    type_param = state.type_param if state.type_param else ""
    target_type = state.target_type if state.target_type else ""
    
    if chosen_script == 1 or chosen_script == "1":
        file_stamp += "mrq.csv"
        print("Executing most related query.")
        result, relations = mrq(question=ques, entity_name=entity_param, entity_type=type_param, target_type=target_type, file_path=file_stamp)

    elif chosen_script == 2 or chosen_script == "2":
        file_stamp += "rq.csv"
        print("Executing relation query.")
        result, relations = rq(question=ques, entity_name=entity_param, entity_type=type_param, file_path=file_stamp)

    elif chosen_script == 3 or chosen_script == "3":
        file_stamp += "tq.csv"
        print("Executing table query.")
        result, relations= tq(question=ques, entity_name=entity_param, query_type=type_param, file_path=file_stamp)

    else :
        result = "No validate script chosen"
        relations = None
    # 将 DataFrame 转为字典以便序列化
    relations_dict = relations.to_dict('records') if relations is not None and hasattr(relations, 'to_dict') else None
    updated_state = state.model_copy(
        update={
            "final_res": result,
            "relation_frame": relations_dict,
            "output": "agent finished the task,the response is that \n" + result
        }
    )
    print("final_res:", result)
    print("executor output: agent finished the task,the response is that \n" + result)
    print("executor finished")
    return updated_state

def general_responder(state: QueryState) -> QueryState:
    print("general_responder started")
    agent = create_agent(
        "general_responder",
        "general_responder",
        [],
        "general_responder",
        state.model_dump()
    )
    result = agent.invoke({"messages": [HumanMessage(content="Please generate a response to the query.")]})
    if "messages" in result:
        result = result["messages"][-1].content
        print("general_responder result:", result)

    try:
        result = json.loads(result)# type: ignore
        print("general_responder json success.")
        print("general_responder output:", result["output"])
    except json.JSONDecodeError:
        print("general_responder json fail.")
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    updated_state = state.model_copy(
        update={
            "final_res": "Agent can not finish the task.",
            "output": result.get("output")
        }
    )
    
    print("general_responder finished")
    return updated_state

