import os
import subprocess

from scanpy import read_10x_h5
from .Debugger_simplified import debug
from tools.mock_adata import generate_enhanced_simulated_data_tool
from tools.summarize_data import summarize_data_tool
from tools.valid_h5ad import validate_h5ad_structure_tool
from tools.write_file import write_file_tool
from tools.read_table import read_table_tool
from Query_Agent import run_query
from .create_agent import create_agent
from .state import BrickState
from langchain_core.messages import SystemMessage, HumanMessage
import json
import tempfile
from datetime import datetime
from pathlib import Path
from utils.sandbox_client import DockerSandbox
from utils.exec_res_tool import (
    format_stdout, format_stderr, format_result,
    format_images, format_error, has_error
)
from utils.docker_path_transform import transform_path
from utils.save_dir_name import get_save_dir
from utils.gene_id_converter import convert_var_names_to_symbol
from tools.extract_code import extract_python_blocks_tool
from tools.read_notebook import read_notebook_tool
from tools.read_notebook import extract_notebook_text
from tools.brick_rag_searcher import perform_rag_search_tool
from langchain_core.output_parsers import JsonOutputParser, StrOutputParser
from .output import safe_stream_model
from .llm import basic_llm as model
from .RAG import RAG, RAG_func
from .notebook_finder import find_notebook
output_parser = JsonOutputParser()

def env_checker(state: BrickState) -> BrickState:

    data_path = state.data_path
    if state.update_data_info:
        state.messages.append(HumanMessage(content=state.update_data_info))

    agent = create_agent(
        "env_checker",
        "env_checker",
        [
            validate_h5ad_structure_tool,
            summarize_data_tool,
            read_table_tool,
        ],
        "env_checker",
        state.model_dump()
    )
    invoke_message = state.messages[-2:]
    print("message list:",state.messages)
    result = agent.invoke({"messages":invoke_message})
    #result = safe_stream_agent_json(agenmessages[-1]t, {"messages":state.messages[-1]})
    

    if "messages" in result:
        result = result["messages"][-1].content
        print("env_checker result:", result)
        result = json.loads(result)
    
    print("env_checker result:", result)
    if isinstance(result, dict):
        print("[thought]", result["thought"])
        print("[output]", result["output"])

        
            
        new_message = SystemMessage(content=f"thought: {result['thought']}, output: {result['output']}, status: {result['status']}")
        
        # 先根据原始文件获取 data_info
        if data_path:
            if data_path.endswith('.h5ad'):
                data_info = summarize_data_tool.invoke({"data_path": data_path})["data_info"]
            elif data_path.endswith('.csv') or data_path.endswith('.tsv'):
                data_info = read_table_tool.invoke({"file_path": data_path})["description"]
            else:
                data_info = f"Unsupported file type: {data_path}"
        else:
            data_info = "No target data found."
        
        
        # 转换 docker 路径
        docker_data_path = transform_path(data_path) if data_path else None
        
        print("di",data_info)
        updated_state = state.model_copy(
            update={
                "messages": state.messages + [new_message],
                "thought": result["thought"],
                "output": result["output"],
                "next": result["next"],
                "status": result["status"],
                "data_path": data_path,
                "docker_data_path": docker_data_path,
                "data_info": data_info, 
                "make_env_check" : True,
                "agent":"env_checker"
            },
        )

    else:
        error_content = f"Agent returned an unexpected type: {type(result)}. Full result: {result}"
        new_message = SystemMessage(content=error_content)
        updated_state = state.model_copy(
            update={
                "messages": state.messages + [new_message],
                "output": error_content,
                "status": "ERROR",
                "next": "env_checker" # 或者设置为一个错误处理节点
            }
        )

    return updated_state

def supervisor(state: BrickState) -> BrickState:

    box = DockerSandbox(server_url="127.0.0.1:8888", workspace="/workspace")
    
    if not box.start():
        print("沙箱启动失败")
        updated_state = state.model_copy(
        update={
            "thought": "Failed to start the Docker sandbox",
            "output": "Please check whether the Docker image is running and then ask again",
            "agent": "supervisor",
            "next": "general_responder"
        })
        return updated_state
        exit(1)
    
    print(f"\n沙箱 ID: {box.get_id()}")
    sandbox_id = box.get_id()
    agent = create_agent(
        "supervisor",
        "supervisor",
        [],
        "supervisor",
        state.model_dump()
    )
    
    print("supervisor invoke message")
    result = agent.invoke({"messages": state.messages[-1]})
    print("supervisor result")
    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        # Fallback: keep raw content so UI can display something
        result = {"output": result, "status": "NOT_FINISHED"}
    
    if not isinstance(result, dict):
        result = {}
    if state.data_path:
        docker_data_path = transform_path(state.data_path)
    else:
        docker_data_path = None
    
    sd = get_save_dir(state.question)
    docker_save_dir = transform_path(sd)
    print("[thought]", result.get("thought", []))
    print("[output]", result.get("output", []))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', [])}, output: {result.get('output', [])}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", []),
            "next": result.get("next", []),
            "status": "NOT_FINISHED",
            "docker_data_path": docker_data_path,
            "sandbox_id": sandbox_id,
            "output": result.get("output", []),
            "agent": "supervisor",
            "language": result.get("language", "English"),
            "save_dir": sd,
            "docker_save_dir": docker_save_dir, 
        }
    )

    return updated_state

def translator(state: BrickState) -> BrickState:

    
    agent = create_agent(
        "translator", 
        "translator", 
        [], 
        "translator",
        state.model_dump()
    )
    invoke_message = state.messages[-1]
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("translator result:", result)
    print("translator result type:", type(result))
    
    if "messages" in result:

        result = result["messages"][-1].content
    result = json.loads(result)
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", []))
    print("[output]", result.get("output", []))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', [])}, question: {result.get('output', [])}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", []),
            "translated_question": result.get("output", []),
            "question": result.get("original_question", ""),
            "next": result.get("next", []),
            "status": "NOT_FINISHED",
            "agent": "translator"
        }
    )
    

    return updated_state

def data_analyzer(state: BrickState) -> BrickState:
    print("data_analyzer received state:")

    if state.update_data_repo:
        state.messages.append(HumanMessage(content=state.update_data_repo))

    agent = create_agent("data_analyzer", "data_analyzer", [], "data_analyzer", state.model_dump())

    print("data_analyzer invoke message:")
    result = agent.invoke({"messages":[state.messages[-1]]})
    print("data_analyzer result:")

    if "messages" in result:
        result = result["messages"][-1].content
        print("json loads:",result)
    result = json.loads(result)

    if isinstance(result, dict):
        print("[thought]", result["thought"])
        print("[output]", result["output"])

        new_message = SystemMessage(content=f"thought: {result['thought']}, data_repo: {result['output']}")
        if result["status"] == "Revise":
            updated_state = state.model_copy(update={
                "messages": state.messages + [new_message],
                "thought": result["thought"],
                "data_repo": result["output"],
                "output": result["output"],
                "status": result["status"],
                "next": "data_analyzer",
                "agent": "data_analyzer"
            })
        else:
            updated_state = state.model_copy(update={
                "messages": state.messages + [new_message],
                "thought": result["thought"],
                "data_repo": result["output"],
                "output": result["output"],
                "make_data_repo": True, 
                "final_result": result["output"],
                "status": "VALIDATED",
                "next": "analyze_planner",
                "agent": "data_analyzer"
            })
    else:
        error_content = f"Agent returned an unexpected type: {type(result)}. Full result: {result}"
        new_message = SystemMessage(content=error_content)
        updated_state = state.model_copy(
            update={
                "messages": state.messages + [new_message],
                "output": error_content,
                "status": "ERROR",
                "next": "data_analyzer" # 或者设置为一个错误处理节点
            }
        )


    #updated_state.messages.append(new_message)
    return updated_state

def analyze_planner(state: BrickState) -> BrickState:

    agent = create_agent(
        "analyze_planner", 
        "analyze_planner", 
        [], 
        "analyze_planner",
        state.model_dump()
    )

    result = agent.invoke({"messages": [HumanMessage(content="请为用户问题{user_question}生成分析计划")]})
    print("analyze_planner result:", result)
    notebook_content = read_notebook_tool.invoke({
    "agentname": "analyze_planner",
    "notebook_path": state.reference_notebook_path
    }).cells_text
    
    if "messages" in result:
        result = result["messages"][-1].content

    print("analyze_planner result2:", result)
    result = json.loads(result)

    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    
    if result.get("status") == "ASK_USER":
        new_message = SystemMessage(content=f"thought: {result.get('thought', '')}")
        updated_state = state.model_copy(
            update={
                "messages": state.messages + [new_message],
                "thought": result.get("thought", ""),
                "output": result.get("output", ""),
                "a_plan": result.get("output", ""),
                "status": "ASK_USER",
                "next": "analyze_planner",
                "agent": "analyze_planner"
            }
        )
    else:
        new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, a_plan: {result.get('output', '')}")
        updated_state = state.model_copy(
            update={
                "messages": state.messages + [new_message],
                "thought": result.get("thought", ""),
                "output": result.get("output", ""),
                "a_plan": result.get("output", ""),
                "status": "VALIDATED",
                "next": "planner",
                "make_analysis": True,
                "agent": "analyze_planner"
            }
        )

    return updated_state

def planner(state: BrickState) -> BrickState:

    
    agent = create_agent(
        "planner", 
        "planner", 
        [perform_rag_search_tool], 
        "planner",
        state.model_dump()
    )
    print("planner invoke message:")
    result = agent.invoke({"messages": state.messages})
    
    print("planner result:")
    print("planner result type:")
    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", []))
    step_num = len(result.get("output"))
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, current_plan: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "current_plan": result.get("output", ""),
            "next": result.get("next", ""),
            "status": "NOT_FINISHED",
            "make_plan": True,
            "agent": "planner",
            "step_num": step_num
        }
    )
    

    return updated_state

def code_debugger(state: BrickState) -> BrickState:

    agent = create_agent(
        "code_debugger", 
        "code_debugger", 
        [], 
        "code_debugger",
        state.model_dump()
    )
    result = agent.invoke({"messages": [HumanMessage(content="请按照模板要求，根据错误信息修改代码")]})

    print("code_debugger result:", result)
    print("code_debugger result type:", type(result))
    if "messages" in result:
        result = result["messages"][-1].content
    print("Content to parse:", repr(result))
    print("Content length:", len(result))
    
    # 尝试提取JSON部分

    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        # Fallback: keep raw content so UI can display something
        result = {"output": result, "status": "NOT_FINISHED"}
    
    if not isinstance(result, dict):
        result = {}
    print("code_debugger result2:", result)
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    new_full_code = state.full_code
    new_full_code[-1] = result.get("output", "")
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, full_code: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "full_code": new_full_code,
            "code":result.get("output", ""),
            "error_code": result.get("output", ""),
            "next": "code_runner",
            "status": "NOT_FINISHED",
            "agent": "code_debugger"
        }
    )
    

    return updated_state


def plan_reviewer(state: BrickState) -> BrickState:

    
    if state.update_instruction:
        state.messages.append(HumanMessage(content=state.update_instruction[-1]))
    
    agent = create_agent(
        "plan_reviewer", 
        "plan_reviewer", 
        [], 
        "plan_reviewer",
        state.model_dump()
    )
    print("plan_reviewer invoke message:")
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("plan_reviewer result:")
    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, re_plan: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "re_plan": result.get("output", ""),
            "next": result.get("next", ""),
            "status": result.get("status", ""),
            "agent": "plan_reviewer"
        }
    )
    

    return updated_state

def responder(state: BrickState) -> BrickState:

    agent = create_agent(
        "responder", 
        "responder", 
        [], 
        "responder",
        state.model_dump()
    )
    print("responder invoke message:")
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("responder result:")

    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[answer]", result.get("output", ""))
    
    new_message = SystemMessage(content=f"final_answer: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "final_answer": result.get("output", ""),
            "status": "FINISHED",
            "agent": "responder"
        }
    )
    return updated_state
def plan_executor(state: BrickState) -> BrickState:
    plan = state.current_plan
    step_num = state.step_num
    current_step = state.current_step
    is_end = current_step >= step_num
    step_content = None
    rag_res = None
    op = ""
    print("step_num:", step_num, "current_step:", current_step)
    if not is_end:
        step_content = plan[current_step]
        print("== starting RAG ==")
        rag_res = RAG_func(query=json.dumps(step_content),vectorstore=state.default_vectorstore)
        print("== rag_res ==", rag_res)
        current_step += 1
        next_agent = "coder"
        op = "Executing step " + str(current_step) + " of " + str(step_num) + "..."
    if is_end:
        op = "All steps executed."
        next_agent = "responder"

    updated_state = state.model_copy(
        update={
            "thought": "pe thinking",
            "output": op,
            "next": next_agent,
            "step_content": step_content,
            "find_function": rag_res,
            "status": "NOT_FINISHED",
            "agent": "plan_executor",
            "current_step": current_step
        }
    )
    return updated_state
def plan_executor2(state: BrickState) -> BrickState:

    
    agent = create_agent(
        "plan_executor", 
        "plan_executor", 
        [], 
        "plan_executor",
        state.model_dump()
    )
    print("plan_executor invoke message:")
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("plan_executor result:")
    plan = state.current_plan
    print("tipe",type(plan))


    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    print("plan_executor result2:", result)
    
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    print("== starting RAG ==")
    if result.get("step_content", ""):
        rag_res = RAG_func(query=json.dumps(result.get("step_content", "")),vectorstore=state.default_vectorstore)
    else:
        rag_res = ""
    print("== rag_res ==", rag_res)
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, output: {result.get('output', '')}")
    updated_checked_step = state.checked_step
    if len(updated_checked_step) == 0:
        updated_checked_step.append(1)
    else:
        updated_checked_step.append(updated_checked_step[-1] + 1)

    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "checked_step": updated_checked_step,
            "next": result.get("next", ""),
            "step_content": result.get("step_content", ""),
            "find_function": rag_res,
            "status": "NOT_FINISHED",
            "agent": "plan_executor"
        }
    )
    

    return updated_state

def searcher(state: BrickState) -> BrickState:

    
    agent = create_agent(
        "searcher", 
        "searcher", 
        [], 
        "searcher",
        state.model_dump()
    )
    print("searcher invoke message:")
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("searcher result:", result)
    print("searcher result type:", type(result))
    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, output: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "next": result.get("next", ""),
            "status": result.get("status", ""),
            "agent": "searcher"
        }
    )
    

    return updated_state

def bkg_searcher(state: BrickState) -> BrickState:

    
    agent = create_agent(
        "bkg_searcher", 
        "bkg_searcher", 
        [], 
        "bkg_searcher",
        state.model_dump()
    )
    print("bkg_searcher invoke message:", state.messages[-1])
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("bkg_searcher result:", result)
    print("bkg_searcher result type:", type(result))
    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, output: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "next": result.get("next", ""),
            "status": result.get("status", ""),
            "agent": "bkg_searcher"
        }
    )

    return updated_state

def coder(state: BrickState) -> BrickState:

    
    agent = create_agent(
        "coder", 
        "coder", 
        [], 
        "coder",
        state.model_dump()
    )
    invoke_message = state.messages
    print("coder invoke message:")
    result = agent.invoke({"messages": invoke_message})
    
    print("coder result:")

    
    if "messages" in result:

        result = result["messages"][-1].content
    # Robust JSON parsing with fallback to raw text
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {"thought": "", "output": str(result), "status": "CODE_CONFIRMATION", "next": "code_runner"}
    if isinstance(result, dict):
        print("[thought]", result.get("thought", ""))
        print("[output]", result.get("output", ""))
        new_full_code = state.full_code
        new_full_code.append(result["output"])
        new_message = SystemMessage(role="assistant", content=f"thought: {result['thought']}, code: {result['output']}")
        updated_state = state.model_copy(update={
            "messages": state.messages + [new_message],
            "thought": result["thought"],
            "output": result.get("output", ""),
            "full_code": new_full_code,
            "code": result["output"],
            "error_code": result["output"],
            "next": "code_runner",
            "status": "NOT_FINISHED",
            "agent": "coder"
        })
    else:
        error_content = f"Agent returned an unexpected type: {type(result)}. Full result: {result}"
        new_message = SystemMessage(role="assistant", content=error_content)
        updated_state = state.model_copy(
            update={
                "messages": state.messages + [new_message],
                "output": error_content,
                "status": "ERROR",
                "next": "coder" # 或者设置为一个错误处理节点
            }
        )

    return updated_state
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, output: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "next": result.get("next", ""),
            "status": result.get("status", "")
        }
    )

    return updated_state

def code_controller(state: BrickState) -> BrickState:

    
    agent = create_agent(
        "code_controller", 
        "code_controller", 
        [], 
        "code_controller",
        state.model_dump()
    )
    print("code_controller invoke message:", state.messages[-1])
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("code_controller result:", result)
    print("code_controller result type:", type(result))
    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, output: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "next": result.get("next", ""),
            "status": result.get("status", ""),
            "agent": "code_controller"
        }
    )
    

    return updated_state

def code_evaluator_old(state: BrickState) -> BrickState:

    agent = create_agent(
        "code_evaluator", 
        "code_evaluator", 
        [], 
        "code_evaluator",
        state.model_dump()
    )
    print("code_evaluator invoke message:", state.messages[-1])
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("code_evaluator result:", result)
    print("code_evaluator result type:", type(result))
    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, output: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "next": result.get("next", ""),
            "status": result.get("status", "")
        }
    )
    

    return updated_state

def code_evaluator(state: BrickState):
    if state.status == "CODE_CONFIRMATION":

        #logger.info(f"[code_evaluator] update_code: {state.update_code}")
        with open('code_evaluator.md', 'r', encoding='utf-8') as f:
            template = f.read().format(**state.model_dump())
        #chain = model_r | output_parser
        #result = safe_invoke(chain, update_temp)
        result_str = ""
        for chunk in safe_stream_model(model, template):
            if chunk.content:
                print(chunk.content,end="",flush=True)
                result_str += chunk.content
        result = output_parser.invoke(result_str)
        #logger.info(f"[code_evaluator] result: {result}")

    if result["status"] == "CODE_CONFIRMATION":
        new_message = SystemMessage(role="assistant", content=f"thought: evaluate current code")
        updated_state = state.model_copy(update={
            "messages": state.messages + [new_message],
            "code": result["code"],
            "next": "code_evaluator",
            "status": result["status"],
            "agent": "code_evaluator"
        })
        
    else:
        #logger.info(f"[code_evaluator] now execute function")
        state.code = result["code"]
        print(f"[DEBUG] result code: {result['code']}") 
        state.code = extract_python_blocks_tool(state.code)
        print(f"[DEBUG] extracted code: {state.code}")  
        if len(state.checked_step) > 0:
            #logger.info(f"[code_evaluator] step id: {state.checked_step[-1]}")
            sid = str(state.checked_step[-1]) + "_" + datetime.now().strftime("%Y%m%d_%H%M%S")
        else:
            sid = "0_" + datetime.now().strftime("%Y%m%d_%H%M%S")
        result_template = code_executer(step=state.step_output, code=state.code, step_id=sid, keep_files=True, save_dir=state.save_dir)
        code_output = result_template["stdout"]
        state.code_output.append(code_output)
        #logger.info(f"[code_evaluator] code_output: {code_output}")
        
        if state.update_code != "":
            with open('code_up_plan.md', 'r', encoding='utf-8') as f:
                template = f.read().format(**state.model_dump())
            #result = safe_invoke_model.invoke(model_r,template)
            result_str = ""
            for chunk in safe_stream_model(model, template):
                if chunk.content:
                    print(chunk.content,end="",flush=True)
                    result_str += chunk.content
            result = output_parser.invoke(result_str)
            #logger.info(f"[code_evaluator] Update plan: {result}")

        state.functions.append(state.code)
        new_message = SystemMessage(role="assistant", content=f"thought: The python code has been saved, execution: The python code has been saved")
        updated_state = state.model_copy(update={
            "messages": state.messages + [new_message],
            "thought": "The python code has been saved",
            "execution": "The python code has been saved",
            "functions": state.functions,
            "find_function": None,
            "find_tutorial": None,
            "code": None,
            "code_output": state.code_output,
            "next": "plan_executor",
            "update_plan": result,
            "status": "NOT_FINISHED",
            "agent": "code_evaluator"
        })
    return updated_state

def code_executer_old(state: BrickState) -> BrickState:
    
    agent = create_agent(
        "code_executer", 
        "code_executer", 
        [], 
        "code_executer",
        state.model_dump()
    )
    print("code_executer invoke message:", state.messages[-1])
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("code_executer result:", result)
    print("code_executer result type:", type(result))
    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, output: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "next": result.get("next", ""),
            "status": result.get("status", "")
        }
    )

    return updated_state



def general_responder(state: BrickState) -> BrickState:
    
    agent = create_agent(
        "general_responder", 
        "general_responder", 
        [], 
        "general_responder",
        state.model_dump()
    )
    print("general_responder invoke message:")
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("general_responder result:")

    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", ""))
    print("[output]", result.get("output", ""))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, output: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "status": "FINISHED",
            "agent": "general_responder"
        }
    )
    

    return updated_state

def notebook_searcher(state: BrickState) -> BrickState:
    fn = find_notebook(state.question,state.notebooks_path)
    nt = extract_notebook_text(fn)
    
    print("find notebook:", fn)
    updated_state = state.model_copy(
        update={
            "reference_notebook_path": find_notebook(state.question,state.notebooks_path),
            "notebook_text": nt,
            "status": "NOT_FINISHED"
        }
    
    )
    return updated_state

def code_runner(state: BrickState) -> BrickState:
    """代码执行节点 - 通过沙箱运行代码
    
    从 state 中获取 sandbox_id，连接到对应的沙箱，执行 state.code
    """
    print("code_runner received state:")
    step_content = state.step_content
    if state.agent == "code_debugger":
        mdata = "Step"+ str(step_content.get("step")) + ":" + step_content.get("type") + " debugging" + "."
    else:
        mdata = "Step"+ str(step_content.get("step")) + ":" + step_content.get("type") + "."
    thought = state.thought
    code = state.code
    notebook_cells = state.notebook_cells

    # 检查是否有 sandbox_id
    if not state.sandbox_id:
        error_msg = "沙箱未初始化，state 中缺少 sandbox_id"
        print(f"❌ {error_msg}")
        return state.model_copy(
            update={
                "process_flag": 0,
                "error_message": error_msg,
                "next": "code_debugger",
                "status": "ERROR",
                "agent": "code_runner"
            }
        )
    
    # 检查是否有代码
    if not state.code:
        error_msg = "没有要执行的代码，state.code 为空"
        print(f"❌ {error_msg}")
        return state.model_copy(
            update={
                "process_flag": 0,
                "error_message": error_msg,
                "next": "coder",
                "status": "NOT_FINISHED",
                "agent": "code_runner"
            }
        )
    
    try:
        # 1. 重建沙箱客户端，连接到已存在的内核
        print("connecting sandbox")
        box = DockerSandbox(
            server_url="127.0.0.1:8888",
            workspace="/workspace"
        )
        box.kernel_id = state.sandbox_id

        max_retries = 60        # 最多重试 60 次
        retry_interval = 1.0    # 每次间隔 1 秒

        attempt = 0
        while True:
            print(f"🔌 尝试连接沙箱（第 {attempt + 1} 次）...")
            if box._connect_ws():
                break
            attempt += 1
            if attempt >= max_retries:
                raise RuntimeError(f"在重试 {max_retries} 次后仍无法连接沙箱 ID: {state.sandbox_id}")
            time.sleep(retry_interval)

        print(f"✅ 已连接到沙箱 ID: {state.sandbox_id}")
            
        
        # 2. 执行代码
        print(f"🚀 开始执行代码 ({len(state.code)} chars)...")
        result = box.run(state.code, verbose=True)
        
        # 3. 解析执行结果
        stdout_lines = format_stdout(result)
        result_value = format_result(result)
        error_msg = format_error(result)
        is_error = has_error(result)
        notebook_cells.append(
            {
                "mdata": mdata,
                "thought": thought,
                "code": code,
                "outputs": result,
            }
        )
        # 5. 根据是否有错误决定下一步
        if is_error:
            print(f"❌ 代码执行出错")
            print(f"错误信息: {error_msg}")
            
            updated_state = state.model_copy(
                update={
                    "process_flag": 0,
                    "error_message": error_msg,
                    "code_output": stdout_lines,
                    "next": "code_debugger",
                    "status": "NOT_FINISHED",
                    "agent": "code_runner",
                    "complete_output": result,
                    "notebook_cells": notebook_cells
                }
            )
        else:
            print(f"✅ 代码执行成功")
            
            updated_state = state.model_copy(
                update={
                    "process_flag": 1,
                    "code_output": stdout_lines,
                    "next": "plan_executor",
                    "status": "NOT_FINISHED",
                    "agent": "code_runner",
                    "complete_output": result,
                    "error_code": "",
                    "notebook_cells": notebook_cells
                }
            )
        
        return updated_state
        
    except Exception as e:
        error_msg = f"沙箱连接或执行异常: {str(e)}"
        print(f"❌ {error_msg}")
        import traceback
        traceback.print_exc()
        
        return state.model_copy(
            update={
                "process_flag": 0,
                "error_message": error_msg,
                "next": "code_debugger",
                "status": "ERROR",
                "agent": "code_runner"
            }
        )

def query_agent(state: BrickState) -> BrickState:
    print("query_agent received state:")
    ques = state.question
    # 确保 save_dir 是绝对路径
    save_dir = os.path.abspath(state.save_dir) if state.save_dir else None
    print(f"[DEBUG] query_agent save_dir: {save_dir}")
    res = run_query(ques, save_dir)
    res_text = res.get("final_res", "")
    res_table = res.get("relation_frame", "")
    return state.model_copy(
        update={
            "output": res_text,
            "status": "FINISHED",
            "relation_frame": res_table,
            "agent": "responder"
        }
    )