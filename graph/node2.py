import os
import subprocess

from scanpy import read_10x_h5
from Debugger_simplified import debug
from tools.mock_adata import generate_enhanced_simulated_data_tool
from tools.summarize_data import summarize_data_tool
from tools.valid_h5ad import validate_h5ad_structure_tool
from tools.write_file import write_file_tool
from create_agent import create_agent
from state import BrickState
from langchain_core.messages import SystemMessage, HumanMessage
import json
import tempfile
from datetime import datetime
from pathlib import Path
from tools.extract_code import extract_python_blocks_tool
from tools.read_notebook import read_notebook_tool
from tools.brick_rag_searcher import perform_rag_search_tool
from langchain_core.output_parsers import JsonOutputParser, StrOutputParser
from output import safe_stream_model
from llm import basic_llm as model
from RAG import RAG, RAG_func
from notebook_finder import find_notebook
output_parser = JsonOutputParser()

def env_checker(state: BrickState) -> BrickState:
    print("env_checker received state:", state)
    print("now msg:",state.messages)
    #state.messages.append(SystemMessage(role="user", content=state.question))

    if state.update_data_info:
        state.messages.append(HumanMessage(content=state.update_data_info))

    agent = create_agent(
        "env_checker",
        "env_checker",
        [
            validate_h5ad_structure_tool,
            generate_enhanced_simulated_data_tool,
            summarize_data_tool,
            write_file_tool,
        ],
        "env_checker",
        state.model_dump()
    )
    invoke_message = state.messages[-2:]
    print("env_checker invoke message:",invoke_message)
    print("message list:",state.messages)
    result = agent.invoke({"messages":invoke_message})
    #result = safe_stream_agent_json(agenmessages[-1]t, {"messages":state.messages[-1]})
    
    print("env_checker result:", result)
    print("env_checker result type:", type(result))

    if "messages" in result:
        for msg in result["messages"]:
            print("msg:",msg)
        result = result["messages"][-1].content
    result = json.loads(result)

    if isinstance(result, dict):
        print("[thought]", result["thought"])
        print("[output]", result["output"])
            
        new_message = SystemMessage(content=f"thought: {result['thought']}, output: {result['output']}, data_info: {result['data_info']}, status: {result['status']}")
        updated_state = state.model_copy(
            update={
                "messages": state.messages + [new_message],
                "thought": result["thought"],
                "output": result["output"],
                "next": result["next"],
                "status": result["status"],
                "valid_data": result["valid_data"],
                "data_path": result["data_path"],
                "data_info": result["data_info"],
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

    print("env_checker updated state:", updated_state)
    return updated_state

def supervisor(state: BrickState) -> BrickState:
    print("supervisor received state:", state)
    print("now msg:", state.messages)

    agent = create_agent(
        "supervisor",
        "supervisor",
        [],
        "supervisor",
        state.model_dump()
    )
    
    print("supervisor invoke message:", state.messages[-1])
    result = agent.invoke({"messages": state.messages[-1]})
    print("supervisor result:", result)
    print("supervisor result type:", type(result))
    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        # Fallback: keep raw content so UI can display something
        result = {"output": result, "status": "NOT_FINISHED"}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[thought]", result.get("thought", []))
    print("[output]", result.get("output", []))
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', [])}, output: {result.get('output', [])}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", []),
            "next": result.get("next", []),
            "status": "NOT_FINISHED",
            "data_path": result.get("data_path", ""),
            "output": result.get("output", []),
            "agent": "supervisor"
        }
    )
    
    print("supervisor updated state:", updated_state)
    return updated_state

def translator(state: BrickState) -> BrickState:
    print("translator received state:", state)
    print("now msg:", state.messages)
    
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
        for msg in result["messages"]:
            print("msg:", msg)
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
    
    print("translator updated state:", updated_state)
    return updated_state

def data_analyzer(state: BrickState) -> BrickState:
    print("data_analyzer received state:", state)

    if state.update_data_repo:
        state.messages.append(HumanMessage(content=state.update_data_repo))

    agent = create_agent("data_analyzer", "data_analyzer", [], "data_analyzer", state.model_dump())
    #result = safe_stream_agent_json(agent, state.model_dump())

    print("data_analyzer invoke message:",state.messages[-1])
    result = agent.invoke({"messages":[state.messages[-1]]})
    print("data_analyzer result:", result)

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
                "status": "NOT_FINISHED",
                "next": "supervisor",
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

    print("updated state: ",updated_state)
    #updated_state.messages.append(new_message)
    return updated_state

def analyze_planner(state: BrickState) -> BrickState:
    print("analyze_planner received state:", state)
    print("now msg:", state.messages)
  
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
                "status": result.get("status", ""),
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
                "status": "NOT_FINISHED",
                "next": "supervisor",
                "make_analysis": True,
                "agent": "analyze_planner"
            }
        )
    
    print("analyze_planner updated state:", updated_state)
    return updated_state

def planner(state: BrickState) -> BrickState:
    print("planner received state:", state)
    print("now msg:", state.messages)
    
    agent = create_agent(
        "planner", 
        "planner", 
        [perform_rag_search_tool], 
        "planner",
        state.model_dump()
    )
    print("planner invoke message:", state.messages[-1])
    result = agent.invoke({"messages": state.messages})
    
    print("planner result:", result)
    print("planner result type:", type(result))
    
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
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, current_plan: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "current_plan": result.get("output", ""),
            "next": result.get("next", ""),
            "status": result.get("status", ""),
            "make_plan": True,
            "agent": "planner"
        }
    )
    
    print("planner updated state:", updated_state)
    return updated_state

def code_debugger(state: BrickState) -> BrickState:
    print("code_debugger received state:", state)
    print("now msg:", state.messages)
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
    
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, full_code: {result.get('output', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "thought": result.get("thought", ""),
            "output": result.get("output", ""),
            "full_code": result.get("output", ""),
            "next": result.get("next", ""),
            "status": result.get("status", ""),
            "agent": "code_debugger"
        }
    )
    
    print("code_debugger updated state:", updated_state)
    return updated_state


def plan_reviewer(state: BrickState) -> BrickState:
    print("plan_reviewer received state:", state)
    print("now msg:", state.messages)
    
    if state.update_instruction:
        state.messages.append(HumanMessage(content=state.update_instruction[-1]))
    
    agent = create_agent(
        "plan_reviewer", 
        "plan_reviewer", 
        [], 
        "plan_reviewer",
        state.model_dump()
    )
    print("plan_reviewer invoke message:", state.messages[-1])
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("plan_reviewer result:", result)
    print("plan_reviewer result type:", type(result))
    
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
    
    print("plan_reviewer updated state:", updated_state)
    return updated_state

def responder(state: BrickState) -> BrickState:
    print("responder received state:", state)
    print("now msg:", state.messages)
    
    agent = create_agent(
        "responder", 
        "responder", 
        [], 
        "responder",
        state.model_dump()
    )
    print("responder invoke message:", state.messages[-1])
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("responder result:", result)
    print("responder result type:", type(result))
    
    if "messages" in result:
        result = result["messages"][-1].content
    
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {}
    
    if not isinstance(result, dict):
        result = {}
    
    print("[answer]", result.get("answer", ""))
    
    new_message = SystemMessage(content=f"final_answer: {result.get('answer', '')}")
    updated_state = state.model_copy(
        update={
            "messages": state.messages + [new_message],
            "final_answer": result.get("answer", ""),
            "status": "FINISHED",
            "agent": "responder"
        }
    )
    
    print("responder updated state:", updated_state)
    return updated_state

def plan_executor(state: BrickState) -> BrickState:
    print("plan_executor received state:", state)
    print("now msg:", state.messages)
    
    agent = create_agent(
        "plan_executor", 
        "plan_executor", 
        [], 
        "plan_executor",
        state.model_dump()
    )
    print("plan_executor invoke message:", state.messages[-1])
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("plan_executor result:", result)
    print("plan_executor result type:", type(result))


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
    rag_res = RAG_func(query=json.dumps(result.get("step_content", "")),vectorstore=state.default_vectorstore)
    print("== rag_res ==", rag_res)
    new_message = SystemMessage(content=f"thought: {result.get('thought', '')}, output: {result.get('output', '')}")
    updated_checked_step = result.get("checked_step", state.checked_step)

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
    
    print("plan_executor updated state:", updated_state)
    return updated_state

def searcher(state: BrickState) -> BrickState:
    print("searcher received state:", state)
    print("now msg:", state.messages)
    
    agent = create_agent(
        "searcher", 
        "searcher", 
        [], 
        "searcher",
        state.model_dump()
    )
    print("searcher invoke message:", state.messages[-1])
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
    
    print("searcher updated state:", updated_state)
    return updated_state

def bkg_searcher(state: BrickState) -> BrickState:
    print("bkg_searcher received state:", state)
    print("now msg:", state.messages)
    
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
    
    print("bkg_searcher updated state:", updated_state)
    return updated_state

def coder(state: BrickState) -> BrickState:
    print("coder received state:", state)
    print("now msg:", state.messages)
    
    agent = create_agent(
        "coder", 
        "coder", 
        [extract_python_blocks_tool], 
        "coder",
        state.model_dump()
    )
    invoke_message = state.messages
    print("coder invoke message:", invoke_message)
    result = agent.invoke({"messages": invoke_message})
    
    print("coder result:", result)
    print("coder result type:", type(result))
    
    if "messages" in result:
        for msg in result["messages"]:
            print("msg:",msg)
        result = result["messages"][-1].content
    # Robust JSON parsing with fallback to raw text
    try:
        result = json.loads(result)
    except json.JSONDecodeError:
        result = {"thought": "", "output": str(result), "status": "CODE_CONFIRMATION", "next": "code_runner"}
    if isinstance(result, dict):
        print("[thought]", result.get("thought", ""))
        print("[output]", result.get("output", ""))

        new_message = SystemMessage(role="assistant", content=f"thought: {result['thought']}, code: {result['output']}")
        updated_state = state.model_copy(update={
            "messages": state.messages + [new_message],
            "thought": result["thought"],
            "output": result.get("output", ""),
            "full_code": state.full_code + "\n" + result["output"],
            "code": result["output"],
            "next": "code_runner",
            "status": "CODE_CONFIRMATION",
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
    
    print("coder updated state:", updated_state)
    return updated_state

def code_controller(state: BrickState) -> BrickState:
    print("code_controller received state:", state)
    print("now msg:", state.messages)
    
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
    
    print("code_controller updated state:", updated_state)
    return updated_state

def code_evaluator_old(state: BrickState) -> BrickState:
    print("code_evaluator received state:", state)
    print("now msg:", state.messages)
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
    
    print("code_evaluator updated state:", updated_state)
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
    print("code_executer received state:", state)
    print("now msg:", state.messages)
    
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
    
    print("code_executer updated state:", updated_state)
    return updated_state

def code_executer(step: str, code: str, step_id: str, keep_files: bool = False, save_dir: str = None) -> dict:
    """
    修正后的执行函数（解决文件保存问题）
    """
    result_template = {
        "stdout": "", 
        "stderr": "",
        "returncode": -1,
        "file_path": None,
        "error": None
    }
    
    script_path = None  # 显式初始化
    use_tempfile = False  # 默认值
    
    try:
         # 检查并创建 plots 目录
        plots_dir = "plots"
        Path(plots_dir).mkdir(parents=True, exist_ok=True)
        
        file_id = f"BRICK_{step_id}"
        if save_dir:
            save_dir = os.path.expanduser(save_dir)
            Path(save_dir).mkdir(parents=True, exist_ok=True)
            script_path = Path(save_dir) / f"generated_{file_id}.py"
            use_tempfile = False
        else:
            temp_dir = tempfile.gettempdir()
            script_path = Path(temp_dir) / f"generated_{file_id}.py"
            use_tempfile = True
        
        fixed_code = """
import sys
import os
os.chdir("/home/lyt/checker_finallap") # your working path
sys.path.append("/home/lyt/checker_finallap")

import scanpy as sc
from scipy.sparse.csgraph import connected_components
import numpy as np
import scipy.sparse
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import itertools
from datetime import datetime
from pathlib import Path

import BRICK

#base_url = ""
#api_key = ""

api_key = "sk-kpsteSkpDGl1xBmDEcC7D51b968e43499092826f17286b55"
base_url = "http://10.224.28.80:3000/v1"
llm_params = {"model_name": "deepseek-v3-hs", "temperature": 0.7}

BRICK.config_llm(modeltype='ChatOpenAI', base_url=base_url, api_key=api_key, llm_params=llm_params)

url = "neo4j://10.224.28.66:7687"
auth = ("neo4j", "bmVvNGpwYXNzd29yZA==")  

BRICK.config(url=url, auth=auth)

plots_dir = "plots"
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir, exist_ok=True)  # 自动创建目录
    print(f"创建目录: {{plots_dir}}")
else:
    print(f"目录已存在: {{plots_dir}}")

current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
"""
        full_script = f"""
{fixed_code}

# ----------User code start----------
{code}
# ----------User code end------------
if __name__ == "__main__":
    try:
        print("EXECUTION_START")
        exec('''{code}''')
    except Exception as e:
        print(f"EXECUTION_ERROR:{{str(e)}}")
        raise
"""
        """"""
        # 调试输出
        print(f"生成脚本路径: {script_path}")  # 添加调试信息
        
        # 写入文件（处理权限）
        try:
            with open(script_path, "w", encoding="utf-8") as f:  # 指定编码
                f.write(full_script)
            script_path.chmod(0o600)
            print(f"文件写入成功: {script_path.stat().st_size} 字节")  # 验证写入
            #print("生成的脚本内容:\n", full_script)
            print("Now running code...")
            debug_state, code, output = debug(step,full_script,work_dir=save_dir)
            print("debug_state:", debug_state)
            #print('output: ', output)
            result_template["stdout"] = output
            return result_template
        except Exception as e:
            result_template["error"] = f"文件写入失败: {str(e)}"
            return result_template
        
    except subprocess.TimeoutExpired as e:
        result_template["error"] = f"执行超时: {str(e)}"
        return result_template
    except Exception as e:
        result_template["error"] = f"系统错误: {str(e)}"
        return result_template
    finally:
        if script_path and script_path.exists():
            should_delete = (use_tempfile and not keep_files) or (not use_tempfile and not keep_files)
            if should_delete:
                try:
                    script_path.unlink()
                except Exception as e:
                    pass

def general_responder(state: BrickState) -> BrickState:
    print("general_responder received state:", state)
    print("now msg:", state.messages)
    
    agent = create_agent(
        "general_responder", 
        "general_responder", 
        [], 
        "general_responder",
        state.model_dump()
    )
    print("general_responder invoke message:", state.messages[-1])
    result = agent.invoke({"messages": [state.messages[-1]]})
    
    print("general_responder result:", result)
    print("general_responder result type:", type(result))
    
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
            "agent": "general_responder"
        }
    )
    
    print("general_responder updated state:", updated_state)
    return updated_state

def notebook_searcher(state: BrickState) -> BrickState:
    print("find notebook:", find_notebook(state.question,state.notebooks_path))
    updated_state = state.model_copy(
        update={
            "reference_notebook_path": find_notebook(state.question,state.notebooks_path)
        }
    
    )
    return updated_state

def code_runner(state: BrickState) -> BrickState:
    success, code, output = debug(state.step,state.full_code,work_dir=state.save_dir)
    if success:
        updated_state = state.model_copy(
            update={
                "code_output": state.code_output+[output],
                "next": "plan_executor",
                "status": "VALIDATED",
                "agent": "code_runner"
            }
        )
    else:
        updated_state = state.model_copy(
            update={
                "error_message": output,
                "next": "code_debugger",
                "status": "NOT_FINISHED",
                "agent": "code_runner"
            }
        )
    return updated_state