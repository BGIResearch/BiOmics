from langgraph.prebuilt import create_react_agent
from langchain_core.runnables import RunnableLambda

from state import BrickState
from llm import get_llm_by_type
from agent import AGENT_LLM_MAP
from load_template import apply_prompt_template

""" def apply_prompt_template(template_name, state):
    from load_template import apply_prompt_template as _apply
    print("my state: ",state)
    return _apply(template_name, state) """



from langchain_core.messages import AIMessage
from typing import Dict, Any
from langchain_core.output_parsers import PydanticOutputParser

def create_agent(agent_name: str, agent_type: str, tools: list, prompt_template_name: str, state_dict: Dict[str, Any]):
    """Factory function to create agents with consistent configuration."""

    #output_parser = PydanticOutputParser(pydantic_object=EnvCheckerOutput)
    prompt = apply_prompt_template(prompt_template_name, state_dict)
    #print("prompt: ",prompt)
    react_agent_core = create_react_agent(
        name=agent_name,
        model=get_llm_by_type(AGENT_LLM_MAP[agent_type]),
        tools=tools,
        #prompt=lambda x: (print("Prompt state:", x) or x.get("messages", [])),
        prompt=prompt
    )
    return react_agent_core
'''
    def _prepare_agent_input(state_dict: Dict[str, Any]) -> Dict[str, Any]:
        """
        这个函数接收完整的 AgentState 字典，并将其转换为 react_agent_core 所需的输入格式。
        它会调用 apply_prompt_template 来生成包含 System Prompt 的消息列表。
        """
        print("RunnableLambda input_dict (from _prepare_agent_input):", state_dict) # 打印 LangGraph 传递的完整状态
        
        rendered_messages_list = apply_prompt_template(prompt_template_name, state_dict)
        print("apply_prompt_template return value (rendered messages list):", rendered_messages_list) # 打印 apply_prompt_template 的返回值

        rendered_messages_list = apply_prompt_template(prompt_template_name, state_dict)


        return {
            "messages": rendered_messages_list,
            **state_dict
        }

    def _extract_final_ai_message_content(agent_output: Dict[str, Any]) -> str:
        """
        从 react_agent_core 的输出字典中提取最后一个 AIMessage 的内容。
        这个内容应该是 LLM 按照 Pydantic 格式指令生成的 JSON 字符串。
        """
        print("Extracting final AI message content from agent_output:", agent_output)
        messages = agent_output.get('messages', [])
        
        # 遍历消息列表，找到最后一个 AIMessage
        for msg in reversed(messages):
            if isinstance(msg, AIMessage):
                print("return msg content: ",msg.content)
                return msg.content
        
        print("Warning: No AIMessage found in agent_output for content extraction.")
        return "" 

    # === 构建完整的 Runnable 链 ===
    full_agent_runnable = (
        RunnableLambda(_prepare_agent_input) # 1. 准备 LLM 输入
        | react_agent_core                    # 2. 运行 ReAct 代理，输出一个复杂字典 (包含 messages, tool_calls, intermediate_steps等)
        #| RunnableLambda(_extract_final_ai_message_content) # 3. 从复杂字典中提取最终的 LLM 文本内容（期望是 JSON 字符串）
        #| output_parser             # 4. PydanticOutputParser 接收字符串并尝试解析为 EnvCheckerOutput
    )

    return full_agent_runnable'''

def create_agent3(agent_name: str, agent_type: str, tools: list, prompt_template_name: str):
    """Factory function to create agents with consistent configuration."""

    react_agent_core = create_react_agent(
        name=agent_name,
        model=get_llm_by_type(AGENT_LLM_MAP[agent_type]),
        tools=tools,
        prompt=lambda x: (print("prompt state:",x) or x.get("messages", [])), 
    )

    full_agent_runnable = (
        RunnableLambda(lambda state_dict: (
            print("RunnableLambda input_dict:", state_dict) or
                {
                    "messages": apply_prompt_template(prompt_template_name, state_dict),
                    **state_dict 
                }
            )
        )
        | react_agent_core 
    )

    return full_agent_runnable

def create_agent2(agent_name: str, agent_type: str, tools: list, prompt_template: str):
    """Factory function to create agents with consistent configuration."""
    #print(get_llm_by_type(AGENT_LLM_MAP[agent_type]))
    print("create agent", agent_name)
    agent = create_react_agent(
        name=agent_name,
        model=get_llm_by_type(AGENT_LLM_MAP[agent_type]),
        tools=tools,
        #prompt=lambda state: apply_prompt_template(prompt_template, state),
        prompt=lambda state: (
            print("Prompt state:", state) or
            apply_prompt_template(prompt_template, state)
        ),
    )
    #return RunnableLambda(lambda input_dict: {**input_dict, "messages": input_dict.get("messages", []), "question": input_dict.get("question", ""), "data_info": input_dict.get("data_info", "")}) | agent
    return RunnableLambda(
        lambda input_dict: (
                print("RunnableLambda input_dict:", input_dict) or
                {**input_dict,
                "messages": input_dict.get("messages", []),
                "question": input_dict.get("question", ""),
                "data_info": input_dict.get("data_info", "")}
            )
    ) | agent

'''
def extract_state_from_messages(messages):
    merged_state = {}
    for msg in messages:
        # 从消息的 additional_kwargs 中查找 state_data
        state_data = None
        if hasattr(msg, "additional_kwargs") and isinstance(msg.additional_kwargs, dict):
            state_data = msg.additional_kwargs.get("state_data", None)
        if state_data:
            if isinstance(state_data, dict):
                merged_state.update(state_data)
            elif isinstance(state_data, AgentState):
                merged_state.update(state_data.dict())
    return merged_state
    


def create_agent(agent_name: str, agent_type: str, tools: list, prompt_template: str):
    print("create agent", agent_name)

    def custom_prompt(state):
        # state这里是消息集合对象或类似结构
        print("Prompt received state:", state)

        messages = getattr(state, "messages", [])
        merged_state_vars = extract_state_from_messages(messages)
        print("merged state vars", merged_state_vars)

        base_state_dict = {}
        if hasattr(state, "dict"):
            base_state_dict = state.dict(exclude={"messages"}, exclude_none=True)
        elif isinstance(state, dict):
            base_state_dict = {k: v for k, v in state.items() if k != "messages"}

        # merged_state_vars 优先覆盖 base_state_dict
        print("base state dict", base_state_dict)
        prompt_vars = {**base_state_dict, **merged_state_vars, "messages": messages}
        print("prompt vars", prompt_vars)
        return apply_prompt_template(prompt_template, state)

    agent = create_react_agent(
        name=agent_name,
        model=get_llm_by_type(AGENT_LLM_MAP[agent_type]),
        tools=tools,
        prompt=custom_prompt,
    )

    def runnable_lambda_fn(input_dict):
        print("RunnableLambda input_dict:", input_dict)
        return {
            **input_dict,
            "messages": input_dict.get("messages", []),
            "question": input_dict.get("question", ""),
            "data_info": input_dict.get("data_info", ""),
        }

    return RunnableLambda(runnable_lambda_fn) | agent
'''

'''
from langgraph.prebuilt import create_react_agent
from langchain_core.runnables import RunnableLambda
from state import GraphState
from llm import get_llm_by_type
from agent import AGENT_LLM_MAP
from load_template import apply_prompt_template
from langchain_core.messages import BaseMessage

def create_agent(agent_name: str, agent_type: str, tools: list, prompt_template: str):
    print("Creating agent:", agent_name)

    
    def custom_prompt(state):
        print("Prompt received GraphState:", state)

        if isinstance(state, dict):
        # 👇 将 messages 中的 BaseMessage 对象转为 dict
            if "messages" in state and isinstance(state["messages"], list):
                new_messages = []
                for msg in state["messages"]:
                    if isinstance(msg, BaseMessage):
                        new_messages.append(msg.model_dump())  # 或者 msg.model_dump() 取决于你用的是哪个库
                    elif isinstance(msg, dict):
                        new_messages.append(msg)
                    else:
                        raise ValueError(f"Unsupported message type: {type(msg)}")
                state["messages"] = new_messages
            state = GraphState(**state)

        # 提取扁平化字段用于 prompt 渲染
        prompt_vars = state.model_dump(exclude_none=True)
        print("Prompt vars (flattened state):", prompt_vars)

        return apply_prompt_template(prompt_template, prompt_vars)


    agent = create_react_agent(
        name=agent_name,
        model=get_llm_by_type(AGENT_LLM_MAP[agent_type]),
        tools=tools,
        prompt=custom_prompt,
    )

    # 简单 passthrough：GraphState -> GraphState dict（用于后续节点处理）
    def passthrough(state_dict: dict):
        print("RunnableLambda passthrough input:", state_dict)
        return state_dict

    return RunnableLambda(passthrough) | agent
'''