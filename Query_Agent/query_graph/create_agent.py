from langgraph.prebuilt import create_react_agent
from langchain_core.runnables import RunnableLambda
from langchain_core.prompts import ChatPromptTemplate

from .state import QueryState
from .llm import get_llm_by_type
from .agent import AGENT_LLM_MAP
from .load_template import apply_prompt_template

from langchain_core.messages import AIMessage
from typing import Dict, Any
from langchain_core.output_parsers import PydanticOutputParser

def create_agent(agent_name: str, agent_type: str, tools: list, prompt_template_name: str, state_dict: Dict[str, Any]):
    """Factory function to create agents with consistent configuration."""

    #output_parser = PydanticOutputParser(pydantic_object=EnvCheckerOutput)
    messages = apply_prompt_template(prompt_template_name, state_dict)
    # Convert list of messages to ChatPromptTemplate
    prompt = ChatPromptTemplate.from_messages(messages)
    #print("prompt: ",prompt)
    react_agent_core = create_react_agent(
        name=agent_name,
        model=get_llm_by_type(AGENT_LLM_MAP[agent_type]),
        tools=tools,
        #prompt=lambda x: (print("Prompt state:", x) or x.get("messages", [])),
        prompt=prompt
    )
    return react_agent_core
