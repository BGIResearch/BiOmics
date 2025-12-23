import os
import dataclasses
from datetime import datetime
from typing import Union
from jinja2 import Environment, FileSystemLoader, select_autoescape
from langchain_core.messages import SystemMessage, BaseMessage

from configuration import Configuration
from state import BrickState

# Initialize Jinja2 environment
env = Environment(
    loader=FileSystemLoader(os.path.dirname(__file__)),
    autoescape=select_autoescape(),
    trim_blocks=True,
    lstrip_blocks=True,
)


def apply_prompt_template(
    prompt_name: str,
    state: Union[BrickState, dict],
    configurable: Configuration = None
) -> list[BaseMessage]:
    """
    Apply template variables to a prompt template and return formatted messages.

    Args:
        prompt_name: Name of the template file (without `.md` extension)
        state: A GraphState object or dictionary representing the agent state
        configurable: Optional Configuration object to inject into the prompt

    Returns:
        A list of LangChain messages (starting with a SystemMessage)
    """
    if isinstance(state, dict):
        state_vars = state
    elif hasattr(state, "model_dump"):
        state_vars = state.model_dump(exclude_none=True)
    else:
        raise ValueError("Invalid state format: must be BrickState or dict")

    print("init state vars:",state_vars)
    # 添加通用变量
    state_vars["CURRENT_TIME"] = datetime.now().strftime("%a %b %d %Y %H:%M:%S %z")

    if configurable:
        state_vars.update(dataclasses.asdict(configurable))

    try:
        template = env.get_template(f"{prompt_name}.md")
        #print("apply_prompt_template vars:", state_vars)
        system_prompt = template.render(**state_vars)
        #print("system prompt:",system_prompt)
        #return [SystemMessage(content=system_prompt)]
        message = state_vars.get("messages", [])[-1]
        return system_prompt + str(message)
    except Exception as e:
        raise ValueError(f"Error applying template {prompt_name}: {e}")



