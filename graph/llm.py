from pathlib import Path
from typing import Any, Dict
from langchain_openai import ChatOpenAI
import json

from agent import LLMType
from llm_loader import load_yaml_config, load_env_config
#from src.utils.timer import timed_node
from output import safe_stream_model_str

# Cache for LLM instances
_llm_cache: dict[LLMType, ChatOpenAI] = {}

#@timed_node("_create_llm_use_conf")
def _create_llm_use_conf(llm_type: LLMType, conf: Dict[str, Any]) -> ChatOpenAI:
    llm_type_map = {
        "reasoning": conf.get("REASONING_MODEL"),
        "basic": conf.get("BASIC_MODEL"),
        "vision": conf.get("VISION_MODEL"),
        "embedding": conf.get("EMBEDDING_MODEL"),
        "supervisor": conf.get("SUPERVISOR_MODEL"),
        "coder":conf.get("CODE_MODEL"),
        "claude":conf.get("CLAUDE_MODEL"),
        "qwen_coder":conf.get("QWEN_CODER_MODEL"),
    }
    
    llm_conf = llm_type_map.get(llm_type)
    #print("llm_conf:", llm_conf)
    #print(type(llm_conf))
    if not llm_conf:
        raise ValueError(f"Unknown LLM type: {llm_type}")
    if isinstance(llm_conf, str):
        llm_conf = json.loads(llm_conf)
    if not isinstance(llm_conf, dict):
        raise ValueError(f"Invalid LLM Conf type (not dict): {llm_type}")
    return ChatOpenAI(**llm_conf)


def get_llm_by_type(llm_type: LLMType, config_source: str="env") -> ChatOpenAI:
    """
    Get LLM instance by type. Returns cached instance if available.
    """
    if llm_type in _llm_cache:
        return _llm_cache[llm_type]

    if config_source == "env":
        # Load env config from the current project directory
        conf = load_env_config(str((Path(__file__).parent / "brick_test_config.env").resolve()))
    else:
        # Load YAML config from the current project directory (if present)
        conf = load_yaml_config(str((Path(__file__).parent / "conf.yaml").resolve()))
    #print("conf:", conf)
    llm = _create_llm_use_conf(llm_type, conf)
    _llm_cache[llm_type] = llm
    return llm


# Initialize LLMs for different purposes - now these will be cached
basic_llm = get_llm_by_type("basic", config_source="env")
embedding_model = get_llm_by_type("embedding", config_source="env")
# In the future, we will use reasoning_llm and vl_llm for different purposes
# reasoning_llm = get_llm_by_type("reasoning")
# vl_llm = get_llm_by_type("vision")


if __name__ == "__main__":
    result = safe_stream_model_str(basic_llm, "Hello, who are you")
    print("Agent:", result)
    #print(basic_llm.invoke("Hello, who are you"))
