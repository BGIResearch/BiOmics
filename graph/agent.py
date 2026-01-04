from typing import Literal

# Define available LLM types
LLMType = Literal["basic", "reasoning", "vision", "embedding","supervisor","claude","coder"]

# Define agent-LLM mapping
AGENT_LLM_MAP: dict[str, LLMType] = {
    "env_checker": "basic",
    #"env_checker": "basic",
    "translator": "basic",
    "data_analyzer": "reasoning",
    "analyze_planner": "reasoning",
    "planner": "basic",
    "searcher": "basic",
    "plan_reviewer": "basic",
    "plan_executor": "basic",
    "coder": "coder",
    #"coder": "basic",
    "code_debugger": "claude",
    #"code_debugger": "basic",
    "code_evaluator": "basic",
    "responder": "reasoning",
    "general_responder": "basic",
    "supervisor": "reasoning",
    #"supervisor": "basic"
}



