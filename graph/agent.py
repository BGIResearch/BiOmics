from typing import Literal

# Define available LLM types
LLMType = Literal["basic", "reasoning", "vision", "embedding","supervisor"]

# Define agent-LLM mapping
AGENT_LLM_MAP: dict[str, LLMType] = {
    "env_checker": "claude",
    "translator": "basic",
    "data_analyzer": "basic",
    "analyze_planner": "basic",
    "planner": "basic",
    "searcher": "basic",
    "plan_reviewer": "basic",
    "plan_executor": "basic",
    "coder": "claude",
    "code_debugger": "claude",
    "code_evaluator": "basic",
    "responder": "basic",
    "general_responder": "basic",
    "supervisor": "claude"
}



