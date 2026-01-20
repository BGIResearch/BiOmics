from typing import Literal

# Define available LLM types
LLMType = Literal["basic","claude","embedding","reasoning"]

# Define agent-LLM mapping
AGENT_LLM_MAP: dict[str, LLMType] = {
    "translator": "reasoning",
    "optant": "reasoning",
    "extractor": "reasoning",
    "reranker": "reasoning",
    "general_responder": "reasoning"
}



