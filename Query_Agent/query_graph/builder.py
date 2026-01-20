from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.memory import MemorySaver

from .state import QueryState
from .nodes import (
    translator,
    optant,
    extractor,
    hybrid_searcher,
    reranker,
    executor,
    general_responder
)

def route_next_step(state: QueryState) -> str:
    if state.chosen_script == "1" or state.chosen_script == 1 or state.chosen_script == "2" or state.chosen_script == 2 or state.chosen_script == "3" or state.chosen_script == 3:
        return "extractor"
    else:
        return "general_responder"

def _build_base_graph():
    """Build and return the base state graph with all nodes and edges."""
    builder = StateGraph(QueryState)
    builder.add_node("translator", translator)
    builder.add_node("optant", optant)
    builder.add_node("extractor", extractor)
    builder.add_node("hybrid_searcher", hybrid_searcher)
    builder.add_node("reranker", reranker)
    builder.add_node("executor", executor)
    builder.add_node("general_responder", general_responder)

    builder.add_edge(START, "translator")
    builder.add_edge("translator", "optant")
    builder.add_conditional_edges(
        "optant",
        route_next_step,
        {
            "extractor": "extractor", 
            "general_responder": "general_responder"
        }
    )
    builder.add_edge("extractor", "hybrid_searcher")
    builder.add_edge("hybrid_searcher", "reranker")
    builder.add_edge("reranker", "executor")
    builder.add_edge("executor", END)
    builder.add_edge("general_responder", END)
    return builder

def build_graph_with_memory():
    """Build and return the agent workflow graph with memory."""
    memory = MemorySaver()
    builder = _build_base_graph()
    return builder.compile(checkpointer=memory) 

def build_graph_with_interaction(interrupt_before: list[str], interrupt_after: list[str]):
    """Build and return the agent workflow graph with memory and human in loop."""
    memory = MemorySaver()
    builder = _build_base_graph()
    return builder.compile(checkpointer=memory,interrupt_before=interrupt_before,interrupt_after=interrupt_after)


def build_graph():
    """Build and return the agent workflow graph without memory and human in loop."""
    builder = _build_base_graph()
    return builder.compile()

#graph = build_graph()

if __name__ == "__main__":
    graph = build_graph()
    print(graph.get_graph(xray=True).draw_mermaid())