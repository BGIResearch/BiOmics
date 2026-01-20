from langgraph.graph import StateGraph, START, END
from langgraph.checkpoint.memory import MemorySaver

from .state import BrickState
from .node2 import (
    env_checker,
    supervisor,
    translator,
    data_analyzer,
    analyze_planner,
    planner,
    plan_reviewer,
    searcher,
    plan_executor,
    coder,
    code_runner,
    code_debugger,
    code_evaluator,
    responder,
    general_responder,
    notebook_searcher,
    query_agent
)

def route_next_step(state: BrickState) -> str:
    routing_map = {
        "env_checker": "env_checker",
        "supervisor": "supervisor",
        "translator": "translator",
        "planner": "planner",
        "plan_reviewer": "plan_reviewer",
        "plan_checker": "plan_checker",
        "plan_executor": "plan_executor",
        "searcher": "searcher",
        "coder": "coder",
        "code_runner": "code_runner",
        "code_debugger": "code_debugger",
        "code_evaluator": "code_evaluator",
        "responder": "responder",
        "data_analyzer": "data_analyzer",
        "analyze_planner": "analyze_planner",
        "general_responder": "general_responder",
        "notebook_searcher": "notebook_searcher",
        "query_agent":"query_agent",
        "END": "END"
    }
    if state.status == "FINISHED":
        return "END"
    elif state.next in routing_map:
        print("choose next step:", state.next)
        return routing_map[state.next]
    else:
        raise ValueError(f"Invalid next step: {state.next}")


def _build_base_graph():
    """Build and return the base state graph with all nodes and edges."""
    builder = StateGraph(BrickState)
    builder.add_node("supervisor", supervisor)
    builder.add_node("env_checker", env_checker)
    builder.add_node("general_responder", general_responder)
    builder.add_node("translator", translator)
    builder.add_node("data_analyzer", data_analyzer)
    builder.add_node("analyze_planner", analyze_planner)
    builder.add_node("planner", planner)
    builder.add_node("plan_reviewer", plan_reviewer)
    builder.add_node("plan_executor", plan_executor)
    builder.add_node("coder", coder) 
    builder.add_node("code_runner", code_runner)
    builder.add_node("code_debugger", code_debugger)
    builder.add_node("responder", responder) 
    builder.add_node("notebook_searcher", notebook_searcher)
    builder.add_node("query_agent", query_agent)
    
    builder.add_edge(START, "notebook_searcher")
    builder.add_edge("notebook_searcher", "supervisor") 
    builder.add_edge("code_debugger", "code_runner")
    builder.add_edge("responder", END)
    builder.add_edge("general_responder", END)
    builder.add_edge("query_agent", END)
    
    builder.add_conditional_edges("translator",route_next_step,{"supervisor": "supervisor"})
    builder.add_edge("coder", "code_runner")
    builder.add_edge("planner", "plan_executor")

    builder.add_conditional_edges(
        "env_checker",
        route_next_step,
        {
            "env_checker": "env_checker", 
            "data_analyzer": "data_analyzer"
        }
    )

    builder.add_conditional_edges(
        "supervisor",
        route_next_step,
        {
            "env_checker": "env_checker",
            "query_agent": "query_agent",
            "general_responder": "general_responder"
        }
    )

    builder.add_conditional_edges(
        "data_analyzer",
        route_next_step,
        {
            "data_analyzer": "data_analyzer", 
            "analyze_planner": "analyze_planner"
        }
    )

    builder.add_conditional_edges(
        "analyze_planner",
        route_next_step,
        {
            "analyze_planner": "analyze_planner", 
            "planner": "planner"
        }
    )

    
    builder.add_conditional_edges(
        "code_runner",
        route_next_step,
        {
            "code_debugger": "code_debugger",
            "plan_executor": "plan_executor"
        }
    )

    builder.add_conditional_edges(
        "plan_executor",
        route_next_step,
        {
            "coder": "coder",
            "responder": "responder"  
        }
    )
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
    graph = build_graph_with_interaction(interrupt_after=["env_checker","data_analyzer"])
    print(graph.get_graph(xray=True).draw_mermaid())