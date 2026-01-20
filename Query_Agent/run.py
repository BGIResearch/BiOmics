from langchain_core.messages import HumanMessage

from .query_graph.state import QueryState
from .query_graph.builder import build_graph

def run_query(query: str, table_dir: str):
    state_data = {"source_question": query, "table_dir": table_dir}
    initial_state = QueryState(**state_data)
    graph = build_graph()
    initial_state_dict = initial_state.model_dump()
    final_result = graph.invoke(initial_state_dict)


    return final_result

if __name__ == "__main__":
    question = input("请输入查询问题:")
    result = run_query(question,"./tables")
    print("final result: ",result)