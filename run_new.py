from langchain.schema import HumanMessage

from graph.state import BrickState
from graph.builder import build_graph_with_interaction
from utils.save_dir_namesave_dir_name import get_save_dir

def test(
    question: str,
    file_path: str,
    config: dict,
    preview_n: int,
    save_dir: str = None,
):

    state_data = {
        "question": question,
        "messages": [HumanMessage(content=question)],
    }
    if file_path:
        state_data["data_path"] = file_path
    if preview_n:
        state_data["preview_n"] = preview_n
    if save_dir:
        state_data["save_dir"] = save_dir

    initial_state = BrickState(**state_data)
    print("init state: ",initial_state)
    print(type(initial_state))
    graph = build_graph_with_interaction(interrupt_after=["env_checker","data_analyzer"],interrupt_before=[])
    initial_state_dict = initial_state.model_dump()
    events = graph.stream(initial_state_dict, config=config, stream_mode="values")

    for event in events:
        status = event.get("status")
        print("event: ",event)        
    try:
        while status not in ["FINISHED", "Confirm"]:
            if status == "AWAITING_CONFIRMATION":
                thought = event.get('thought', 'Did you forget to upload your data?')
                output = event.get('output', 'Did you forget to upload your data?')
                update_data_info = input(f"System thought: {thought},\n System output:{output} \n")
                graph.update_state(config=config,values={"update_data_info": update_data_info})
                events = graph.stream(None, config, stream_mode="values")

            elif status == "Revise":
                thought = event.get('thought', 'Did you forget to upload your data?')
                output = event.get('data_repo', 'Can not find data_repo')
                update_data_repo = input(f"System thought: {thought},\n System output:{output} \n")
                graph.update_state(config=config,values={"update_data_repo": update_data_repo})
                print("update_data_repo: ",update_data_repo)
                events = graph.stream(None, config, stream_mode="values")

            elif status in ["VALIDATED","NOT_FINISHED"]:
                events = graph.stream(None, config, stream_mode="values")
            
            for event in events:
                status = event.get("status")
                print("status: ",status,type(status))
        
        print("quit while loop")
        if status in ["FINISHED", "Confirm"]:
            final_result = event.get("final_result","Oops! Something went wrong.")

    except StopIteration:
        pass
    
    return final_result

if __name__ == "__main__":
    config = {"configurable": {"thread_id": "22"}, "recursion_limit": 200}
    question = "我想用BRICK进行轨迹推断任务，/Users/shang/Desktop/华大/test/enhanced_simulated.h5ad"
    save_dir = get_save_dir(question)
    result = test(question, file_path="/Users/shang/Desktop/华大/test/enhanced_simulated.h5ad", config=config, preview_n=20, save_dir=save_dir)
    print("final result: ",result)