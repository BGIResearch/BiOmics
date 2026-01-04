import os
import time
from pathlib import Path
from dotenv import load_dotenv
from graph.llm import get_llm_by_type
from langchain_core.messages import SystemMessage, HumanMessage
from langchain_core.prompts import ChatPromptTemplate

# 加载环境变量
config_file = Path(__file__).parent.parent / 'graph' / 'brick_test_config.env'
load_dotenv(dotenv_path=str(config_file))
PROJECT_ROOT = os.getenv('PROJECT_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def summarize_question(question: str) -> str:
    model = get_llm_by_type("basic")
    template = """
<role>
You are a professional summary generation assistant that only returns highly summarized problem summaries.
</role>

<rules>
1. Answer in the same language as the question
2. The Chinese abstract should not exceed five characters
3. Abstracts in other languages should not exceed five words
4. Use "_" to separate words
5. Only the summary content is returned, without any other text
</rules>

<examples>
Input: "I want to use BRICK for trajectory inference tasks."
Output: "Trajectory_Inference"

Input: "How to calculate 3 + 2"
Output: "Add_calculation"
</examples>
    """
    #template = ChatPromptTemplate.from_template(template.format(question=question))
    #print("system prompt: ",[SystemMessage(content=template)])
    #print("huamn: ",[HumanMessage(content=question)])
    result = model.invoke([SystemMessage(content=template)]+[HumanMessage(content=question)])
    #print("new question:", result)
    if result.content: 
        output = result.content
    return output 

def get_save_dir(question: str, save_path: str = None) -> str:
    if save_path is None:
        save_path = os.path.join(PROJECT_ROOT, "data", "run_file")
    summary = summarize_question(question)
    timestamp = time.strftime("%a_%b_%d_%Y_%H_%M_%S", time.localtime())
    folder_name = f"{summary}_{timestamp}"
    save_dir = os.path.join(save_path, folder_name)
    os.makedirs(save_dir, exist_ok=True)
    return save_dir

if __name__ == "__main__":
    question = "What is the capital of France?"
    save_dir = get_save_dir(question)
    print(f"Save directory: {save_dir}")