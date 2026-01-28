import os
from pydantic import BaseModel, Field, model_validator
from langchain_core.messages import AnyMessage, BaseMessage
from langgraph.graph.message import add_messages
from datetime import datetime
from typing import Annotated, Literal, Union, List, Dict, Any, Optional

from pynndescent.pynndescent_ import process_candidates

# 获取项目根目录
PROJECT_ROOT = os.getenv('PROJECT_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

class BrickState(BaseModel):
    # 分析LLM生成的plan
    a_plan: Optional[Union[str, dict, list]] = None
    # AGENT
    agent: Optional[Union[str, dict, list]] = None
    # BRICK说明
    brick_info: list = []
    # BRICK说明路径
    brick_info_path: Optional[str] = None
    # 用户对步骤的更新
    change_step: Optional[str] = None
    # 环境报告
    check_md: str = ""
    # checked step list
    checked_step: list[Union[str, dict, list, int]] = []
    # 合并后的code
    code: Optional[str] = None
    # code运行输出
    code_output: list[Union[str, dict, list]] = []
    # 完整的代码运行结果
    complete_output: dict = {}
    # 用户id
    customer_id: str = ""
    # 生成的cypher
    cypher: Optional[str] = None
    # planner生成的当前的plan
    current_plan: Optional[Union[str, dict, list]] = None
    # 当前第几步
    current_step: int = 0
    # 数据的总览
    data_info: Optional[Union[str, Dict[str, Union[str, int, float, list, dict]]]] = None
    # 上传数据路径
    data_path: Optional[str] = None
    # 数据的报告
    data_repo: Optional[Union[str,dict,list]] = None
    # docker中的存储文件夹
    docker_save_dir : Optional[str] = None
    # docker中的数据路径
    docker_data_path: Optional[str] = None
    # debug历史记录
    debug_history: List[Dict[str, Any]] = []
    # 默认的向量库
    default_vectorstore: dict = {
        "Code": os.path.join(PROJECT_ROOT, "vectorstore/BRICK_code.faiss"),
        "Notebook": os.path.join(PROJECT_ROOT, "vectorstore/BRICK_notebook.faiss")
    }
    # LLM对step的执行判断
    execution: Optional[Union[str, dict]] = None
    # 函数执行的结果
    execution_result: Optional[dict] = None
    # 代码报错信息
    error_message: Optional[str] = None
    # 过去生成的错误代码
    error_code: Optional[str] = None
    # 最终答案
    final_answer: Optional[str] = None
    # 最终结果
    final_result: Optional[Union[str, dict, list]] = None
    # 关系查询结果表格 (DataFrame类型)
    relation_frame: Optional[Any] = None
    # 单次返回的未验证的function
    find_function: Optional[Union[str,list]] = None
    # 单次返回的tutorial
    find_tutorial: Optional[Union[str,list]] = None
    # 可用的function
    functions: list[Union[str,list]] = []
    # 所有生成的代码
    full_code: Optional[list[str]] = []
    # KG的schema
    kg_schema: dict = {"nodes": [], "edges": []}
    # 语言
    language: str = "English"
    # 是否生成了分析报告
    make_analysis: bool = False
    # 是否检查了数据
    make_env_check: bool = False
    # 是否生成了数据报告
    make_data_repo: bool = False
    # 是否生成了plan
    make_plan: bool = False
    # 是否进行了plan review
    make_review: bool = False
    # 记忆：每一轮agent的对话
    messages: Annotated[list[AnyMessage], add_messages]
    # 代码执行的notebook结果
    notebook_cells: list = []
    # notebook默认库
    notebooks_path: str = os.path.join(PROJECT_ROOT, "notebooks")
    # notebook的文本内容
    notebook_text: Optional[str] = None
    # 下一个agent
    next: Literal["env_checker", "supervisor", "translator", "data_analyzer", "analyze_planner", "analyse_planner", "planner", "planner_stepwise", "searcher", "BRICK_searcher", "bkg_searcher", "plan_reviewer", "plan_executor", "plan_checker", "plan_strategy", "coder", "code_runner","code_debugger","code_evaluator", "code_controller", "code_executer", "responder", "general_responder", "step_checker", "step_spliter", "parse_interact", "parse_update", "test", "verify", "END","query_agent"] = "supervisor"
    # Agent的系统输出
    output: Optional[Union[str, dict, list]] = None
    # 代码运行是否成功
    process_flag: int = -1
    # 预定义的plan
    predefined_plans: dict = {
        "trajectory_inference": [
            {"step": 1, "type": "preprocess", "details": "Use Scanpy to read the user's inputted h5ad data, visualize UMAP colored by cell annotation label (e.g. 'cell_type'), calculate PAGA graph based on cell type clustering, use 'connectivities_tree' with threshold 0.05."},
            {"step": 2, "type": "split", "details": "Extract connected components from the PAGA 'connectivities_tree', assign cells to different subgroups ('paga_cluster'), and separately subset Anndata objects for each subgroup."},
            {"step": 3, "type": "preprocess", "details": "For each subgroup, perform standard preprocessing steps again (neighbors, PCA, UMAP) to better represent the internal structure of each subgroup."},
            {"step": 4, "type": "retrieval", "details": "Use BRICK.qr.query_shortest_path function to retrieve shortest paths between pairs of nodes (cell types) within each subgroup, enriching the developmental relationships."},
            {"step": 5, "type": "filter", "details": "Use BRICK.pp.filter_results function to filter the retrieved paths and remove irrelevant or low-quality trajectories."},
            {"step": 6, "type": "integration", "details": "Use BRICK.pp.complete_results function to integrate the filtered retrieval results into the original PAGA-derived graph to reconstruct a more complete developmental trajectory."},
            {"step": 7, "type": "visualization", "details": "Visualize both the original and the enriched developmental graphs using BRICK.pl.static_visualize_network, and save the network plots in spring layout and tree layout formats."},
            {"step": 8, "type": "pseudotime_inference", "details": "Identify the root cell type (highest degree node) and perform pseudotime inference using Scanpy's DPT algorithm, then visualize pseudotime distribution on UMAP."},
            {"step": 9, "type": "interpretation", "details": "Use BRICK.inp.interpret_results function to interpret the cell developmental trajectory graph results."}
        ],
        'cell_annotation': [
            {"step": 1, "type": "preprocess", "details": "Use Scanpy to read user's inputted h5ad data and only do necessary preprocess steps for cell annotation."},
            {"step": 2, "type": "retrieval", "details": "Use BRICK.qr.query_neighbor function to find the neighbor of a node, and it can be used to annotate the cell type in the omics data."},
            {"step": 3, "type": "rank", "details": "Use BRICK.rk.rank_results function to rank the results after retrieval step for the cell type annotation."},
            {"step": 4, "type": "interpretation", "details": "Use BRICK.inp.interpret_results function to interpret the results after intergration step for the cell type annotation."}
        ],
        'undefined': []
    }
    # 展示的数据数量
    preview_n: int = 20
    # 输入的用户问题
    question: Optional[str] = None
    # 评审LLM生成的plan建议
    re_plan: Optional[Union[str, dict, list, int]] = None
    # 剩余的步骤
    remaining_steps: int = 200
    # 参考notebook路径
    reference_notebook_path: Optional[str] = None
    # 保存目录
    save_dir: str = "./"
    # 沙箱id
    sandbox_id: Optional[str] = None
    # 状态
    status: str = "NOT_FINISHED"
    # 单独的step
    step: list[Union[str, dict, list]] = []
    # 单独的step
    step_output: Optional[Union[str, dict]] = None
    # 单独的step的输出
    step_content:Optional[dict] = {} 
    # 计划步骤数
    step_num: int = 0
    # LLM的思考过程
    thought: Optional[str]= None
    # 翻译后的用户问题
    translated_question: Optional[str] = None
    # 更新的代码指令
    update_code: str = "Don't need to modify the current code."
    # 更新的数据/环境指令(env_checker)
    update_data_info: Optional[str] = None
    # 更新的数据报告(data_analyzer)
    update_data_repo: Optional[str] = None
    # 更新的方案(reviewer)
    update_instruction: list = []
    # 更新的方案
    update_plan: list = []
    # 用户的回答(analyze_planner)
    user_update_detail: Optional[str] = None
    # 是否是标准h5ad文件
    valid_data: bool = False

    class Config:
        arbitrary_types_allowed = True  # 允许 DataFrame 等任意类型
