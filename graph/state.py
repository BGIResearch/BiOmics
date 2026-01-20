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
    notebook_text: Optional[str] = """
    
=== Cell 1 (markdown) ===

# Celltype Annotation

This tutorial shows how to use BRICK for cell type annotation based on the clustering results of scanpy. Through the differentially expressed genes of each cell cluster, we locate its corresponding cell type.

You can use this turtorial when the user ask : 'Annotate this dataset.'

## Load packages and data 



=== Cell 2 (code) ===

import BRICK
import scanpy as sc

url = "neo4j://10.224.28.66:7687"
auth = ("neo4j", "bmVvNGpwYXNzd29yZA==")  

BRICK.config(url=url, auth=auth)
BRICK.config_llm(modeltype='ChatOpenAI', 
                 api_key="sk-kpsteSkpDGl1xBmDEcC7D51b968e43499092826f17286b55",  
                 base_url='http://10.224.28.80:3000/v1', 
                 llm_params={'model_name': 'qwen-max'})


[stdout output]:
Graph database has been configured and initialized successfully.
LLM has been configured and initialized successfully.



=== Cell 3 (code) ===

adata = sc.read('../../../KG_annotation/adata_new1.h5ad')
adata.X = adata.layers['lognorm'].copy() # use log normalized dataset to calculated differential expression gene



=== Cell 4 (code) ===

adata


[output]:
AnnData object with n_obs × n_vars = 822 × 14821
    obs: 'sample_name', 'n_genes', 'n_counts', 'annotation', 'leiden', 'KGannotator2', 'delta_specific', 'gamma_specific', 'leiden_combined', 'KGannotator_refinement'
    var: 'Ensembl_id', 'Symbol', 'NCBI_id', 'MGI_id', 'mean', 'std'
    uns: 'KGannotator2_colors', 'KGannotator_refinement_colors', 'leiden', 'log1p', 'neighbors', 'pca', 'rank_genes_groups', 'umap'
    obsm: 'X_pca', 'X_umap', 'annotation_au'
    varm: 'PCs'
    layers: 'lognorm'
    obsp: 'connectivities', 'distances'


=== Cell 5 (markdown) ===

## Cell cluster



=== Cell 6 (code) ===

sc.tl.leiden(adata)



=== Cell 7 (markdown) ===

## Differential expressed gene



=== Cell 8 (code) ===

sc.tl.rank_genes_groups(adata, groupby='leiden', pts=True)
BRICK.pp.rank_genes_groups2df(adata)



=== Cell 9 (markdown) ===

BRICK has implement a function to make the rank_genes_groups result as a dictionary where key is the cell cluster and value is a dataframe to record DEG for each cluster



=== Cell 10 (markdown) ===

## Annotation according to the DEG

Take cell cluster 0 as example, we will show how the BRICK works to transform a gene list into a cell type



=== Cell 11 (code) ===

genelist = BRICK.pp.filter_marker(adata.uns['rank_genes_groups_df']['0'], topgenenumber=10)
genelist


[output]:
['Meg3',
 'Prss53',
 'Ero1b',
 'Ptprn',
 'Neat1',
 'Chga',
 'Prlr',
 'Syt7',
 'Scg3',
 'Igf1r']


=== Cell 12 (markdown) ===

### query graph to get possible celltypes



=== Cell 13 (code) ===

query_df = BRICK.qr.query_neighbor(genelist, source_entity_type='Gene', relation='marker_of', target_entity_type='Cell')
query_df.head(3)


[output]:
                                                                                         path.0.def  \
0           secretogranin III<loc>:9 D|9 42.32 cM<xref>ENSEMBL:ENSMUSG00000032181|MGI:103032</xref>   
1  maternally expressed 3<loc>:12 F1|12 60.25 cM<xref>ENSEMBL:ENSMUSG00000021268|MGI:1202886</xref>   
2         synaptotagmin VII<loc>:19 A|19 6.58 cM<xref>MGI:1859545|ENSEMBL:ENSMUSG00000024743</xref>   

    path.0.id path.0.name  \
0  NCBI:20255        Scg3   
1  NCBI:17263        Meg3   
2  NCBI:54525        Syt7   

                                               path.0.synonym path.0.type  \
0                                           1B1075|Chgd|SgIII        Gene   
1  2900016C05Rik|3110050O07Rik|6330408G06Rik|D12Bwg1266e|Gtl2        Gene   
2                                               B230112P13Rik        Gene   

      path.1         path.1.condition  \
0  marker_of  [Undef, UBERON:0000955]   
1  marker_of         [UBERON:0002107]   
2  marker_of         [UBERON:0003126]   

                                                                                                                     path.1.info_source  \
0                               [Cell_Taxonomy:Human Cell Landscape:PMID:32214235, Cell_Taxonomy:CellMatch:PMID:29775597, SCT000000946]   
1                                                                                                                        [SCT000000976]   
2  [Cell_Taxonomy:CellMarker:PMID:30069046, Cell_Taxonomy:CellMarker:PMID:30069044, CellMarker:PMID:30069046, CellMarker:PMID:30069044]   

   path.1.info_source_length                      path.1.original_relation  \
0                          3                [marker_of, marker_of, DEG_of]   
1                          1                                      [DEG_of]   
2                          4  [marker_of, marker_of, marker_of, marker_of]   

  path.1.relation path.1.relation_confidence  \
0       marker_of          [1.0, 1.0, 0.617]   
1       marker_of                    [0.836]   
2       marker_of               [1, 1, 1, 1]   

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     path.2.def  \
0  A class of large neuroglial (macroglial) cells in the central nervous system - the largest and most numerous neuroglial cells in the brain and spinal cord. Astrocytes (from 'star' cells) are irregularly shaped with many long processes, including those with 'end feet' which form the glial (limiting) membrane and directly and indirectly contribute to the blood-brain barrier. They regulate the extracellular ionic and chemical environment, and 'reactive astrocytes' (along with microglia) respond to injury.<xref>CALOHA:TS-0060|BTO:0000099|FMA:54537</xref>   
1                                                                                                                                                                                                                                                                                          A connective tissue cell which secretes an extracellular matrix rich in collagen and other macromolecules. Flattened and irregular in outline with branching processes; appear fusiform or spindle-shaped.<xref>CALOHA:TS-0362|VHOG:0001482|BTO:0000452|NCIT:C12482|FMA:63877</xref>   
2                                                                                                                                                                                                                                                                                                                                                                                                                                           A neuron that is capable of some hormone secretion in response to neuronal signals.<xref>BTO:0002691|FMA:83810|FBBT:00005130</xref>   

    path.2.id          path.2.name                             path.2.synonym  \
0  CL:0000127            astrocyte                            astrocytic glia   
1  CL:0000057           fibroblast                                      Undef   
2  CL:0000165  neuroendocrine cell  neurosecretory cell|neurosecretory neuron   

  path.2.type  
0        Cell  
1        Cell  
2        Cell  


=== Cell 14 (markdown) ===

### Prune parent celltype
This could be optional, we can delete parent node to make the annotation more specific. Here we only deleted the those parent nodes with only one child node



=== Cell 15 (code) ===

all_possible_cell = list(query_df['path.2.name'].unique())
parent2children = BRICK.qr.query_relation(all_possible_cell, target_entity_type='Cell', relation='is_a', directed=True)
parent_cells = [x for x, y in parent2children.groupby('path.2.name') if y.shape[0] == 1]
only_child_parent = parent2children.loc[parent2children['path.2.name'].isin(parent_cells)]
only_child_parent = dict(zip(only_child_parent['path.2.name'], only_child_parent['path.0.name']))
query_df['path.2.name'] = [ only_child_parent[x] if x in only_child_parent else x for x in query_df['path.2.name'] ]



=== Cell 16 (markdown) ===

### rank target entity to get the most possible celltype



=== Cell 17 (code) ===

target_df = BRICK.rk.rank_voting(query_df, metrics=['match_count', 'match_probability', 'info_source_count'])



=== Cell 18 (code) ===

target_df.sort_values('path.2.rank_voting', ascending = True).head(3)


[output]:
                                                                                                                                                                                                          path.0.name  \
142                                                                                     [Ptprn, Prlr, Chga, Scg3, Prss53, Syt7, Ero1b, Meg3, Neat1, Igf1r, Chga, Igf1r, Meg3, Prlr, Ptprn, Scg3, Prss53, Syt7, Ero1b]   
147  [Chga, Prlr, Ptprn, Scg3, Prss53, Syt7, Ero1b, Meg3, Meg3, Chga, Prlr, Ptprn, Scg3, Prss53, Syt7, Ero1b, Meg3, Chga, Prlr, Ptprn, Scg3, Prss53, Syt7, Ero1b, Meg3, Chga, Prlr, Ptprn, Scg3, Prss53, Syt7, Ero1b]   
44                                                                                                                                                               [Prss53, Syt7, Scg3, Prlr, Chga, Meg3, Ero1b, Ptprn]   

                                                                                                                                                                                                                                                                                                                                                      path.1.relation  \
142                                                                                                                                                 [marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of]   
147  [marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of]   
44                                                                                                                                                                                                                                                                           [marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of, marker_of]   

                                                                            path.1.info_source_length  \
142                                         [2, 3, 2, 2, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2, 2, 1, 3]   
147  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]   
44                                                                           [1, 1, 5, 1, 8, 1, 1, 3]   

                                                                                                                                                                                                                                                           path.1.relation_confidence  \
142  [[1.0, 0.5449999999999999], [1.0, 0.6146, 0.6937], [1.0, 0.598], [1.0, 0.611], [1.0, 0.7399, 0.7978], [1.0, 0.43899999999999995], [1], [0.5070000000000001], [0.44399999999999995], [0.524], [1, 1], [1], [1], [1], [1.0, 0.9508], [1.0, 0.9156], [1.0, 0.9386], [1], [1, 1, 1]]   
147                                                                                                                  [[1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1]]   
44                                                                                                                                                                  [[1], [1], [1, 1, 1, 1, 1], [1], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.7430000000000001, 0.4481], [1], [1], [1, 1, 1]]   

      path.2.id                  path.2.name path.2.type  path.2.match_count  \
142  CL:0000168       type B pancreatic cell        Cell                  10   
147  CL:0000506  type G enteroendocrine cell        Cell                   8   
44   CL:0000164         enteroendocrine cell        Cell                   8   

     path.2.match_probability  path.2.info_source_count  path.2.rank_voting  
142                  0.044759                        33            1.000000  
147                  0.011426                        32           13.833333  
44                   0.011426                        21           14.166667  


=== Cell 19 (code) ===

query_target_df = target_df.sort_values('path.2.rank_voting', ascending = True).head(5)
query_target_df = query_target_df[['path.0.name', 'path.2.name', 'path.2.match_count', 'path.2.match_probability', 'path.2.info_source_count', 'path.2.rank_voting']]




=== Cell 20 (code) ===

print(BRICK.inp.cell_annotation_analysis(query_target_df))


[stdout output]:
Based on the queried table, the most reasonable annotation for this cell cluster is **type B pancreatic cell**. Here's why:

1. **Rank Voting and Match Count**: The "type B pancreatic cell" has the highest rank (rank_voting = 1.000000) and a relatively high match count (10), indicating a strong match with the differential gene set.

2. **Subset-Superset Relations**: The other candidate cell types (e.g., type G enteroendocrine cell, enteroendocrine cell, type EC2 enteroendocrine cell, type EC1 enteroendocrine cell) are more general or specific subtypes of enteroendocrine cells. Since "type B pancreatic cell" is a distinct and more specific cell type, it is the preferred annotation.

3. **Mixture of Sub-Cell Types**: There is no strong evidence from the data to suggest that this cluster is a mixture of two sub-cell types. The high rank and match count for "type B pancreatic cell" indicate a single, well-defined cell type.

In summary, the most reasonable annotation for this cell cluster is **type B pancreatic cell**. No further sub-cell type refinement is necessary based on the current data.



=== Cell 21 (markdown) ===

According to the result, the most possible annotation for cluster 0 is type B pancreatic cell



=== Cell 22 (markdown) ===

## run cell annotation for all cell clusters 



=== Cell 23 (code) ===

import tqdm
cellcluster2celltype = {}
cellcluster2annoation_result = {}

for x, y in tqdm.tqdm(adata.uns['rank_genes_groups_df'].items()):
    genelist = BRICK.pp.filter_marker(y, topgenenumber=10)
    query_df = BRICK.qr.query_neighbor(genelist, source_entity_type='Gene', relation='marker_of', target_entity_type='Cell')
    all_possible_cell = list(query_df['path.2.name'].unique())
    parent2children = BRICK.qr.query_relation(all_possible_cell, target_entity_type='Cell', relation='is_a', directed=True)
    parent_cells = [x for x, y in parent2children.groupby('path.2.name') if y.shape[0] == 1]
    only_child_parent = parent2children.loc[parent2children['path.2.name'].isin(parent_cells)]
    only_child_parent = dict(zip(only_child_parent['path.2.name'], only_child_parent['path.0.name']))
    query_df['path.2.name'] = [ only_child_parent[x] if x in only_child_parent else x for x in query_df['path.2.name'] ]
    target_df = BRICK.rk.rank_voting(query_df, metrics=['match_count', 'match_probability', 'info_source_count'])
    cellcluster2annoation_result[x] = target_df.sort_values('path.2.rank_voting', ascending = True).head(5)
    cellcluster2celltype[x] = list(cellcluster2annoation_result[x]['path.2.name'])[0]


[stderr output]:
100%|██████████████████████████████████████████████████████████████████████████████████| 11/11 [00:08<00:00,  1.28it/s]



=== Cell 24 (code) ===

cellcluster2celltype


[output]:
{'0': 'type B pancreatic cell',
 '1': 'type B pancreatic cell',
 '2': 'type D enteroendocrine cell',
 '3': 'pancreatic ductal cell',
 '4': 'epithelial cell',
 '5': 'type B pancreatic cell',
 '6': 'endothelial cell',
 '7': 'pancreatic ductal cell',
 '8': 'pancreatic stellate cell',
 '9': 'macrophage',
 '10': 'neuroendocrine cell'}


=== Cell 25 (code) ===

adata.obs['celltype'] = [cellcluster2celltype[x] for x in adata.obs['leiden']]
sc.pl.umap(adata, color = 'celltype')
"""
    # 下一个agent
    next: Literal["env_checker", "supervisor", "translator", "data_analyzer", "analyze_planner", "analyse_planner", "planner", "planner_stepwise", "searcher", "BRICK_searcher", "bkg_searcher", "plan_reviewer", "plan_executor", "plan_checker", "plan_strategy", "coder", "code_runner","code_debugger","code_evaluator", "code_controller", "code_executer", "responder", "general_responder", "step_checker", "step_spliter", "parse_interact", "parse_update", "test", "verify", "END","query_agent"] = "supervisor"
    # Agent的系统输出
    output: Optional[Union[str, dict, list]] = None
    #代码运行是否成功
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
    user_update_detail: list = []
    # 是否是标准h5ad文件
    valid_data: bool = False

    class Config:
        arbitrary_types_allowed = True  # 允许 DataFrame 等任意类型
