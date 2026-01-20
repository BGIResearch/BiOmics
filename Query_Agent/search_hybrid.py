import os
import sys
from tkinter import N
from typing import List
import pandas as pd
import time
from typing import Annotated, Literal, Union, List, Dict, Any, Optional, Tuple
from neo4j.graph import Entity
# 使用相对路径，适配任何工作目录
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if project_root not in sys.path:
    sys.path.append(project_root)
brick_path = os.path.join(project_root, 'BRICK')
if brick_path not in sys.path:
    sys.path.insert(0, brick_path)


import BRICK

# Configure BRICK
url = "neo4j://10.224.28.66:7687"
auth = ("neo4j", "bmVvNGpwYXNzd29yZA==")

BRICK.config(url=url, auth=auth)
BRICK.config_llm(modeltype='ChatOpenAI', 
                 api_key="sk-kpsteSkpDGl1xBmDEcC7D51b968e43499092826f17286b55",  
                 base_url='http://10.224.28.80:3000/v1', 
                 llm_params={'model_name': 'qwen-max'})

client = BRICK.se.BRICKSearchClient()
print(f"字符串索引: {client.search_config.string_index_name}")
print(f"向量索引: {client.search_config.vector_index_name or '已禁用'}\n")

def search_entity(query_entities:List[Tuple[str, str]]):
    searched_res = []
    for query_entity in query_entities:
        time_stamp= str(int(time.time())) 
        print(f"Query Entity: {query_entity}")
        name, type = query_entity
        if type == "Gene" or type == "Protein":
            search_payload = {
                "query_id": time_stamp + "Gene|Protein",
                "options": {"top_k": 5, "return_diagnostics": True},
                "Gene|Protein": [name],
            }
            std, diag, resp = client.search_hybrid(
                search_payload,
                type_mix_override={"Gene|Protein": 1.0},
            )
            searched_res.append(diag["primary_name"].to_list())
        elif type == "Mutation":
            search_payload = {
                "query_id": time_stamp + "Mutation",
                "options": {"top_k": 5, "return_diagnostics": True},
                "Mutation": [name],
            }
            std, diag, resp = client.search_hybrid(
                search_payload,
                type_mix_override={"Mutation": 1.0},
            )
            searched_res.append(diag["primary_name"].to_list())
        elif type == "Chemical":
            search_payload = {
                "query_id": time_stamp + "Chemical",
                "options": {"top_k": 5, "return_diagnostics": True},
                "Chemical": [name],
            }
            std, diag, resp = client.search_hybrid(
                search_payload,
                type_mix_override={"Chemical": 1.0},
            )
            searched_res.append(diag["primary_name"].to_list())
        elif type == "Disease" or type == "Phenotype":
            search_payload = {
                "query_id": time_stamp + "Disease|Phenotype",
                "options": {"top_k": 5, "return_diagnostics": True},
                "Disease|Phenotype": [name],
            }
            std, diag, resp = client.search_hybrid(
                search_payload,
                type_mix_override={"Disease|Phenotype": 1.0},
            )
            searched_res.append(diag["primary_name"].to_list())
        elif type == "Process" or type == "Function" or type == "Pathway" or type == "Cell_Component":
            search_payload = {
                "query_id": time_stamp + "Process|Function|Pathway|Cell_Component",
                "options": {"top_k": 5, "return_diagnostics": True},
                "Process|Function|Pathway|Cell_Component": [name],
            }
            std, diag, resp = client.search_hybrid(
                search_payload,
                type_mix_override={"Process|Function|Pathway|Cell_Component": 1.0},
            )
            searched_res.append(diag["primary_name"].to_list())
        elif type == "Species":
            search_payload = {
                "query_id": time_stamp + "Species",
                "options": {"top_k": 5, "return_diagnostics": True},
                "Species": [name],
            }
            std, diag, resp = client.search_hybrid(
                search_payload,
                type_mix_override={"Species": 1.0},
            )
            searched_res.append(diag["primary_name"].to_list())
        elif type == "Cell" or type == "Tissue":
            search_payload = {
                "query_id": time_stamp + "Cell|Tissue",
                "options": {"top_k": 5, "return_diagnostics": True},
                "Cell|Tissue": [name],
            }
            std, diag, resp = client.search_hybrid(
                search_payload,
                type_mix_override={"Cell|Tissue": 1.0},
            )
            searched_res.append(diag["primary_name"].to_list())
        else:
            payload_hybrid_auto = {
                "query_id": time_stamp + "unknown",
                "options": {"top_k": 5, "return_diagnostics": True, "query_text": name},
            }
            std, diag, resp = client.search_hybrid(payload_hybrid_auto)
            searched_res.append(diag["primary_name"].to_list())
    return searched_res

if __name__ == "__main__":
    query_entities = [("crab", "Species"), ("Pigment cell", "Cell")]
    searched_res = search_entity(query_entities)
    print(searched_res)
