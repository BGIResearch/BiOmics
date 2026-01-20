import os
import sys
from typing import List
import pandas as pd

from neo4j.graph import Entity
# 使用相对路径，适配任何工作目录
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
if project_root not in sys.path:
    sys.path.append(project_root)



def most_related_query(question, entity_name, entity_type, file_path, target_type) -> tuple:
    
    import BRICK
    import scanpy as sc

    # Configure BRICK
    url = "neo4j://10.224.28.66:7687"
    auth = ("neo4j", "bmVvNGpwYXNzd29yZA==")

    BRICK.config(url=url, auth=auth)
    BRICK.config_llm(modeltype='ChatOpenAI', 
                    api_key="sk-kpsteSkpDGl1xBmDEcC7D51b968e43499092826f17286b55",  
                    base_url='http://10.224.28.80:3000/v1', 
                    llm_params={'model_name': 'qwen3-max'})
    rel_frame = BRICK.qr.query_neighbor(
        source_entity_set=entity_name, 
        source_entity_type=entity_type,
        target_entity_type=target_type,
        return_type="dataframe")
    most_rel_frame = BRICK.rk.info_source_count(rel_frame)
    most_rel_frame.to_csv(file_path)
    ans = BRICK.inp.interpret_query(question, most_rel_frame)
    return ans, most_rel_frame

if __name__ == "__main__":
    question = "What are the most related diseases to Isl1?"
    entity_name = "Isl1"
    entity_type = "Gene"
    target_type = "Disease"
    print(most_related_query(question, entity_name, entity_type, "output.csv", target_type))
