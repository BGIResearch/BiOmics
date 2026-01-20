import os
import sys
from typing import List
import pandas as pd

from neo4j.graph import Entity
# 使用相对路径，适配任何工作目录
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
if project_root not in sys.path:
    sys.path.append(project_root)


def relation_query(question: str, entity_name: List[str], entity_type: str, file_path: str) -> tuple:
    import BRICK
    import scanpy as sc

    # Configure BRICK
    url = "neo4j://10.224.28.66:7687"
    auth = ("neo4j", "bmVvNGpwYXNzd29yZA==")

    BRICK.config(url=url, auth=auth)
    BRICK.config_llm(modeltype='ChatOpenAI', 
                    api_key="sk-kpsteSkpDGl1xBmDEcC7D51b968e43499092826f17286b55",  
                    base_url='http://10.224.28.80:3000/v1', 
                    llm_params={'model_name': 'qwen-max'})
    relation_frame = BRICK.qr.query_relation(source_entity_set=entity_name, source_entity_type=entity_type, return_type="dataframe")  # type: ignore
    relation_frame.to_csv(file_path)  # type: ignore
    ans = BRICK.inp.interpret_query(question, relation_frame)
    return ans, relation_frame

if __name__ == "__main__":
    question = "What is the relationship between Isl1 and pp cell?"
    entity_name = ["Isl1","PP cell"]
    entity_type = "Gene | Cell"
    print(relation_query(question, entity_name, entity_type, "relation.csv"))
