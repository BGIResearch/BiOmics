import os
import sys
from typing import List
import pandas as pd

from neo4j.graph import Entity
# 使用相对路径，适配任何工作目录
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
if project_root not in sys.path:
    sys.path.append(project_root)

# 导入配置加载器
sys.path.insert(0, os.path.join(project_root, 'Query_Agent'))
from config_loader import configure_brick


def most_related_query(question, entity_name, entity_type, file_path, target_type) -> tuple:
    
    import BRICK
    import scanpy as sc

    # Configure BRICK - 从配置文件加载
    configure_brick(model_key="BASIC_MODEL")
    
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
