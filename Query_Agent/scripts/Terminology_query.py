import os
import sys

from neo4j.graph import Entity
# 使用相对路径，适配任何工作目录
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
if project_root not in sys.path:
    sys.path.append(project_root)


def terminology_query(question, entity_name, file_path, query_type) -> tuple:    
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

    
    neighbor_df = BRICK.qr.query_neighbor(
        source_entity_set=entity_name,
        relation=None,
        source_entity_type=query_type,
        target_entity_type=None,
        multi_hop=1,
        directed=True,
        query_attribution="name",
        return_type="dataframe",
    )

    target_df = BRICK.rk.match_count(neighbor_df)
    grouped = target_df.groupby('path.2.type').head(10)
    grouped.to_csv(file_path, index=False)
    print("type",type(grouped))
    ans = BRICK.inp.interpret_query(question, grouped)
    return ans, grouped

if __name__ == "__main__":
    question = "What is the gene BRCA1?"
    entity_name = "BRCA1"
    print(terminology_query(question, entity_name, "output.csv", "Gene"))
