import scanpy as sc
import numpy as np
import scipy.sparse
from langchain_core.tools import tool
from pydantic import BaseModel

class Input(BaseModel):
    data_path: str = ""

class OutputData(BaseModel):
    description: str = ""
    decision: bool = False

class Output(BaseModel):
    output: OutputData = OutputData()

@tool(args_schema=Input)
def validate_h5ad_structure_tool(data_path) -> Output:
    """
    Use this tool to verify whether the h5ad file contains a valid expression matrix (ndarray/csr/csc). 
    
    Args:
        data_path: The path to the h5ad file.

    Returns:
        output: A dictionary containing the verification result description and the boolean decision.
    """
    try:
        print(data_path)
        adata = sc.read_h5ad(data_path)
        # 基础数据结构
        assert isinstance(adata.X, (np.ndarray, 
                                  scipy.sparse.csr_matrix,
                                  scipy.sparse.csc_matrix)), \
               "X must be an ndarray/csr/csc matrix." 
        output = {"description": "The h5ad file contains a valid expression matrix.", "decision": True}
        return output
    except Exception as e:
        output = {"description": f"The h5ad file does not contain a valid expression matrix: {str(e)}", "decision": False}
        return output
    finally:
        if 'adata' in locals():
            del adata  # 释放内存

if __name__ == "__main__":
    result = validate_h5ad_structure_tool.invoke({"data_path": "./enhanced_simulated.h5ad"})
    print(result)
    