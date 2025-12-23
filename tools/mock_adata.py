import anndata as ad
import numpy as np
import scanpy as sc
from langchain_core.tools import tool
from pydantic import BaseModel

class Input(BaseModel):
    n_cells: int = 200
    n_genes: int = 1000
    output_path: str = "./simulated.h5ad"
    cell_type: list = ["T cell", "B cell", "Monocyte"]

class Output(BaseModel):
    output: dict = {}

@tool(args_schema=Input)
def generate_enhanced_simulated_data_tool(n_cells=200, n_genes=1000, output_path="./enhanced_simulated.h5ad", cell_type=["T cell", "B cell", "Monocyte"]) -> Output:
    """  
    Use this tool to generate augmented simulated single-cell data and save as an h5ad file.
    
    Parameters:
        n_cells: Number of cells  
        n_genes: Number of genes  
        output_path: Output file path  
        cell_type: List of cell types  
        
    Returns:
        output: A dictionary containing the path to the generated h5ad file in the format {"data_path": "path/to/file.h5ad"}  
    """
    X = np.random.poisson(lam=1.0, size=(n_cells, n_genes))
    adata = ad.AnnData(X)

    adata.var_names = [f"Gene{i}" for i in range(n_genes)]
    adata.obs_names = [f"Cell{i}" for i in range(n_cells)]

    cell_type = cell_type
    batches = ["Batch1", "Batch2"]
    adata.obs["cell_type"] = np.random.choice(cell_type, size=n_cells)
    adata.obs["batch"] = np.random.choice(batches, size=n_cells)
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1)

    adata.var["gene_id"] = [f"ENSG{i:06d}" for i in range(n_genes)]
    adata.var["highly_variable"] = np.random.choice([True, False], size=n_genes)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata.raw = adata

    sc.pp.pca(adata, n_comps=50)
    adata.obsm["X_pca"] = adata.obsm["X_pca"]

    adata.write_h5ad(output_path)
    output = {"data_path": output_path}
    return output

if __name__ == "__main__":
    path = generate_enhanced_simulated_data_tool.invoke({"output_path":"./test.h5ad"})
    print(path)
    if isinstance(path, dict):
        adata = sc.read_h5ad(path["data_path"])
        print(adata)