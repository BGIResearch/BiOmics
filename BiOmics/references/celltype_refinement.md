# Sub-celltype Refinement with brick CLI

## Task Overview

Cell types are hierarchical concepts. For example, CD4+ T cell is a type of T cell. During single-cell annotation, it is often necessary to determine cell annotation results at different levels. Sub-celltype refinement is an important analysis task—after routine annotation, further subdivision of cell subpopulations is often required.

This tutorial demonstrates how to use brick CLI tools to refine cell subtypes, using pancreatic endocrine cell subdivision into pancreatic PP cell and pancreatic D cell as an example.

## Prerequisites

- BiOmics-BRICK package installed: `pip install BiOmics-BRICK`
- AnnData format single-cell data prepared (containing cell annotation results)

## Analysis Workflow

### Step 1: Query Markers for Candidate Cell Types

First, query marker genes for candidate cell types from the knowledge graph. In this example, we need to query markers for pancreatic D cell and pancreatic PP cell. If you are unsure which candidate cell types to query, you can first use the differential genes from this cell cluster and query candidate cell types once via the produce relation.

```bash
brick qr "pancreatic D cell, pancreatic PP cell" \
    -r produce \
    --st Cell \
    --tt Gene \
    -n 500 \
    -o ./sub_adata_markers.csv
```

**Parameter Description:**
- `"pancreatic D cell, pancreatic PP cell"`: Candidate cell types, comma-separated
- `-r produce`: Query marker evidence edges in BiOmics-KG. The concept alias `marker_of` is accepted and resolved to `produce`.
- `--st Cell`: Source entity type is Cell
- `--tt Gene`: Target entity type is Gene
- `-n 500`: Return up to 500 results
- `-o ./sub_adata_markers.csv`: Output results to CSV file

### Step 2 (Optional): Rank Query Results

Use brick rk to rank query results and filter more reliable markers.

```bash
brick rk -i ./sub_adata_markers.csv \
    --metrics info_source_count \
    -o ./sub_adata_markers_ranked.csv
```

**Parameter Description:**
- `-i ./sub_adata_markers.csv`: Input brick qr output file
- `--metrics info_source_count`: Use information source count as ranking metric
- `-o ./sub_adata_markers_ranked.csv`: Output ranked results

### Step 3: Prepare AnnData Subset

Extract the cell subpopulation to be refined in Python. Assume the cells to be refined are labeled as "somatostatin secreting cell" in the `KGannotator2` column:

```python
import scanpy as sc

# Read raw data
adata = sc.read("../KG_annotation/adata_new1.h5ad")

# Extract cell subpopulation to be refined
adata_sub = adata[adata.obs['KGannotator2'] == 'somatostatin secreting cell', :].copy()

# Ensure correct expression matrix (raw counts or log normalized)
adata_sub.X = adata_sub.layers['lognorm']

# Save subset for subsequent analysis
adata_sub.write("./adata_sub.h5ad")
```

### Step 4: Build Integrated Graph and Compute an Embedding

Use brick emb to combine single-cell data with knowledge graph query results. The default path builds an adjacency/PCA/UMAP representation; GRACE can instead learn a graph-contrastive representation when graph topology should contribute to subtype refinement.

```bash
brick emb -a ./adata_sub.h5ad \
    -i ./sub_adata_markers.csv \
    -o ./adata_graph.h5ad
```

Optional GRACE analysis:

```bash
pip install -e ".[grace]"
brick emb -a ./adata_sub.h5ad \
    -i ./sub_adata_markers.csv \
    --method grace \
    --cell-neighbor \
    --grace-epochs 300 \
    --seed 20241204 \
    -o ./adata_grace.h5ad
```

**Parameter Description:**
- `-a ./adata_sub.h5ad`: Input AnnData file (subset generated in Step 3)
- `-i ./sub_adata_markers.csv`: Input brick qr query results
- `-o ./adata_graph.h5ad`: Output integrated graph file

**Output file contains:**
- `X_pca`: PCA dimensionality reduction results
- `X_umap`: UMAP embedding coordinates
- `uns['KG_Omics_distance']`: Distance matrix between knowledge graph entities and omics data

For `--method grace`, the output instead stores learned representations in `.X` and `.obsm['X_grace']`, training metadata/losses in `.uns['grace']`, integrated adjacency in `.obsp['integrated_adjacency']`, UMAP in `.obsm['X_umap']`, and cosine distance by default in `.uns['KG_Omics_distance']`. See `grace_embedding.md` before interpreting candidate proximity.

### Step 5: Subtype Refinement Based on Distance

Load results in Python and perform subdivision based on cell distance to candidate cell types:

```python
import scanpy as sc

# Read integrated graph and original subset
adata_graph = sc.read("./adata_graph.h5ad")
adata_sub = sc.read("./adata_sub.h5ad")

# Assign cell types based on distance matrix
# The column name with minimum value per row is the final annotation result
adata_sub.obs["refined_celltype"] = adata_graph.uns["KG_Omics_distance"].idxmin(axis=1)

# Visualize refinement results
sc.pl.umap(adata_sub, color=["refined_celltype", "KGannotator2"])
```

### Step 6: Write Refined Annotations Back to the Original h5ad

Map the refined cell types from the subset back to the full AnnData object and save an updated file. Cells that were not part of the refinement keep their original annotations.

```python
import scanpy as sc

# Read the original full dataset and the refined subset from Step 5
adata = sc.read("../KG_annotation/adata_new1.h5ad")
adata_sub = sc.read("./adata_sub.h5ad")

# Map refined labels back to the full object by cell barcode
refined = adata_sub.obs["refined_celltype"]
adata.obs["refined_celltype"] = adata.obs_names.map(refined).astype("category")

# Write refined labels into the next annotation hierarchy level
# Non-refined cells keep their parent label from KGannotator2
parent_col = "KGannotator2"
new_col = "KGannotator3"
adata.obs[new_col] = adata.obs[parent_col].astype(str)
adata.obs.loc[refined.index, new_col] = refined.values
adata.obs[new_col] = adata.obs[new_col].astype("category")

# Save the updated AnnData file
adata.write("./adata_refined.h5ad")
```

**Parameter Description:**
- `../KG_annotation/adata_new1.h5ad`: Original full AnnData file (same as Step 3)
- `refined_celltype`: Refinement result for the subset; mapped to all cells in the full object (NaN for non-refined cells)
- `KGannotator3`: Hierarchical annotation column; refined cells get subtypes, others keep `KGannotator2` labels
- `./adata_refined.h5ad`: Output file with refinement results written back

## Complete Workflow

```bash
# 1. Query markers for candidate cell types
brick qr "pancreatic D cell, pancreatic PP cell" \
    -r produce \
    --st Cell \
    --tt Gene \
    -n 500 \
    -o ./sub_adata_markers.csv

# 2. Rank query results (optional)
brick rk -i ./sub_adata_markers.csv \
    --metrics info_source_count \
    -o ./sub_adata_markers_ranked.csv

# 3. Build integrated graph (using original or ranked results)
brick emb -a ./adata_sub.h5ad \
    -i ./sub_adata_markers.csv \
    -o ./adata_graph.h5ad
```

```python
# 4. Perform subtype refinement in Python
import scanpy as sc

adata_graph = sc.read("./adata_graph.h5ad")
adata_sub = sc.read("./adata_sub.h5ad")

# Assign cell types based on nearest distance
adata_sub.obs["refined_celltype"] = adata_graph.uns["KG_Omics_distance"].idxmin(axis=1)

# Visualize results
sc.pl.umap(adata_sub, color=["refined_celltype", "KGannotator2"])
```

```python
# 5. Write refined annotations back to the original h5ad
import scanpy as sc

adata = sc.read("../KG_annotation/adata_new1.h5ad")
adata_sub = sc.read("./adata_sub.h5ad")

refined = adata_sub.obs["refined_celltype"]
adata.obs["refined_celltype"] = adata.obs_names.map(refined).astype("category")

parent_col = "KGannotator2"
new_col = "KGannotator3"
adata.obs[new_col] = adata.obs[parent_col].astype(str)
adata.obs.loc[refined.index, new_col] = refined.values
adata.obs[new_col] = adata.obs[new_col].astype("category")

adata.write("./adata_refined.h5ad")
```

## Result Interpretation

- **refined_celltype**: Refined cell type (pancreatic D cell or pancreatic PP cell) for cells in the refined subset; NaN for all other cells in the full object
- **KGannotator3**: Updated hierarchical annotation column in the full h5ad; refined cells receive subtype labels, others retain `KGannotator2` values
- **adata_refined.h5ad**: Full AnnData file with refinement results written back
- **KG_Omics_distance**: Distance matrix, can be used to evaluate refinement confidence
- Compare UMAP plots before and after refinement to assess subdivision quality

## Notes

1. **Data Preparation**: Ensure AnnData `.X` or `layers` contains the correct expression matrix
2. **Candidate Cell Type Selection**: Choose reasonable candidate cell types based on prior knowledge
3. **Query Parameter Adjustment**: Adjust `-n` parameter according to data scale to obtain sufficient markers
4. **Distance Matrix Interpretation**: Smaller values in `KG_Omics_distance` matrix indicate higher similarity
5. **Write-Back Strategy**: Step 6 only updates cells present in `adata_sub`; adjust `parent_col` / `new_col` if your annotation columns use different names
