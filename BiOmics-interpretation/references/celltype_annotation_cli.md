# Cell Type Annotation with brick CLI

## Task Overview

Cell type annotation is a core step in single-cell RNA sequencing data analysis. This tutorial demonstrates how to use brick CLI tools to complete an end-to-end cell type annotation workflow, including data preprocessing, differential expression analysis, knowledge graph querying, and cell type inference.

**Note: This workflow is for reference only. Please adjust according to user needs and real-time conditions.**

## Prerequisites

- BiOmics-BRICK package installed: `pip install BiOmics-BRICK`
- scanpy installed: `pip install scanpy`
- AnnData format single-cell data prepared (containing raw counts and normalized data)

## Analysis Workflow

### Step 1: Data Preprocessing and Clustering

Use scanpy to complete data preprocessing, dimensionality reduction, clustering, and differential expression analysis.

> **Note**: The following code is for demonstration only. Please adjust preprocessing parameters (neighbor count, PC count, resolution, etc.) according to actual data characteristics.

```python
import scanpy as sc
import pandas as pd

adata = sc.read("./your_data.h5ad")
adata.X = adata.layers['lognorm'].copy()

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
sc.tl.rank_genes_groups(adata, groupby='leiden', pts=True, use_raw=False)

adata.write("./adata_clustered.h5ad")
```

### Step 2: Extract Marker Genes for Each Cluster

Extract top marker genes for each cluster and save as text files for subsequent brick queries.

> **Note**: The following code demonstrates basic marker extraction logic. Please adjust filtering conditions (pval, logFC thresholds) and extraction quantity according to data conditions.

```python
clusters = adata.obs['leiden'].cat.categories

for cluster in clusters:
    df = pd.DataFrame({
        'gene': adata.uns['rank_genes_groups']['names'][cluster],
        'pval_adj': adata.uns['rank_genes_groups']['pvals_adj'][cluster],
        'logfc': adata.uns['rank_genes_groups']['logfoldchanges'][cluster]
    })
    # Filter significant differentially expressed genes (adjust thresholds based on data)
    filtered = df[(df['pval_adj'] < 0.05) & (df['logfc'] > 0)]
    top_genes = filtered.head(10)['gene'].tolist()

    with open(f"./cluster_{cluster}_markers.txt", 'w') as f:
        f.write('\n'.join(top_genes))
```

### Step 3: Query Cell Types Using brick CLI

Use brick qr to query possible cell types for each cluster's marker genes.

```bash
# Query cell types for each cluster (loop processing)
for cluster in 0 1 2 3 4 5 6 7 8; do
    echo "Querying cluster $cluster..."
    brick qr -i ./cluster_${cluster}_markers.txt \
        --st Gene \
        --tt Cell \
        -r produce \
        -n 500 \
        -o ./cluster_${cluster}_query.csv
done
```

**Parameter Description:**
- `-i ./cluster_${cluster}_markers.txt`: Input marker genes file
- `--st Gene`: Source entity type is Gene
- `--tt Cell`: Target entity type is Cell
- `-r produce`: Query marker evidence edges in BiOmics-KG. The concept alias `marker_of` is also accepted by `brick qr`, but it is resolved to KG relation `produce`.
- `-n 500`: Return up to 500 results
- `-o ./cluster_${cluster}_query.csv`: Output query results

### Step 4: Rank Cell Type Candidates Using brick rk

Rank query results to obtain the most likely cell types.

```bash
# Rank query results for each cluster
for cluster in 0 1 2 3 4 5 6 7 8; do
    echo "Ranking cluster $cluster..."
    brick rk -i ./cluster_${cluster}_query.csv \
        --rank \
        --metrics "match_count,info_source_count" \
        --prune-parent \
        -o ./cluster_${cluster}_ranked.csv
done
```

**Parameter Description:**
- `-i ./cluster_${cluster}_query.csv`: Input brick qr output file
- `--rank`: Perform ranking analysis
- `--metrics "match_count,info_source_count"`: Use match count and information source count as ranking metrics
- `--prune-parent`: Prune parent nodes, keeping the most specific cell types
- `-o ./cluster_${cluster}_ranked.csv`: Output ranking results

### Step 5: Integrate Results and Assign Cell Types

Read brick rk ranking results and assign the most likely cell type to each cluster.

> **Note**: The following code demonstrates basic result integration. Please adjust assignment strategy (e.g., take top 1 or combine multiple metrics) according to actual conditions.

```python
import pandas as pd

cluster_annotations = {}

for cluster in adata.obs['leiden'].unique():
    cluster_str = str(cluster)
    try:
        df = pd.read_csv(f"./cluster_{cluster_str}_ranked.csv")
        if len(df) > 0:
            top_hit = df.sort_values('path.2.rank_voting').iloc[0]
            cluster_annotations[cluster_str] = top_hit['path.2.name']
        else:
            cluster_annotations[cluster_str] = "Unknown"
    except:
        cluster_annotations[cluster_str] = "Unknown"

adata.obs['celltype'] = adata.obs['leiden'].map(cluster_annotations)
adata.write("./adata_annotated.h5ad")
```

### Step 6: Visualize Annotation Results

Generate UMAP plots comparing clustering results with brick annotation results.

> **Note**: The following code is for visualization demonstration. Please adjust plotting parameters according to data scale and presentation needs.

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='Leiden Clusters')
sc.pl.umap(adata, color='celltype', ax=axes[1], show=False, title='brick Annotation')
plt.tight_layout()
plt.savefig("./celltype_annotation.png", dpi=150, bbox_inches='tight')
```


## Result Interpretation

- **cluster_X_query.csv**: brick qr raw query results, containing all candidate cell types
- **cluster_X_ranked.csv**: brick rk ranking results, containing ranking scores and statistics
- **adata_annotated.h5ad**: Complete data with cell type annotations
- **annotation_umap.png**: Visualization result plot

## Notes

1. **Data Quality**: Ensure input data has undergone appropriate quality control and normalization
2. **Clustering Resolution**: Leiden algorithm's resolution parameter affects cluster number; adjust according to data characteristics
3. **Marker Quantity**: Default uses top 10 markers; adjust as needed
4. **Query Parameters**: `-n` parameter controls returned result quantity; increase appropriately for complex data
5. **Ranking Metrics**: brick rk supports multiple ranking metrics; combine according to needs
6. **Result Validation**: Validate annotation results using domain knowledge
