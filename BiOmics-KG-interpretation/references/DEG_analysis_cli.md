# Differential Expression Gene Analysis with brick CLI

## Task Overview

Differential expression gene (DEG) analysis is a core step in single-cell RNA sequencing data analysis. This tutorial demonstrates how to use brick CLI tools combined with scanpy to perform DEG analysis, and compare known markers from the brick knowledge graph with experimental data through visualization.

**Note: This workflow is for reference only. Please adjust according to user needs and real-time conditions.**

## Prerequisites

- BiOmics-BRICK package installed: `pip install BiOmics-BRICK`
- scanpy installed: `pip install scanpy`
- AnnData format single-cell data prepared (containing expression matrix and cell type annotations)

## Analysis Workflow

### Step 1: Compute DEGs

Use scanpy to perform differential expression analysis on the target cell type to obtain logFC and p-values.

> **Note**: The following Python code is for demonstration only. Please adjust the computation method according to your data characteristics and research needs.

```python
import scanpy as sc
import pandas as pd
import numpy as np

# Load data and compute DEGs
adata = sc.read("./your_data.h5ad")
adata.X = adata.layers['lognorm'].copy()
sc.tl.rank_genes_groups(adata, groupby='celltype_annotation', pts=True, use_raw=False)

# Extract DEG results for the target cell type
cell_type = "your_cell_type"
df = pd.DataFrame({
    'gene': adata.uns['rank_genes_groups']['names'][cell_type],
    'pvals': adata.uns['rank_genes_groups']['pvals'][cell_type],
    'logfoldchanges': adata.uns['rank_genes_groups']['logfoldchanges'][cell_type]
})
df['x'] = df['logfoldchanges']
df['y'] = -1 * np.log10(df['pvals'])
```

### Step 2: Query brick Markers

Use `brick qr` to retrieve known markers for the target cell type from the knowledge graph.

```bash
brick qr "your_cell_type" \
    --st Cell \
    --tt Gene \
    -r produce \
    -n 500 \
    -o ./brick_markers.csv
```

Associate brick query results with DEG data, filtering markers present in the dataset and assigning logFC and p-values:

> **Note**: The following code demonstrates basic data merging logic. Please adjust according to your data structure.

```python
query_df = pd.read_csv("./brick_markers.csv")
query_df = query_df[query_df['path.2.name'].isin(adata.var_names)]
query_df.index = query_df['path.2.name']
for key in ['x', 'y']:
    query_df[key] = df.set_index('gene')[key]
query_df = query_df.dropna(subset=['x', 'y'])
```

### Step 3: Generate Focused Volcano Plot

Create a volcano plot focused on the brick markers region to avoid label overlap in the full plot.

> **Note**: The following Python code is for visualization demonstration. Please adjust plotting parameters (colors, thresholds, label density, etc.) according to actual data distribution and presentation needs.

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12, 10))

# Focus on the brick markers region (adjust based on data distribution)
focus_x_min = -1
focus_x_max = query_df['x'].max() + 1
focus_y_max = query_df['y'].max() + 2

# Plot background genes
focus = df[(df['x'] >= focus_x_min) & (df['x'] <= focus_x_max)]
ax.scatter(focus['x'], focus['y'], s=5, color='lightgrey', alpha=0.6)

# Distinguish significant and non-significant markers
sig = query_df[(query_df['y'] > -np.log10(0.05)) & (abs(query_df['x']) > 0.5)]
ns = query_df[~query_df.index.isin(sig.index)]
ax.scatter(ns['x'], ns['y'], s=50, facecolors='none', edgecolors='coral', linewidths=1.5)
ax.scatter(sig['x'], sig['y'], s=150, facecolors='lightcoral', edgecolors='darkred', linewidths=2)

# Add gene labels (alternating placement to avoid overlap)
for i, (gene, row) in enumerate(sig.sort_values('y', ascending=False).iterrows()):
    offset = 0.1 if i % 2 == 0 else -0.1
    ha = 'left' if i % 2 == 0 else 'right'
    ax.text(row['x'] + offset, row['y'], gene, fontsize=9, ha=ha,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

ax.axhline(y=-np.log10(0.05), color='blue', linestyle='--', alpha=0.5)
ax.set_xlabel('Log2 Fold Change')
ax.set_ylabel('-Log10(p-value)')
ax.set_title(f'brick Markers: {cell_type}')
ax.set_xlim([focus_x_min, focus_x_max])
ax.set_ylim([0, focus_y_max])
ax.legend()
plt.tight_layout()
plt.savefig("./DEG_volcano.png", dpi=150)
```

## Visualization Key Points

1. **Focused Display**: Only show the region where brick markers are located, avoiding wasted space on irrelevant regions
2. **Layered Display**:
   - Significant markers (dark red filled circles)
   - Non-significant markers (light red hollow circles)
   - Other genes (light grey dots)
3. **Smart Labels**: Alternating left-right placement of gene names with white background boxes for clarity
4. **Threshold Line**: Blue dashed line indicating the p = 0.05 threshold

## Notes

1. **Data Adaptation**: Adjust annotation column names, cell type names, and other parameters according to actual data
2. **Focus Region**: Adjust `focus_x_min` and `focus_x_max` based on actual brick markers distribution
3. **Label Density**: If there are too many markers, only label the top N most significant ones
4. **Result Validation**: Validate consistency between brick markers and experimental data using domain knowledge
