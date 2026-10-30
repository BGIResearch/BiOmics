# brick CLI Gene Enrichment Analysis Tutorial

## Task Overview

Gene enrichment analysis is used to interpret the functional significance of differentially expressed genes (DEGs). Gene sets obtained from single-cell analysis are input into brick to query their relationships with various entities (cell types, diseases, pathways, etc.) in the knowledge graph. Significantly associated biological entities are identified through enrichment analysis and displayed via radial plots.

## Analysis Workflow

```
Extract gene set ->Query brick ->Enrichment analysis ->Filter significant results ->Visualization
```

## Prerequisites

- brick installed: `pip install BiOmics-BRICK`
- Single-cell data file (.h5ad format, containing rank_genes_groups results)

## Step-by-step Guide

### Step 1: Extract Gene Set

Extract Top N genes for the target cell type from single-cell analysis results.

```python
import scanpy as sc

# Read data
adata = sc.read('your_data.h5ad')
adata.X = adata.layers['lognorm'].copy()

# Compute DEGs (if not already done)
sc.tl.rank_genes_groups(adata, groupby='your_annotation_column', pts=True)

# Extract Top 100 genes
cell_type = 'your_cell_type'
top_genes = list(adata.uns['rank_genes_groups']['names'][cell_type])[:100]

# Save to file
with open('gene_list.txt', 'w') as f:
    f.write('\n'.join(top_genes))
```

### Step 2: Query brick

Use `brick qr` to batch query relationships between genes and various entity types (avoid overly large single queries).

```bash
# Define target entity type list
target_types=(Cell Tissue Disease Phenotype Chemical Process Function Cell_Component Pathway)

# Define relationship types
relations="produce,participate_in,causally_related_to,capable_of,spatially_related_to,processly_related_to,functionally_related_to,disease_relation,pharmaco_effect_on"

# Loop query each entity type
for tt in "${target_types[@]}"; do
    brick qr -i gene_list.txt \
             --st Gene \
             --tt "$tt" \
             -r "$relations" \
             -n 1000 \
             -o "query_${tt,,}.csv"
done
```

**Parameter Description:**
- `-i`: Gene list file
- `--st`: Source entity type (Gene)
- `--tt`: Target entity type
- `-r`: Relationship types (comma-separated)
- `-n`: Maximum number of returned results
- `-o`: Output CSV file

### Step 3: Merge Query Results

Merge query results from all entity types into a single file.

```python
import pandas as pd
import glob

# Read all query results
csv_files = glob.glob('query_*.csv')
dfs = [pd.read_csv(f) for f in csv_files if pd.read_csv(f).shape[0] > 0]

# Merge
combined = pd.concat(dfs, ignore_index=True)
combined.to_csv('combined_query.csv', index=False)
```

### Step 4: Enrichment Analysis

Use `brick rk` to perform enrichment analysis on query results.

```bash
brick rk -i combined_query.csv --enrich -o enrichment_results.csv
```

**Output Column Description:**
- `path.2.name`: Target entity name
- `path.2.type`: Target entity type
- `path.2.enrich_pvalue`: Enrichment p-value (smaller = more significant)
- `path.2.enrich_statistic`: Enrichment statistic
- `path.2.enrich_fdr`: FDR-corrected p-value

### Step 5: Filter Significant Results

Filter significantly enriched entities with p < 0.05.

```python
import pandas as pd

enrich_df = pd.read_csv('enrichment_results.csv')
significant_df = enrich_df[enrich_df['path.2.enrich_pvalue'] < 0.05]
significant_df.to_csv('significant_enrichment.csv', index=False)

print(f"Significantly enriched entities: {len(significant_df)}")
```

### Step 6: Generate Radial Enrichment Plot

Use `brick pl radial` to generate radial enrichment visualization.

```bash
brick pl radial -i significant_enrichment.csv -o radial_enrichment.pdf --top-k 15
```

**Parameter Description:**
- `-i`: Input CSV file (enrichment analysis results)
- `-o`: Output image file (supports .png/.pdf)
- `--top-k`: Display top K entities (default 7)


## Output Files

| File | Description |
|------|-------------|
| `gene_list.txt` | Input Top N gene list |
| `query_*.csv` | Raw query results for each entity type |
| `combined_query.csv` | Merged query results |
| `enrichment_results.csv` | Complete enrichment analysis results |
| `significant_enrichment.csv` | Significant enrichment results (p < 0.05) |
| `radial_enrichment.pdf` | Radial enrichment visualization |

## Result Interpretation

The radial enrichment plot displays Top K significantly enriched entities in polar coordinates:
- **Radius**: Represents enrichment significance (-log10(p-value)), larger radius = more significant
- **Angle**: Grouped by entity type
- **Color**: Distinguishes different entity types (Cell/Process/Disease/Tissue)

## Notes

1. **Batch Querying**: When the number of genes is large, it is recommended to batch query different entity types to avoid server timeout
2. **P-value Threshold**: Default significance threshold is p < 0.05; adjust as needed
3. **Top-K Selection**: It is recommended to set `--top-k` to 10-20 for radial plots; too many will cause label overlap
4. **Memory Management**: Pay attention to memory usage when merging large query results; process in batches if necessary
