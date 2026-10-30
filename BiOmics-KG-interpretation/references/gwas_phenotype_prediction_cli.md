# GWAS Phenotype Prediction Analysis with brick CLI

## Task Overview

GWAS phenotype prediction analysis aims to systematically predict the phenotypic spectrum that these genetic variants may affect through multi-dimensional association analysis of the knowledge graph, starting from significant loci of genome-wide association studies. This tutorial demonstrates how to use brick CLI tools combined with gene annotation to complete the full workflow from SNP extraction, gene mapping to multi-type entity association query and ranking enrichment.

**Note: This workflow is for reference only. Please adjust according to user needs and real-time conditions.**

## Analysis Workflow

```
Extract significant SNPs ->Gene annotation ->brick qr multi-type query ->Merge results ->brick rk ranking/enrichment ->Visualization and interpretation
```

## Prerequisites

- BiOmics-BRICK package installed: `pip install BiOmics-BRICK`
- pandas, matplotlib, networkx installed (for data processing and visualization)
- GWAS summary statistics prepared (at least containing chromosome, genomic position, statistical significance columns)

## Step-by-step Guide

### Step 1: Extract Significant SNPs

Filter genome-wide significant loci (typically p < 5e-8) from GWAS summary data, and extract Top SNPs by chromosome for subsequent gene annotation.

> **Note**: The following Python code is for demonstration only. Please adjust actual SNP filtering strategy (e.g., LD clumping, independent signal detection, etc.) according to data characteristics and research needs.

```python
import pandas as pd

# Load GWAS summary data
df = pd.read_csv("your_gwas_data.csv")

# Extract Top N SNPs per chromosome
chromosomes = sorted(df['CHROM'].unique())
top_snps_list = []
for chrom in chromosomes:
    chrom_data = df[df['CHROM'] == chrom].sort_values('LOG10P', ascending=False)
    top_n = min(3, len(chrom_data))  # Top 3 per chromosome
    top_snps_list.append(chrom_data.head(top_n))

top_snps = pd.concat(top_snps_list).sort_values('LOG10P', ascending=False)
top_snps.to_csv('top_snps_per_chr.csv')
```

### Step 2: Gene Annotation

Annotate nearby genes for extracted Top SNPs, mapping genomic coordinates to gene symbols. MyGene.info online API is recommended; professional annotation tools like ANNOVAR can also be used.

> **Note**: Annotation window size (default 50kb) can be adjusted according to linkage disequilibrium decay distance; prioritize protein-coding genes. The following code demonstrates basic MyGene.info calling; adjust request frequency according to data scale.

```python
import requests
import time

def annotate_gene(chrom, pos, window=50000):
    """Query MyGene.info for genes near SNP"""
    url = 'https://mygene.info/v3/query'
    params = {
        'q': f'chr{chrom}:{pos-window}-{pos+window}',
        'species': 'human',
        'fields': 'symbol,name,genomic_pos',
        'size': 5
    }
    try:
        resp = requests.get(url, params=params, timeout=15)
        if resp.status_code == 200:
            hits = resp.json().get('hits', [])
            return [h.get('symbol', '') for h in hits if h.get('symbol')]
    except Exception as e:
        print(f"Query failed chr{chrom}:{pos}: {e}")
    return []

# Annotate each Top SNP
all_genes = set()
for _, row in top_snps.iterrows():
    genes = annotate_gene(int(row['CHROM']), int(row['GENPOS']))
    all_genes.update(genes)
    time.sleep(0.3)  # Avoid too many requests

# Save candidate gene list (for brick qr)
gene_list = sorted(all_genes)
with open('candidate_genes.txt', 'w') as f:
    f.write('\n'.join(gene_list))
```

### Step 3: brick qr Multi-type Query

Use `brick qr` to perform multi-dimensional knowledge graph queries on candidate genes, covering Disease, Phenotype, Pathway, Chemical, Process, and other entity types.

> **Note**: When the number of genes is large or there are many query types, it is recommended to batch queries to avoid timeout from overly large single queries.

```bash
# Define target entity type list
target_types=(Disease Phenotype Pathway Chemical Process)

# Loop query each entity type
for tt in "${target_types[@]}"; do
    brick qr -i candidate_genes.txt \
             --st Gene \
             --tt "$tt" \
             -n 500 \
             -o "query_${tt,,}.csv"
done
```

**Parameter Description:**
- `-i`: Candidate gene list file (one gene symbol per line)
- `--st`: Source entity type (Gene)
- `--tt`: Target entity type
- `-n`: Maximum number of returned results
- `-o`: Output CSV file

### Step 4: Merge Query Results

Merge query results from all entity types into a single file for subsequent ranking and enrichment analysis.

> **Note**: The following code demonstrates basic data merging logic. Please adjust according to file naming and data structure.

```python
import pandas as pd
import glob

# Read all query results
csv_files = glob.glob('query_*.csv')
dfs = [pd.read_csv(f) for f in csv_files if pd.read_csv(f).shape[0] > 0]

# Merge
combined = pd.concat(dfs, ignore_index=True)
combined.to_csv('combined_query.csv', index=False)
print(f"Merge complete: {len(combined)} association results")
```

### Step 5: Ranking and Enrichment Analysis

Use `brick rk` to perform ranking and enrichment analysis on merged query results, identifying phenotypes and pathways most strongly associated with the candidate gene set.

#### 5.1 Ranking Analysis

```bash
brick rk -i combined_query.csv --rank -o ranked_results.csv
```

**Key Output Fields:**
- `path.2.name`: Target entity name
- `path.2.type`: Target entity type
- `path.2.rank_voting`: Ranking voting score (lower = more significant)
- `path.2.match_count`: Number of matched source genes
- `path.2.info_source_count`: Number of information sources (literature evidence volume)

#### 5.2 Enrichment Analysis

```bash
brick rk -i combined_query.csv --enrich -o enriched_results.csv
```

**Key Output Fields:**
- `path.2.enrich_pvalue`: Enrichment p-value (smaller = more significant)
- `path.2.enrich_fdr`: FDR-corrected p-value
- `path.2.enrich_statistic`: Enrichment statistic

### Step 6: Visualization

Generate visualization charts based on ranking and enrichment results, showing the association network between candidate genes and predicted phenotypes.

> **Note**: The following Python code is for visualization demonstration. Please adjust plotting parameters (colors, thresholds, layout, etc.) according to actual data distribution and presentation needs.

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ── Ranking bar chart ──
rank_df = pd.read_csv('ranked_results.csv')
top_rank = rank_df.nsmallest(20, 'path.2.rank_voting')

fig, ax = plt.subplots(figsize=(10, 9))
ax.barh(range(len(top_rank)), top_rank['path.2.rank_voting'].values)
ax.set_yticks(range(len(top_rank)))
ax.set_yticklabels(top_rank['path.2.name'], fontsize=9)
ax.set_xlabel('Rank Voting Score (lower = more significant)')
ax.invert_yaxis()
plt.tight_layout()
plt.savefig('ranking_barplot.png', dpi=150)

# ── Enrichment bar chart ──
enrich_df = pd.read_csv('enriched_results.csv')
sig_enrich = enrich_df[enrich_df['path.2.enrich_pvalue'] < 0.05]
top_enrich = sig_enrich.nsmallest(15, 'path.2.enrich_pvalue')
top_enrich['-log10(p)'] = -np.log10(top_enrich['path.2.enrich_pvalue'].clip(lower=1e-300))

fig, ax = plt.subplots(figsize=(10, 7))
ax.barh(range(len(top_enrich)), top_enrich['-log10(p)'].values)
ax.set_yticks(range(len(top_enrich)))
ax.set_yticklabels(top_enrich['path.2.name'], fontsize=9)
ax.set_xlabel('-log10(enrichment p-value)')
ax.invert_yaxis()
plt.tight_layout()
plt.savefig('enrichment_barplot.png', dpi=150)
```

### Step 7: Build Association Network

Use networkx to build a gene-phenotype association network, with the candidate gene set as the central node and associated entities as peripheral nodes.

```python
import networkx as nx
import matplotlib.pyplot as plt

# Read ranking results, take Top N for network construction
rank_df = pd.read_csv('ranked_results.csv')
top_network = rank_df.nsmallest(15, 'path.2.rank_voting')

G = nx.Graph()
center_node = 'Candidate Genes'
G.add_node(center_node)

for _, row in top_network.iterrows():
    target = row['path.2.name']
    G.add_node(target)
    G.add_edge(center_node, target)

fig, ax = plt.subplots(figsize=(12, 10))
pos = nx.spring_layout(G, k=0.5)
nx.draw(G, pos, with_labels=True, node_size=2000, font_size=9, edge_color='gray', ax=ax)
ax.set_title('GWAS Phenotype Prediction Network')
plt.tight_layout()
plt.savefig('phenotype_network.png', dpi=150)
```

## Output Files

| File | Description |
|------|-------------|
| `top_snps_per_chr.csv` | Top SNP list per chromosome |
| `candidate_genes.txt` | Candidate gene list (brick qr input) |
| `query_disease.csv` | Disease type raw query results |
| `query_phenotype.csv` | Phenotype type raw query results |
| `query_pathway.csv` | Pathway type raw query results |
| `query_chemical.csv` | Chemical type raw query results |
| `query_process.csv` | Process type raw query results |
| `combined_query.csv` | Merged multi-dimensional query results |
| `ranked_results.csv` | brick rk ranking results |
| `enriched_results.csv` | brick rk enrichment results |
| `ranking_barplot.png` | Ranking analysis bar chart |
| `enrichment_barplot.png` | Enrichment analysis bar chart |
| `phenotype_network.png` | Gene-phenotype association network |

## Result Interpretation

1. **Ranking Analysis**: The lower the `path.2.rank_voting` score, the stronger the association between the entity and candidate gene set. Pay attention to top-ranked entities with high `path.2.info_source_count`, as the number of information sources reflects the adequacy of literature evidence.

2. **Enrichment Analysis**: The smaller the `path.2.enrich_pvalue`, the more significant the enrichment of the entity in the candidate gene set. It is recommended to filter results with p < 0.05 for focused interpretation.

3. **Pleiotropy Identification**: When Top results from multiple types (Disease, Phenotype, Pathway) point to the same biological theme (e.g., metabolism, inflammation, tumor), it supports the inference of genetic pleiotropy.

4. **Network Centrality**: Entities with high connectivity to the central node in the network graph are core predicted phenotypes; causal hypotheses of gene-phenotype relationships can be constructed accordingly.

## Notes

1. **Gene Annotation Tool Selection**: MyGene.info online API is convenient but network-dependent; ANNOVAR, VEP, and other local tools are more stable but require configuring reference genomes; pyensembl and other Python libraries can also be used
2. **Annotation Window**: Default 50kb window is suitable for most cases, but regulatory elements may be farther away; adjust according to research needs
3. **Batch Querying**: When the number of candidate genes is large, it is recommended to batch query different entity types to avoid server timeout
4. **LD Deduplication**: Multiple highly significant SNPs in strong linkage disequilibrium regions may annotate to the same gene; pay attention to deduplication during interpretation
5. **Unannotated Loci**: Some highly significant loci may not annotate to known genes, possibly involving non-coding RNAs or long-distance regulatory elements, worthy of further functional validation
6. **Result Validation**: Phenotype prediction results should be validated with domain knowledge and independent cohort data to avoid over-interpreting co-occurrence relationships

For direct Mutation queries, store exactly one canonical rsID per line and use `-i`. Preserve complete HGVS literals—including spaces, colons, commas, alleles, and parentheses—in the upstream rsID-mapping input. The CLI rejoins adjacent unquoted words, but quoting or file input remains the unambiguous choice; for disease names containing commas, use `--target-file` rather than a comma-separated `-t` value.

```bash
brick qr -i rsids.txt --st Mutation --tt "Disease|Phenotype" -n 500 -o mutation_neighbors.csv
```

