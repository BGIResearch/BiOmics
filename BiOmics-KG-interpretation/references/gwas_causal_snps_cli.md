# GWAS Causal SNP Analysis with brick CLI

## Task Overview

GWAS causal SNP analysis aims to identify SNP loci with potential causal effects from genome-wide association study summary statistics, and verify their biological associations with target diseases/phenotypes through the knowledge graph. This tutorial demonstrates how to combine brick CLI tools to complete the full workflow from SNP extraction, gene annotation to knowledge graph association analysis.

**Note: This workflow is for reference only. Please adjust according to user needs and real-time conditions.**

## Analysis Workflow

```
GWAS summary data ->Extract Top SNP ->Gene annotation ->brick qr query association ->brick rk ranking/enrichment ->Visualization and interpretation
```

## Prerequisites

- BiOmics-BRICK package installed: `pip install BiOmics-BRICK`
- GWAS summary statistics prepared (at least containing chromosome, genomic position, p-value or -log10(p) columns)
- matplotlib, pandas, numpy installed (for visualization and data processing)

## Step-by-step Guide

### Step 1: Extract Significant SNPs

Filter genome-wide significant loci (typically p < 5e-8) from GWAS summary data, and extract Top SNPs by chromosome.

> **Note**: The following Python code is for demonstration only. Please adjust actual SNP filtering strategy (e.g., LD clumping, independent signal detection, etc.) according to data characteristics and research needs.

```python
import pandas as pd
import numpy as np

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

### Step 2: Generate Manhattan Plot

> **Note**: The following Python code is for visualization demonstration. Please adjust plotting parameters (colors, labels, threshold lines, etc.) according to actual data distribution and presentation needs.

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(16, 5))
colors = ['#1f77b4', '#d62728']

# Compute chromosome offsets for horizontal arrangement
chrom_offset = {}
offset = 0
for chrom in chromosomes:
    chrom_data = df[df['CHROM'] == chrom]
    chrom_offset[chrom] = offset - chrom_data['GENPOS'].min()
    offset += chrom_data['GENPOS'].max() - chrom_data['GENPOS'].min() + 5e6

for i, chrom in enumerate(chromosomes):
    chrom_data = df[df['CHROM'] == chrom]
    pos = chrom_data['GENPOS'] + chrom_offset[chrom]
    ax.scatter(pos, chrom_data['LOG10P'], c=colors[i % 2], s=8, alpha=0.6)

ax.axhline(y=-np.log10(5e-8), color='grey', linestyle='--', label='p=5e-8')
ax.set_ylabel('-log10(P-value)')
ax.set_title('GWAS Manhattan Plot')
plt.tight_layout()
plt.savefig('manhattan_plot.png', dpi=150)
```

### Step 3: Gene Annotation

Annotate nearby genes for extracted Top SNPs, mapping genomic coordinates to gene symbols. MyGene.info online API is recommended; local tools such as pyensembl can also be used.

> **Note**: Annotation window size (default 50kb) can be adjusted according to linkage disequilibrium decay distance; prioritize protein-coding genes.

```python
import requests
import time

def annotate_gene_mygene(chrom, pos, window=50000):
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
    genes = annotate_gene_mygene(int(row['CHROM']), int(row['GENPOS']))
    all_genes.update(genes)
    time.sleep(0.3)  # Avoid too many requests

# Save candidate gene list (for brick qr)
gene_list = sorted(all_genes)
with open('candidate_genes.txt', 'w') as f:
    f.write('\n'.join(gene_list))
```

### Step 4: Query brick Knowledge Graph

Use `brick qr` to query associations between candidate genes and diseases. Two query modes are supported:

#### 4.1 Neighbor Mode: Query All Diseases Associated with Genes

```bash
brick qr -i candidate_genes.txt \
    --st Gene \
    --tt Disease \
    -n 200 \
    -o genes_diseases_query.csv
```

**Parameter Description:**
- `-i`: Candidate gene list file (one gene symbol per line)
- `--st Gene`: Source entity type is Gene
- `--tt Disease`: Target entity type is Disease
- `-n`: Maximum number of returned results
- `-o`: Output CSV file

#### 4.2 Target Mode: Query Direct Association Between Genes and Specific Disease

If the research target is clear (e.g., Type 2 diabetes), you can directly specify the target entity:

```bash
brick qr "GENE1,GENE2,GENE3" \
    -t "your_target_disease" \
    --st Gene \
    --tt Disease \
    -n 200 \
    -o genes_disease_direct.csv
```

**Parameter Description:**
- Positional parameter: comma-separated gene list
- `-t`: Target entity name (only returns relationships with this entity when specified)
- Other parameters same as above

### Step 5: Ranking and Enrichment Analysis

Use `brick rk` to perform ranking and enrichment analysis on query results, identifying the most significant disease associations.

#### 5.1 Ranking Analysis

```bash
brick rk -i genes_diseases_query.csv --rank -o genes_diseases_ranked.csv
```

**Key Output Fields:**
- `path.2.name`: Target entity name
- `path.2.rank_voting`: Ranking voting score (lower = more significant)
- `path.2.rank_information_source`: Number of information sources

#### 5.2 Enrichment Analysis

```bash
brick rk -i genes_diseases_query.csv --enrich -o genes_diseases_enriched.csv
```

**Key Output Fields:**
- `path.2.enrich_pvalue`: Enrichment p-value (smaller = more significant)
- `path.2.enrich_fdr`: FDR-corrected p-value
- `path.2.enrich_statistic`: Enrichment statistic

### Step 6: Visualization

Use `brick pl` to generate visualizations, or use Python to plot directly.

#### Method A: brick pl Radial Plot

```bash
brick pl radial -i genes_diseases_enriched.csv --top-k 10 -o enrichment_radial.png
```

**Parameter Description:**
- `-i`: Input CSV file (enrichment analysis results)
- `--top-k`: Display top K entities (default 7)
- `-o`: Output image file (supports .png/.pdf)


## Output Files

| File | Description |
|------|-------------|
| `top_snps_per_chr.csv` | Top SNP list per chromosome |
| `manhattan_plot.png` | Genome-wide Manhattan plot |
| `candidate_genes.txt` | Candidate gene list (brick qr input) |
| `genes_diseases_query.csv` | brick qr neighbor mode query results |
| `genes_disease_direct.csv` | brick qr target mode query results (direct associations) |
| `genes_diseases_ranked.csv` | brick rk ranking results |
| `genes_diseases_enriched.csv` | brick rk enrichment results |
| `enrichment_barplot.png` | Enrichment analysis visualization |

## Result Interpretation

1. **Manhattan Plot**: Identify significant signal regions across the genome; pay attention to strong signal loci far exceeding the genome-wide significance threshold
2. **Ranking Analysis**: The lower the ranking voting score, the stronger the association between the disease and candidate gene set; more information sources indicate more sufficient evidence
3. **Enrichment Analysis**: The smaller the p-value, the more significant the enrichment of the disease in the candidate gene set; focus on results with FDR < 0.05
4. **Direct Association**: Target mode queries that reveal `functionally_related_to` or `causally_related_to` gene-disease pairs are the strongest knowledge graph evidence
5. **Multi-locus Synergy**: Genes on different chromosomes may jointly point to the target disease through different pathways (e.g., Wnt signaling, insulin secretion, cAMP), forming a multi-pathway genetic susceptibility network

## Notes

1. **Gene Annotation Tool Selection**: MyGene.info online API is convenient but network-dependent; pyensembl runs locally more stably but requires downloading data packages; professional annotation tools like ANNOVAR can also be used
2. **Annotation Window**: Default 50kb window is suitable for most cases, but regulatory elements may be farther away; adjust according to research needs
3. **Query Entity Types**: In addition to Disease, Gene-Cell, Gene-Pathway, and other relationships can also be queried; choose flexibly according to research objectives
4. **Batch Querying**: When the number of candidate genes is large, batch query different target entity types to avoid excessively large single queries
5. **LD Deduplication**: Multiple highly significant SNPs in strong linkage disequilibrium regions may annotate to the same gene; pay attention to deduplication during interpretation

#### 4.0 Preserve SNP and Disease Literals

Use newline-delimited files for Mutation queries. Put exactly one canonical rsID on each line, and put a disease name containing commas in its own target file. Keep complete HGVS literals intact in the upstream rsID-mapping input. This avoids shell word splitting and prevents punctuation inside identifiers from being treated as formatting.

`rsids.txt`:

```text
rs7756992
rs4506565
```

`t2d_target.txt`:

```text
Diabetes Mellitus, Type 2
```

```bash
brick qr -i rsids.txt \
    --target-file t2d_target.txt \
    --st Mutation \
    --tt "Disease|Phenotype" \
    -n 200 \
    -o snp_t2d_query.csv
```

The current CLI rejoins adjacent unquoted source words, but files or shell quoting remain preferred because they preserve the intended entity boundary unambiguously. If the starting value is an HGVS genomic position rather than a canonical rsID, preserve the full colon-containing literal during the mapping step and map it to an rsID before a name-based `brick qr` Mutation query; do not shorten it to the reference accession prefix.
