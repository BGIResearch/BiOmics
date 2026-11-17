# brick CLI Gene Regulatory Network Analysis Tutorial

## Task Overview

Gene regulatory network (GRN) analysis is used to resolve regulatory relationships between transcription factors (TFs) and their target genes (Regulons). This tutorial uses brick CLI commands as the core, combined with minimal Python auxiliary processing, to complete the full workflow from pySCENIC output to TF regulatory network visualization.

## Analysis Workflow

```
pySCENIC output ──->Data preprocessing (Python auxiliary) ──->brick pl tf visualization ──->brick qr query ──->Result interpretation
```

## Prerequisites

- brick installed: `pip install BiOmics-BRICK`
- pySCENIC analysis completed `_modules.tsv` file

---

## Step 1: Prepare brick pl tf Input Data

`brick pl tf` requires a CSV file in a specific format, containing a `Regulon` column (with gene lists as string representations).

### Input Format Requirements

```csv
TF,Regulon,IntersectionSize
Fos,"['Gene1', 'Gene2', 'Gene3']",3
Stat2,"['Gene4', 'Gene5']",2
```

### Auxiliary Processing Script (Example)

The following Python script converts pySCENIC `_modules.tsv` to the standard input format for brick pl tf:

```python
import pandas as pd
import re
import os

# Configure input/output paths
module_file = "your_modules.tsv"   # pySCENIC output module file
output_dir = "grn_output"
os.makedirs(output_dir, exist_ok=True)

def extract_genes(target_genes_str):
    """Extract gene names from pySCENIC format"""
    return set(re.findall(r"'([^']+)'", target_genes_str))

# Read pySCENIC output (skip first 2 comment lines)
df = pd.read_csv(module_file, sep="\t", skiprows=2)

# Locate target gene column and extract (column name may be Unnamed: 8, etc.; adjust according to actual file)
target_col = [c for c in df.columns if 'Gene' in c or 'Unnamed' in c][-1]
df = df[["TF", target_col]].dropna()
df["Genes"] = df[target_col].apply(extract_genes)

# Merge multiple motifs for the same TF, deduplicate
tf_dict = {}
for _, row in df.iterrows():
    tf_dict.setdefault(row["TF"], set()).update(row["Genes"])

# Generate standard input format
records = []
for tf, genes in tf_dict.items():
    records.append({
        "TF": tf,
        "Regulon": str(list(genes)),
        "IntersectionSize": len(genes)
    })

out_df = pd.DataFrame(records)
out_file = os.path.join(output_dir, "tf_regulon_input.csv")
out_df.to_csv(out_file, index=False)
print(f"Input file generated: {out_file}")
```

**Key Points:**
- The `Regulon` column must be a string representation of a Python list (e.g., `"['Gene1', 'Gene2']"`), not an actual list object, otherwise brick pl tf will parse the entire string as a single node
- The first column is the TF name, used as an index
- Multiple motifs for the same TF need to be deduplicated and merged

---

## Step 2: Visualize TF Regulatory Network

Use the `brick pl tf` command to plot the regulatory network between TFs and their target genes.

### Basic Usage

```bash
brick pl tf -i tf_regulon_input.csv --tf "TF_NAME" -o tf_network.png
```

**Parameter Description:**

| Parameter | Short | Description | Required |
|-----------|-------|-------------|----------|
| `--input` | `-i` | TF-Regulon data CSV (containing Regulon column) | Yes |
| `--tf` | - | TF name, comma-separated for multiple | Yes |
| `--verified` | - | Verified genes file (one gene per line) | No |
| `--output` | `-o` | Output image path (default tf_regulon.png) | No |

### Example Commands

```bash
# Single TF visualization
brick pl tf -i tf_regulon_input.csv --tf "Fos" -o Fos_regulon.png

# Multiple TFs simultaneous visualization
brick pl tf -i tf_regulon_input.csv --tf "Fos,Stat2,Maf" -o core_tf_network.png

# With verified gene highlighting
brick pl tf -i tf_regulon_input.csv --tf "Fos" --verified verified_genes.txt -o Fos_regulon.png
```

### Generate Verified Genes File (Optional)

Verified genes are regulatory relationship genes supported by literature or databases, obtainable through brick knowledge graph queries:

```bash
# Query known relationships between TF and other genes
brick qr "TF_NAME" \
    --st Gene \
    --tt Gene \
    -r "regulated,interact_with" \
    -n 500 \
    -o tf_gene_relations.csv
```

Then extract target gene names and save as a verification file:

```python
import pandas as pd

df = pd.read_csv("tf_gene_relations.csv")
genes = df["path.2.name"].dropna().unique()
with open("verified_genes.txt", "w") as f:
    f.write("\n".join(genes))
```

---

## Step 3: Query TF-related Diseases and Functions

Use `brick qr` to query associations between TFs and diseases, functions, and other entities to obtain biological background information.

```bash
# Query TF-related diseases
brick qr "TF_NAME" \
    --st Gene \
    --tt Disease \
    -r "causally_related_to,functionally_related_to" \
    -n 20 \
    -o tf_diseases.csv

# Query TF-related biological processes
brick qr "TF_NAME" \
    --st Gene \
    --tt Process \
    -r "participate_in" \
    -n 20 \
    -o tf_processes.csv
```

---

## Step 4: Summarize Analysis Results (Optional)

A simple Python script can be used to summarize analysis results:

```python
import pandas as pd

tf_name = "TF_NAME"
regulon_df = pd.read_csv("tf_regulon_input.csv")
tf_row = regulon_df[regulon_df["TF"] == tf_name].iloc[0]
regulon_size = tf_row["IntersectionSize"]

# Read verification results (if available)
try:
    verified = pd.read_csv("tf_gene_relations.csv")
    verified_count = verified["path.2.name"].nunique()
except:
    verified_count = 0

print(f"TF: {tf_name}")
print(f"Regulon size: {regulon_size}")
print(f"Verified gene count: {verified_count}")
```

---

## Output Files

| File | Source | Description |
|------|--------|-------------|
| `tf_regulon_input.csv` | Python preprocessing | Standard input for brick pl tf |
| `*_regulon.png` | brick pl tf | TF regulatory network visualization |
| `tf_gene_relations.csv` | brick qr | Known relationships between TF and genes |
| `tf_diseases.csv` | brick qr | TF-related disease query results |
| `verified_genes.txt` | Python auxiliary | Verified gene list (optional) |

---

## Result Interpretation

In the TF regulatory network plot:
- **Pink square**: Target transcription factor (TF)
- **Grey circles**: Regulon target genes of this TF
- **Green circles** (if present): Verified regulatory relationship genes
- **Lines**: Regulatory relationships between TF and target genes

---

## Notes

1. **Regulon Column Format**: The `Regulon` column must be a **string representation** of a Python list (e.g., `"['Gene1', 'Gene2']"`), not an actual list object, otherwise brick pl tf will parse the entire string as a single node

2. **pySCENIC Column Name Differences**: Different versions of pySCENIC may have different output column names; the target gene column may be `TargetGenes` or `Unnamed: 8`, etc. Adjust according to actual file

3. **Gene Name Case**: The brick knowledge graph typically uses uppercase gene names. If input is lowercase (e.g., mouse gene symbols), verified gene matching may be empty. It is recommended to standardize to uppercase

4. **Multiple TF Visualization**: The `--tf` parameter supports comma separation (e.g., `--tf "TF1,TF2,TF3"`), suitable for comparing core TF regulatory networks

5. **Verified Genes File**: The `--verified` parameter accepts a plain text file with one gene name per line; successfully matched genes will be highlighted in green
