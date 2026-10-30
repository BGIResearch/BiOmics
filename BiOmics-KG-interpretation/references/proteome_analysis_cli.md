# Proteomics Differential and PTM-Aware PPI Analysis with brick

Use this reference to analyze bulk proteomics or peptide/PTM-feature data, select an auditable candidate set, query protein interactions from BiOmics-KG, and build a traceable network. Keep statistical evidence, KG evidence, and observational expression correlation separate.

## Contents

1. Analysis contract
2. Environment and required input
3. Audit data and preserve feature identity
4. Run differential analysis
5. Freeze protein or PTM-state candidates
6. Query an induced PPI subgraph
7. Review and deduplicate KG evidence
8. Build an attributed network
9. Render interactive and publication figures
10. Reproduce legacy Figure 7 deliberately
11. Interpret, validate, and report results

## 1. Analysis Contract

Freeze these decisions before inspecting differential rankings or KG results:

- Contrast, group column, numerator, denominator, and sample counts.
- Whether `adata.X` contains raw non-negative intensity, normalized intensity, or an already transformed matrix.
- Zero and missing-value handling.
- Statistical test, multiple-testing correction, effect-size definition, and primary threshold.
- Any exploratory fallback, including its trigger and label.
- Feature identity and feature-to-protein aggregation.
- Candidate cap and tie-breaking rule.
- KG relation, hop count, direction, target restriction, and evidence-source policy.
- Network node mode: one node per protein or one node per `(protein, modification class)`.
- Random seeds, output tables, and figure semantics.

Do not silently relax thresholds after seeing a sparse PPI network. Treat an empty KG result as missing evidence, not evidence of no biological interaction.

Recommended workflow:

```text
AnnData audit
  -> stable feature identity
  -> frozen differential analysis
  -> admitted feature/protein table
  -> strict one-hop Protein-to-Protein interact_with query
  -> client-side evidence review and canonical-pair aggregation
  -> attributed NetworkX graph
  -> interactive HTML plus PNG/PDF/SVG and source tables
```

## 2. Environment and Required Input

Install the CLI and the differential-analysis dependency:

```bash
python -m pip install BiOmics-BRICK statsmodels
brick --help
brick qr --help
brick ga data2nxg --help
brick pl interact --help
```

Require an AnnData `.h5ad` file with:

- observations as samples and variables as proteins, peptides, or PTM-site features;
- a sample-group column such as `adata.obs["State"]`;
- an explicit protein identifier column such as `adata.var["Proteins"]` or `adata.var["UniProt"]`;
- a gene symbol or other KG-resolvable protein name;
- a modification class when performing PTM-aware analysis.

Never use `var_names_make_unique()` as the biological aggregation policy. Preserve duplicated peptide/PTM measurements and create stable feature IDs instead.

## 3. Audit Data and Preserve Feature Identity

Adapt only the declared column names and scale decision in this example:

```python
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

INPUT = Path("proteome.h5ad")
OUT = Path("proteome_results")
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

GROUP_COLUMN = "State"
CONTROL = "Control"
TREATMENT = "Infected"
PROTEIN_COLUMN = "Proteins"       # UniProt accession or other stable protein ID
MOD_COLUMN = "Modifications"      # use "Unmodified" if no PTM classes exist

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

adata = ad.read_h5ad(INPUT)
X = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
X = X.astype(float, copy=False)

assert X.shape == adata.shape
assert np.isfinite(X).all(), "Define missing/non-finite handling before analysis"
assert GROUP_COLUMN in adata.obs
assert set(adata.obs[GROUP_COLUMN]) == {CONTROL, TREATMENT}
assert PROTEIN_COLUMN in adata.var

var = adata.var.copy()
var["var_index"] = np.arange(adata.n_vars)
var["original_var_name"] = [str(x) for x in adata.var_names]
var["gene_symbol"] = var.get("gene_symbol", var["original_var_name"]).astype(str)
var["UniProt"] = var[PROTEIN_COLUMN].astype(str)
if MOD_COLUMN in var:
    var["Modifications"] = var[MOD_COLUMN].astype(str)
else:
    var["Modifications"] = "Unmodified"
var["feature_id"] = [
    f"{u}__{i:05d}" for i, u in zip(var["var_index"], var["UniProt"])
]

audit = {
    "input_sha256": sha256_file(INPUT),
    "n_samples": int(adata.n_obs),
    "n_features": int(adata.n_vars),
    "n_zero": int((X == 0).sum()),
    "n_negative": int((X < 0).sum()),
    "groups": adata.obs[GROUP_COLUMN].value_counts().to_dict(),
    "n_unique_features": int(var["feature_id"].nunique()),
    "n_unique_proteins": int(var["UniProt"].nunique()),
    "modification_counts": var["Modifications"].value_counts().to_dict(),
}
pd.Series(audit).to_json(OUT / "data_audit.json", indent=2)
var.to_csv(TABLES / "feature_identity_mapping.csv", index=False)
print(audit)
```

Inspect the audit before transforming data. If `X` is raw non-negative intensity, define a finite zero floor and use `log2`. If `X` is already log2 transformed, copy it without transforming again. Record the decision.

```python
ANALYSIS_SCALE = "raw_intensity"  # choose raw_intensity or already_log2

if ANALYSIS_SCALE == "raw_intensity":
    assert (X >= 0).all()
    positive = X[X > 0]
    assert positive.size > 0
    zero_floor = float(positive.min() / 2)
    L = np.log2(np.where(X == 0, zero_floor, X))
else:
    zero_floor = None
    L = X.copy()
```

Do not interpret a zero floor as biological imputation. If zero means missing or below detection in the assay, use and document an assay-appropriate missingness model instead.

## 4. Run Differential Analysis

Use a two-sided Welch t-test on the declared analysis scale, report log2 fold change and uncertainty, and control the false-discovery rate across all tested features. Adjust the design for paired, longitudinal, batch-confounded, or repeated-measures data; do not use this two-group example for those designs.

```python
FDR_Q = 0.05
ABS_LOG2FC = 1.0
ALLOW_EXPLORATORY_FALLBACK = True
NODE_MODE = "ptm_state"           # choose protein or ptm_state
CANDIDATE_CAP = 50
KG_RELATION = "interact_with"
KG_SOURCES = ["BioGRID_hs", "IntAct", "CORUM", "SIGNOR", "StringDB"]
LAYOUT_SEED = 42

analysis_contract = {
    "contrast": {
        "group_column": GROUP_COLUMN,
        "control": CONTROL,
        "treatment": TREATMENT,
    },
    "analysis_scale": ANALYSIS_SCALE,
    "zero_floor": zero_floor,
    "test": "two-sided Welch t-test",
    "multiple_testing": "Benjamini-Hochberg FDR across all features",
    "primary_threshold": {"BH_FDR": FDR_Q, "abs_log2FC": ABS_LOG2FC},
    "allow_exploratory_fallback": ALLOW_EXPLORATORY_FALLBACK,
    "node_mode": NODE_MODE,
    "candidate_cap": CANDIDATE_CAP,
    "kg": {
        "relation": KG_RELATION,
        "sources_requested": KG_SOURCES,
        "hops": 1,
        "directed": False,
        "strict_induced_subgraph": True,
    },
    "layout_seed": LAYOUT_SEED,
}
with open(OUT / "analysis_contract.json", "w", encoding="utf-8") as handle:
    json.dump(analysis_contract, handle, indent=2)

states = adata.obs[GROUP_COLUMN].astype(str).to_numpy()
idx_c = np.flatnonzero(states == CONTROL)
idx_t = np.flatnonzero(states == TREATMENT)
assert len(idx_c) >= 2 and len(idx_t) >= 2

rows = []
for j in range(adata.n_vars):
    control = L[idx_c, j]
    treatment = L[idx_t, j]
    t_stat, p_value = stats.ttest_ind(
        treatment, control, equal_var=False, nan_policy="raise"
    )
    if not np.isfinite(p_value):
        p_value = 1.0
    if not np.isfinite(t_stat):
        t_stat = 0.0
    log2fc = float(treatment.mean() - control.mean())
    se = float(np.sqrt(treatment.var(ddof=1) / len(treatment)
                       + control.var(ddof=1) / len(control)))
    rows.append({
        "var_index": j,
        "feature_id": var.iloc[j]["feature_id"],
        "gene_symbol": var.iloc[j]["gene_symbol"],
        "UniProt": var.iloc[j]["UniProt"],
        "Modifications": var.iloc[j]["Modifications"],
        "mean_log2_control": float(control.mean()),
        "mean_log2_treatment": float(treatment.mean()),
        "log2FC": log2fc,
        "ci95_lo": log2fc - 1.96 * se,
        "ci95_hi": log2fc + 1.96 * se,
        "welch_t": float(t_stat),
        "p_value": float(p_value),
    })

diff = pd.DataFrame(rows)
diff["BH_FDR"] = multipletests(diff["p_value"], method="fdr_bh")[1]
diff["abs_log2FC"] = diff["log2FC"].abs()
diff = diff.sort_values(["p_value", "abs_log2FC"], ascending=[True, False])
diff.to_csv(TABLES / "differential_results_all_features.csv", index=False)
```

Freeze the primary and exploratory rules before reading the ranked table. The values below are an explicit example, not a universal proteomics standard:

```python
primary = diff[
    (diff["BH_FDR"] < FDR_Q) &
    (diff["abs_log2FC"] >= ABS_LOG2FC)
].copy()

if not primary.empty:
    admitted = primary
    admission_label = "primary: BH_FDR<0.05 and |log2FC|>=1"
elif ALLOW_EXPLORATORY_FALLBACK:
    admitted = diff[
        (diff["p_value"] < 0.05) &
        (diff["abs_log2FC"] >= ABS_LOG2FC)
    ].copy()
    admission_label = "exploratory: raw Welch p<0.05 and |log2FC|>=1"
else:
    admitted = diff.iloc[0:0].copy()
    admission_label = "no admitted candidates"

admitted["admission_label"] = admission_label
admitted.to_csv(TABLES / "admitted_features.csv", index=False)
print(admission_label, "| features:", len(admitted))
```

Do not describe an uncorrected fallback as statistically significant. Do not recommend a blanket threshold such as `p < 0.3` without an explicit, pre-registered exploratory rationale.

## 5. Freeze Protein or PTM-State Candidates

Choose one node contract:

- `protein`: create one node per UniProt accession.
- `ptm_state`: create one node per `(UniProt, Modifications)` pair and preserve multiple measured states of the same protein.

Select one representative feature per node deterministically for node color and correlation. Keep all feature-level results in the differential table.

```python
if admitted.empty:
    raise RuntimeError(
        "No candidates passed the frozen admission rule; stop before KG retrieval"
    )

keys = ["UniProt"] if NODE_MODE == "protein" else ["UniProt", "Modifications"]
candidates = (
    admitted.sort_values(
        ["p_value", "abs_log2FC", "var_index"],
        ascending=[True, False, True],
    )
    .drop_duplicates(keys)
    .head(CANDIDATE_CAP)
    .copy()
)

mod_short = {
    "Phospho (STY)": "PHOS",
    "Deamidation (NQ)": "DEAM",
    "Glycosylated": "GLYCO",
    "Unmodified": "UNMOD",
}
if NODE_MODE == "protein":
    candidates["node_id"] = candidates["UniProt"]
else:
    candidates["node_id"] = [
        f"{u}::{mod_short.get(m, 'MOD')}"
        for u, m in zip(candidates["UniProt"], candidates["Modifications"])
    ]

candidates["query_name"] = candidates["gene_symbol"].astype(str)
candidates.to_csv(TABLES / "candidate_nodes.csv", index=False)

query_names = sorted(set(candidates["query_name"]) - {"", "nan"})
if len(query_names) < 2:
    raise RuntimeError("At least two resolved protein names are required for a PPI graph")
with open(OUT / "admitted_protein_names.txt", "w", encoding="utf-8") as handle:
    handle.write("\n".join(query_names) + "\n")
```

Validate the query namespace with a small dry run. In some BiOmics-KG snapshots, gene symbols resolve as `Protein` names while UniProt accessions do not. Record the successful namespace and retain the mapping table.

## 6. Query an Induced PPI Subgraph

Use `interact_with` for the current schema-valid PPI query. Do not use `activate`: the current `brick qr` schema does not recognize it as a query relation.

Pass the same admitted list as both source and target. Use strict target mode and disable fallback so the query cannot silently drop the target or relation constraint.

```bash
brick qr -i proteome_results/admitted_protein_names.txt \
  --target-file proteome_results/admitted_protein_names.txt \
  --st Protein --tt Protein \
  -r interact_with \
  --info-source BioGRID_hs,IntAct,CORUM,SIGNOR,StringDB \
  --hops 1 --target-mode strict --no-fallback \
  -n 5000 --dry-run-plan --show-schema

brick qr -i proteome_results/admitted_protein_names.txt \
  --target-file proteome_results/admitted_protein_names.txt \
  --st Protein --tt Protein \
  -r interact_with \
  --info-source BioGRID_hs,IntAct,CORUM,SIGNOR,StringDB \
  --hops 1 --target-mode strict --no-fallback \
  -n 5000 -o proteome_results/tables/kg_ppi_evidence_raw.csv
```

Treat `--info-source` as a requested server-side filter, then verify the returned `path.1.info_source` values client-side. Preserve raw rows even when the endpoint ignores or broadens the requested source filter.

For more than approximately 50-100 candidates, pre-register a cap or batch source entities. Keep the complete target set in every batch, merge raw evidence, and apply the same endpoint and evidence filters after merging.

## 7. Review and Deduplicate KG Evidence

Keep row-level evidence and create a separate edge table. Canonicalize undirected protein pairs, remove self-loops, exclude weak-trace relations from positive topology, and prioritize curated physical evidence over broader functional association.

```python
import ast

raw = pd.read_csv(TABLES / "kg_ppi_evidence_raw.csv")
required = {
    "path.0.name", "path.2.name", "path.1.relation", "path.1.info_source"
}
missing = required - set(raw.columns)
assert not missing, f"Missing brick columns: {sorted(missing)}"

allowed = set(query_names)
raw.to_csv(TABLES / "kg_ppi_evidence_unfiltered.csv", index=False)
ppi = raw[
    raw["path.0.name"].astype(str).isin(allowed) &
    raw["path.2.name"].astype(str).isin(allowed)
].copy()
ppi = ppi[ppi["path.0.name"].astype(str) != ppi["path.2.name"].astype(str)]
ppi = ppi[ppi["path.1.relation"].astype(str) == "interact_with"]

def parse_sources(value):
    if pd.isna(value):
        return []
    text = str(value).strip()
    try:
        parsed = ast.literal_eval(text)
        if isinstance(parsed, (list, tuple, set)):
            return [str(x) for x in parsed]
    except (SyntaxError, ValueError):
        pass
    return [text]

CURATED = {"BioGRID", "BIOGRID", "BioGRID_hs", "IntAct", "INTACT", "CORUM", "SIGNOR"}
FUNCTIONAL = {"StringDB", "STRINGDB", "StingDB", "StingDB_hs", "StingDB_mm"}

def evidence_tier(value):
    sources = set(parse_sources(value))
    if sources & CURATED:
        return "curated_physical"
    if sources & FUNCTIONAL:
        return "functional_association"
    return "unclassified_excluded"

ppi["edge_tier"] = ppi["path.1.info_source"].map(evidence_tier)
ppi["pair_key"] = ppi.apply(
    lambda r: "||".join(sorted([str(r["path.0.name"]), str(r["path.2.name"])])),
    axis=1,
)
ppi.to_csv(TABLES / "kg_ppi_evidence_admitted.csv", index=False)

edge_rows = []
for pair_key, group in ppi.groupby("pair_key", sort=True):
    tiers = set(group["edge_tier"])
    if "curated_physical" in tiers:
        tier = "curated_physical"
    elif "functional_association" in tiers:
        tier = "functional_association"
    else:
        continue
    source_name, target_name = pair_key.split("||", 1)
    sources = sorted({
        token
        for value in group["path.1.info_source"]
        for token in parse_sources(value)
    })
    edge_rows.append({
        "source_name": source_name,
        "target_name": target_name,
        "relation": "interact_with",
        "edge_tier": tier,
        "info_source": ";".join(sources),
        "evidence_rows": int(len(group)),
    })

edges = pd.DataFrame(edge_rows, columns=[
    "source_name", "target_name", "relation", "edge_tier",
    "info_source", "evidence_rows",
])
edges.to_csv(TABLES / "ppi_edges_canonical.csv", index=False)
```

Do not sort `relation_confidence` as though it were a universal probability. Preserve it with the raw evidence and interpret it according to the originating relation/source contract.

## 8. Build an Attributed Network

Use Python/NetworkX for an attributed Figure 7-style graph. The current `brick ga data2nxg -t df` command is suitable for a quick topology-only graph, but it does not preserve `path.1.*` evidence fields, node effects, PTM shapes, isolates, or correlation overlays.

For a quick topology check only, run:

```bash
brick ga data2nxg \
  -i proteome_results/tables/kg_ppi_evidence_admitted.csv \
  -o proteome_results/ppi_topology_only.graphml \
  -t df
brick ga info -i proteome_results/ppi_topology_only.graphml
```

Do not use this topology-only artifact as the attributed Figure 7c output.

```python
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import networkx as nx

G = nx.Graph()
shape_map = {
    "Phospho (STY)": "triangle",
    "Deamidation (NQ)": "square",
    "Glycosylated": "square",
    "Unmodified": "dot",
}
vmax = max(float(candidates["abs_log2FC"].max()), 1e-6)
norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

for row in candidates.itertuples(index=False):
    G.add_node(
        row.node_id,
        name=str(row.gene_symbol),
        label=str(row.gene_symbol),
        UniProt=str(row.UniProt),
        Modifications=str(row.Modifications),
        feature_id=str(row.feature_id),
        var_index=int(row.var_index),
        log2FC=float(row.log2FC),
        p_value=float(row.p_value),
        BH_FDR=float(row.BH_FDR),
        color=mcolors.to_hex(plt.cm.RdBu_r(norm(float(row.log2FC)))),
        shape=shape_map.get(str(row.Modifications), "dot"),
        size=18,
    )

nodes_by_query = (
    candidates.groupby("query_name")["node_id"].apply(list).to_dict()
)
design = np.column_stack([
    np.ones(adata.n_obs),
    (states == TREATMENT).astype(float),
])

def residualized_corr(node_a, node_b):
    ja = int(G.nodes[node_a]["var_index"])
    jb = int(G.nodes[node_b]["var_index"])
    xa, xb = L[:, ja], L[:, jb]
    ra = xa - design @ np.linalg.lstsq(design, xa, rcond=None)[0]
    rb = xb - design @ np.linalg.lstsq(design, xb, rcond=None)[0]
    value = float(np.corrcoef(ra, rb)[0, 1])
    return value if np.isfinite(value) else 0.0

for row in edges.itertuples(index=False):
    for node_a in nodes_by_query.get(row.source_name, []):
        for node_b in nodes_by_query.get(row.target_name, []):
            if node_a == node_b:
                continue
            corr = residualized_corr(node_a, node_b)
            G.add_edge(
                node_a,
                node_b,
                relation="interact_with",
                edge_tier=str(row.edge_tier),
                info_source=str(row.info_source),
                residualized_correlation=corr,
                color="#d62728" if corr >= 0 else "#1f77b4",
                weight=float(1.0 + 6.0 * abs(corr)),
                arrows="none",
            )

nx.write_graphml(G, OUT / "ppi_network_attributed.graphml")

node_source = pd.DataFrame.from_dict(dict(G.nodes(data=True)), orient="index")
node_source.index.name = "node_id"
node_source.reset_index().to_csv(TABLES / "ppi_network_nodes.csv", index=False)
nx.to_pandas_edgelist(G).to_csv(TABLES / "ppi_network_edges.csv", index=False)
```

Interpret the red/blue edge overlay as observational coordination after removing the group mean. Do not call it activation/inhibition or causal regulation. When using `ptm_state`, record that one protein-level KG edge may expand to multiple measured PTM-state node pairs.

## 9. Render Interactive and Publication Figures

Generate an interactive graph from the attributed GraphML:

```bash
brick pl interact \
  -i proteome_results/ppi_network_attributed.graphml \
  -o proteome_results/ppi_network.html \
  --width 1200 --height 900
```

Use custom Matplotlib drawing for publication figures because the current `brick pl static` renderer does not provide the complete triangle/square/PTM legend and evidence-tier styling needed for Figure 7c. Freeze the layout seed and export source coordinates.

```python
marker_map = {
    "Phospho (STY)": "^",
    "Deamidation (NQ)": "s",
    "Glycosylated": "s",
    "Unmodified": "o",
}
pos = nx.spring_layout(G, seed=LAYOUT_SEED)
pd.DataFrame([
    {"node_id": node, "x": float(xy[0]), "y": float(xy[1])}
    for node, xy in pos.items()
]).to_csv(TABLES / "network_coordinates.csv", index=False)

fig, ax = plt.subplots(figsize=(10, 8))
for a, b, data in G.edges(data=True):
    nx.draw_networkx_edges(
        G, pos, edgelist=[(a, b)], ax=ax,
        edge_color=data["color"], width=data["weight"],
        style="solid" if data["edge_tier"] == "curated_physical" else "dashed",
        alpha=0.8,
    )
for modification, marker in marker_map.items():
    nodes = [n for n, d in G.nodes(data=True) if d["Modifications"] == modification]
    if nodes:
        nx.draw_networkx_nodes(
            G, pos, nodelist=nodes, node_shape=marker,
            node_color=[G.nodes[n]["color"] for n in nodes],
            node_size=500, edgecolors="black", ax=ax,
        )
nx.draw_networkx_labels(G, pos, labels={n: G.nodes[n]["name"] for n in G},
                        font_size=8, ax=ax)
ax.set_axis_off()
for suffix in ("png", "pdf", "svg"):
    fig.savefig(FIGURES / f"ppi_network.{suffix}", dpi=300, bbox_inches="tight")
plt.close(fig)
```

Create the volcano plot from `differential_results_all_features.csv`. Use `x=log2FC` and `y=-log10(raw Welch p)` while encoding FDR-qualified and exploratory points differently. Label the axis as raw p-value when it is raw; do not visually imply that the y-axis is FDR.

## 10. Reproduce Legacy Figure 7 Deliberately

Separate exact legacy reproduction from the recommended workflow. The uploaded `Proreg_new.ipynb` contains history-dependent execution and inconsistent selection rules:

| Component | Legacy notebook behavior | Recommended behavior |
|---|---|---|
| Effect | `log2(mean(Infected)/mean(Control))` | mean difference on a declared log2 scale |
| Test | ordinary independent t-test on the available matrix | Welch test on the declared scale plus BH-FDR |
| Volcano coloring | `|logFC| > 0.2` and raw `p < 0.1` | pre-registered FDR/effect rule; label fallback explicitly |
| Other legacy selection | ratio `<0.6666` or `>1.5`; another description uses `|logFC|>0.5, p<0.3` | use one frozen admission rule |
| PPI relation | direct Neo4j `activate` query | current schema-valid `interact_with` query |
| PPI candidate handoff | the 149-protein `all_pro_set` was not consistently reused | pass the exact admitted list as source and target |
| Node shape | `P_` triangle, `G_` square, otherwise default | derive shape from the saved modification field |
| Node color | logFC with a blue-white-red scale | same, with saved scale limits |
| Edge overlay | Pearson correlation across all 10 samples | within-condition residualized correlation |

To reproduce the legacy figure exactly, pin the original input files, KG snapshot, Cypher schema, thresholds, execution order, color limits, and layout. Label the result as legacy reproduction. Do not claim numerical identity when using the recommended statistics or current `brick qr` schema.

Resolve the glycosylation/deamidation label from the data and assay contract. The legacy notebook names a `glycoproteome_result.txt` branch and prefixes it with `G_`, while a converted AnnData dataset may encode those features as `Deamidation (NQ)`. Preserve the raw annotation and describe any glycopeptide proxy interpretation explicitly; do not silently rename deamidation as direct glycosylation.

## 11. Interpret, Validate, and Report Results

Apply these evidence boundaries:

- Treat BioGRID/IntAct/CORUM/SIGNOR edges as curated or complex/signaling evidence according to the returned source record.
- Treat StringDB edges as functional association unless the originating record demonstrates a physical interaction.
- Treat node log2FC as measured differential abundance for the selected representative feature, not a KG property.
- Treat residualized correlation as observational coordination, not causality.
- Treat hub degree as graph- and query-dependent; do not call a hub a biological master regulator without independent evidence.
- Keep isolates in the candidate-node table even if the KG returns no admitted edge.
- Verify high-priority edges in their originating databases and validate PTM changes with independent cohorts or targeted assays.

Retain at least these outputs:

| Artifact | Purpose |
|---|---|
| `data_audit.json` | Input fingerprint, dimensions, groups, and measurement audit |
| `analysis_contract.json` | Frozen statistical, candidate, KG, and layout decisions |
| `feature_identity_mapping.csv` | Stable feature-to-protein/PTM mapping |
| `differential_results_all_features.csv` | Complete feature-level statistics |
| `admitted_features.csv` | Frozen threshold result and admission label |
| `candidate_nodes.csv` | Protein/PTM node contract and representative features |
| `kg_ppi_evidence_raw.csv` | Unmodified brick response |
| `kg_ppi_evidence_admitted.csv` | Endpoint-filtered row-level evidence |
| `ppi_edges_canonical.csv` | One canonical row per admitted protein pair |
| `ppi_network_nodes.csv` / `ppi_network_edges.csv` | Figure source tables |
| `network_coordinates.csv` | Deterministic layout coordinates |
| `ppi_network_attributed.graphml` | Attributed network |
| `ppi_network.html` | Interactive inspection |
| `ppi_network.png/.pdf/.svg` | Publication outputs |

Before reporting, verify:

1. Every code cell or script runs from a fresh environment without hidden notebook state.
2. Every candidate node maps to one declared representative feature and one KG query name.
3. Every network edge maps to at least one saved KG evidence row.
4. Both edge endpoints belong to the frozen candidate set.
5. No weak-trace or unclassified row enters positive topology by accident.
6. Counts in the report are computed from saved source tables.
7. Seeds, package versions, query parameters, fallbacks, endpoint failures, and file hashes are recorded.
