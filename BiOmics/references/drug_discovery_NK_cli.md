# NK-cell drug discovery with BiOmics-BRICK

Use this reference for an NK-cell adaptation that starts from a complete
single-cell atlas, derives a two-level NK gene union, expands it with
BiOmics-KG, learns a joint GRACE representation, and ranks Chemical entities
against one Disease anchor. This is a procedural contract only: it contains no
expected genes, chemicals, ranks, scores, counts, or historical result paths.

For this task, the package-parity rules in this reference override the generic
query-expansion and source-routing defaults in `SKILL.md`. Do not perform query
expansion, relation narrowing, or `info_source` filtering unless this reference
explicitly requires it. Method parity is assessed from frozen inputs, ordered
intermediate artifacts, and implementation provenance; no intermediate count
or candidate identity is an expected result.

## Contents

1. Integrity boundary
2. Runtime preflight
3. Input and DEG audit
4. One-hop KG retrieval
5. Enrichment and graph construction
6. Primary GRACE fit and Chemical ranking
7. Ranking freeze and reference query
8. Dynamic diagnostic
9. Executed notebook and QA
10. Required artifacts and hard stops

## 1. Integrity boundary

Before reading outcomes, freeze these task parameters in the run provenance:

- input AnnData path and SHA-256;
- target cell label and cell-type column;
- state column and returned state labels;
- top-gene count and exact DEG calls;
- DEG conversion implementation, ordered group labels, and ordered gene union;
- KG source/target types, hop count, physical batch size, logical query limit,
  relation policy, API endpoint, and implementation hash;
- enrichment implementation, target identity, background method, statistic,
  p-value field, and threshold;
- graph implementation, native identity policy, ordered nodes, and ordered
  edges;
- Disease anchor label and type;
- GRACE device, seed policy, epochs, learning rate, and model defaults;
- ranking metric, tie-break rule, and validation cutoffs.

Treat the task prompt and listed data as the only task inputs. Do not inspect
historical notebooks, result folders, candidate lists, rankings, or reference
tables. Do not put expected candidates or values into the prompt, skill, code,
or tests. Preserve the input AnnData as read-only and write all outputs into a
new task directory.

The independent Disease-to-Chemical reference query is evaluation-only. Do
not issue it before the complete primary Chemical ranking has been saved and
hashed.

## 2. Runtime preflight

Record Python and package versions, the `brick` and `brick-gpu` help output,
hostname, CPU/memory, CUDA visibility, GPU model, and random seeds. Verify that
the CLI and the package GRACE model resolve in the same prepared environment.

Run a one-epoch CUDA smoke fit in a separate `smoke/` directory before the
primary fit. The smoke output must never replace or tune the primary graph,
embedding, ranking, seed, or hyperparameters. If DGL/GraphBolt is incompatible
with the installed PyTorch version, use only the task-local compatibility
mechanism supplied by the runtime. Do not install a fake global `dgl` package
or replace the prepared environment.

## 3. Input and DEG audit

Require an AnnData object with a measured expression matrix, a cell-type
column, and a clinical-state column. Inspect shape, matrix representation,
layers/raw, existing UMAP, cell/state counts, donors, batches, and missingness.
Do not normalize, log-transform, scale, or recompute neighbors unless the
task explicitly requests it or the audit proves the required representation
is absent. Never apply `log1p` twice.

For the package-parity NK route, use the exact two-level DEG sequence below.
Do not add a method, reference group, layer, filter, or preprocessing step to
either `rank_genes_groups` call.

```python
from pathlib import Path

import scanpy as sc

adata = sc.read_h5ad("data/input_atlas.h5ad")

celltype_col = "cell_type_fine"
state_col = "Ventilated"
target_celltype = "NK"
top_n = 100

required = {celltype_col, state_col}
missing = sorted(required.difference(adata.obs.columns))
if missing:
    raise ValueError(f"Missing required obs columns: {missing}")

celltype_values = adata.obs[celltype_col].astype(str)
if not celltype_values.eq(target_celltype).any():
    raise ValueError(f"No cells with {celltype_col} == {target_celltype!r}")

selected_gene_sources = {}

sc.tl.rank_genes_groups(adata, celltype_col)
atlas_tables = rank_genes_groups_tables(adata)
atlas_groups = list(atlas_tables)
if target_celltype not in atlas_tables:
    raise ValueError("Target cell type is absent from atlas DEG groups")
for gene in atlas_tables[target_celltype]["names"].astype(str).iloc[:top_n]:
    selected_gene_sources.setdefault(gene, []).append("atlas_NK")

adata_nk = adata[celltype_values.eq(target_celltype)].copy()
if adata_nk.obs[state_col].astype(str).nunique() < 2:
    raise ValueError("NK subset has fewer than two clinical states")

sc.tl.rank_genes_groups(adata_nk, state_col)
state_tables = rank_genes_groups_tables(adata_nk)
state_groups = list(state_tables)
for group, table in state_tables.items():
    for gene in table["names"].astype(str).iloc[:top_n]:
        selected_gene_sources.setdefault(gene, []).append(f"nk_state_{group}")

ordered_genes = sorted(
    gene for gene in selected_gene_sources if gene not in {"", "nan", "None"}
)
if not ordered_genes:
    raise ValueError("Primary gene union is empty")

adata_nk.uns["query_gene_set"] = ordered_genes
adata_nk.write_h5ad("nk_primary.h5ad")
Path("primary_gene_union.txt").write_text("\n".join(ordered_genes) + "\n")
```

`rank_genes_groups_tables` must first use the installed
`BRICK.pp.rank_genes_groups2df` when it is importable. A compatibility
implementation is allowed only to convert the structured arrays into the same
per-group rows. It must compare its `names` values and order against
`adata.uns["rank_genes_groups"]["names"]`, record the implementation path and
hash, and stop on a mismatch. It may not filter, re-sort, clean, map, or curate
genes. In particular, do not add p-value, fold-change, direction, marker,
transcription-factor, symbol-conversion, or expected-union-size filters.

Never pass a Python `set` to a query, AnnData slice, graph constructor, or
serializer. Freeze and hash the sorted `ordered_genes` list before KG access.

Save the complete atlas-level and NK-state DEG tables, every selected top-100
row with its source group, and a gene-union provenance table. Record group
order, duplicate removal, unmapped identifiers, and the biological replicate
limitation; do not treat cells as independent donors in interpretation.

## 4. One-hop KG retrieval

Query the frozen gene union from source type `Gene` to this ordered target
union unless the task overrides it:

```text
Disease, Phenotype, Process, Cell_Component, Pathway,
Function, Tissue, Cell, Chemical, Stage
```

Use one hop and one logical multi-target query. When the package must split the
ordered genes into physical request batches, every batch must carry the same
ordered target-type union. Preserve raw triplets, stable IDs, relation
direction, `original_relation`, `info_source`, condition, confidence, batch
inputs, attempts, responses, and failures. Derive per-target CSVs from the
successful merged logical response; do not issue ten scientifically distinct
queries and do not merge partial results.

```python
target_types = [
    "Disease", "Phenotype", "Process", "Cell_Component", "Pathway",
    "Function", "Tissue", "Cell", "Chemical", "Stage",
]

kg_triplets = query_neighbor_fail_closed(
    ordered_genes,
    source_entity_type="Gene",
    target_entity_type=target_types,
    return_type="triplet",
    multi_hop=1,
    limit_number=None,
)
```

The wrapper must use the installed package/API query primitive, serialize the
target union exactly once per physical batch, and fail the whole stage if any
batch fails after its predeclared retries. It must not silently continue,
substitute target types, expand names, filter relations or sources, or add a
per-target row cap. Record the physical batch size, number of batches, response
rows, de-duplication key, per-target counts, matched-source coverage, explicit
limit (`null` when absent), and whether a hidden server cap can be detected.
Never describe a live response as historically complete when the KG snapshot
or an undocumented server cap is not available for verification.

## 5. Enrichment and graph construction

Run the installed package enrichment primitive once on the merged raw table.
For package parity, retain rows with raw enrichment p-value below `0.01` only
after saving the complete enrichment table. Accept `pvalue` or
`path.2.enrich_pvalue` only as a recorded schema alias; do not change the
statistic. Preserve FDR when available.

```bash
brick rk -i kg_multitype_triplets.csv --enrich -o enrichment_all.csv
```

Record the exact background method, target identity key, package origin, and
package/source hash. The primary route requires an exact KG background. A
query-subset/degraded background may be saved as a diagnostic but must hard-stop
the primary workflow. Do not change the threshold after seeing graph size, GPU
memory, ranking, or validation overlap.

Build the graph with the installed `BRICK.pp.combine_graph` primitive (or a
compatibility wrapper that calls the same native graph adapter) from
`adata_nk[:, ordered_genes]` and the retained enrichment rows. Do not add edges,
inject entities, select a connected component, or redesign node/edge semantics.

The native implementation may key nodes by labels. Before calling it, audit
every target identity by entity type and stable ID. If distinct non-Gene target
identities share a label, or a non-Gene target label collides with an omics or
source node, replace only the internal graph key with the deterministic form
`<display label> [<stable ID>]`. Keep the original display label in the ordered
identity table and in every ranking/report output. This graph-key qualification
must preserve every input row, entity, edge, type, and feature; it may not
select, filter, collapse, inject, or reconnect anything. Hard stop when a
collision cannot be resolved uniquely, when a Gene label maps to multiple
stable IDs, or when a ranked Chemical or Disease anchor does not map to exactly
one stable ID. Record the complete qualification mapping, input-row count, and
that no filtering occurred. Hash the ordered genes, filtered enrichment rows,
ordered node list, ordered edge list, and canonical topology.

Before fitting, require exactly one graph node whose type is `Disease` and
label is the frozen Disease anchor. If it is absent or duplicated, stop. Never
inject a standalone anchor, add synthetic edges, or borrow the later reference
query to repair the primary graph.

## 6. Primary GRACE fit and Chemical ranking

Fit the package `GRACE_model` on the audited native graph, on CUDA, with exactly
the three seeds `20241204`, `20241205`, and `20241206`, each for exactly 5,000
epochs. Every seed must receive the identical graph, ordered node list, feature
matrix, and hyperparameters, with an independent model initialization. Use the
installed model defaults unless this reference freezes an override; record all
resolved defaults and RNG state. Save the complete loss trace, ordered node
identity mapping, embedding matrix, runtime, GPU memory, and model import
origin.

When the local CLI exposes the audited drug route, prefer:

```bash
brick-gpu drug grace-rank \
  -g primary_joint_graph.graphml \
  --identity primary_joint_graph_nodes.csv \
  --disease COVID-19 \
  --epochs 5000 \
  --device cuda \
  -o primary_grace
```

Otherwise call the package model directly through a small checked Python
launcher. Do not replace GRACE with PCA or UMAP. Use PCA/UMAP only to visualize
the saved GRACE embeddings.

Within each seed, rank every `Chemical` node by descending raw cosine similarity
to the unique Disease anchor. The primary consensus is ascending mean raw rank
across all three seeds, with stable KG ID and then label as deterministic tie
breaks. Do not select a seed or use validation to aggregate seeds. Compute
Euclidean distance as a diagnostic only. Stable KG ID is the primary identity;
labels are display fields. Save all candidates, non-finite checks, ties, graph
degree, and deterministic tie-break fields. Do not filter elements,
metabolites, radicals, groups, or other Chemical entities with a manual name
blacklist.

## 7. Ranking freeze and reference query

Write the full Chemical ranking and SHA-256 before reference access:

```bash
brick drug freeze -i primary_chemical_ranking.csv -o ranking_frozen.sha256
```

Only then query the Disease anchor to `Chemical` with
`pharmaco_effect_on`, preserving the complete raw response:

```bash
brick qr COVID-19 \
  --st Disease \
  --tt Chemical \
  -r pharmaco_effect_on \
  -n 500 \
  -o covid_chemical_reference.csv
```

Validate the unchanged ranking at predeclared cutoffs. Match by stable ID
first and report exact/normalized labels as separate diagnostics. Report total
reference IDs, eligible reference IDs in the ranked universe, recovered
eligible IDs, and references absent from the universe. A valid zero-row result
must still be a parseable CSV with expected headers; a blank file is failure.

## 8. Dynamic diagnostic

Run one separate diagnostic GRACE model with learning rate `0.001`. Execute
100 one-epoch updates and capture the Chemical-to-anchor cosine similarity
after each actual update. Save the 100 losses, every per-epoch Chemical score,
and per-Chemical Spearman trend.

Label snapshots by actual update count (`1..100`). Do not relabel one update as
ten epochs, and do not use this diagnostic model to replace, tune, or reorder
the primary 5,000-epoch ranking.

## 9. Executed notebook and QA

Create a new reader-facing notebook from the task prompt and current run only.
It must contain the executable workflow or execute one auditable driver that
can rebuild the workflow from the supplied H5AD. A driver call must store its
command, return code, stdout, stderr, input hashes, output hashes, and stage
status in the notebook. A notebook that only reads final tables, loads saved
embeddings, or displays pre-rendered images is incomplete. The notebook must be
executable after deleting task-generated outputs while preserving the supplied
input, installed skill, task-local compatibility layer, and driver.

Produce these method-output roles without referring to any earlier notebook or
cell numbering:

| Figure ID | Analysis-output role |
| --- | --- |
| F01 | Input data UMAP |
| F02 | GRACE loss |
| F03 | Joint embedding UMAP |
| F04 | Joint embedding alternate view |
| F05 | Enrichment match-count distribution |
| F06 | Target entity-type composition |
| F07 | Full Chemical-Disease similarity heatmap |
| F08 | Enriched Chemical-Disease similarity heatmap |
| F09 | Chemical similarity versus distance |
| F10 | Post-freeze reference-highlight view |
| F11 | Reference versus other Chemical distribution |
| F12 | Post-reference distance diagnostic |
| F13 | Dynamic trend diagnostic |
| F14 | Dynamic-section distribution |

Use relative paths and checked subprocess return codes. Save every figure as
both PNG and PDF, with labels, title, caption, and a data-derived denominator.
Label the dynamic diagnostic by its actual update count `1..100`.

Write `FIGURE_MANIFEST.csv` with at least these fields: `figure_id`, `role`,
`notebook_cell`, `png_path`, `pdf_path`, `source_data_paths`, `denominator`,
`status`, and `deviation_note`. A role passes only when both files exist, are
non-empty, can be decoded, and the executed notebook contains its displayed
output. Derive `FIGURE_QA.json` from the reopened notebook and manifest, not
from expected constants.

Execute once from a fresh kernel and save outputs in the same notebook. Reopen
it with `nbformat`; count a code cell as executed only when
`execution_count is not None`. Require all counts to be non-null, unique, and
monotonically increasing, with zero error outputs. Derive `NOTEBOOK_QA.json`
from the reopened notebook. Never hardcode pass, counts, or booleans.

## 10. Required artifacts and hard stops

At minimum retain:

- input/DEG/gene-union audits and `nk_primary.h5ad`;
- per-target raw KG tables and merged triplets;
- complete/filtered enrichment tables and background audit;
- native graph, node identity table, and graph audit;
- primary loss, embeddings, runtime audit, and complete Chemical ranking;
- ranking freeze, raw reference table, identity audit, and validation table;
- dynamic losses, scores, and trends;
- all analysis figures, `FIGURE_MANIFEST.csv`, `FIGURE_QA.json`, executed
  notebook, `NOTEBOOK_QA.json`, provenance, and report.

Stop with a concrete blocker instead of reporting success when required input
fields/groups are absent, any KG physical batch fails, enrichment is empty or
uses a degraded background, the Disease anchor is absent/non-unique, CUDA or
GRACE fails, the graph or any of the three 5,000-epoch fits cannot finish with
frozen parameters, ranking is not frozen before validation, or notebook QA
fails. Do not convert a hard stop into success by changing the gene set,
threshold, background, graph, anchor, epochs, model, seed set, or QA logic.
