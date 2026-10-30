# GRACE Integrated Graph Embedding

Use this reference when `brick emb --method grace` is part of an analysis.

## Purpose

GRACE learns node representations from two independently perturbed views of the integrated omics/KG graph. Each view removes a fraction of edges and masks a fraction of node features; a shared GCN encoder and InfoNCE objective pull the two representations of the same node together while separating other nodes.

This is an analysis representation, not an evidence-strength upgrade. A short embedding distance can prioritize a cell-to-KG candidate for follow-up, but the biological claim must still be checked against the original KG relation, provenance, and expression evidence.

The implementation is adapted from [`embedgraph/GRACE.py`](https://github.com/caolei2-BGI/BRICK/blob/main/embedgraph/GRACE.py). It keeps the upstream two-view GCN/InfoNCE design and adds a scalable node-feature mode, CLI parameterization, optional dependencies, deterministic seeds, AnnData metadata, and graph-size safeguards.

## Installation

```bash
pip install "BiOmics-BRICK[grace]"
```

The optional extra installs PyTorch and PyTorch Geometric. The GRACE loss and augmentations are implemented inside BRICK, avoiding PyGCL's additional DGL runtime. The legacy `--method umap` path does not require this extra beyond the base BRICK dependencies.

## Recommended Analysis

Build a focused KG query first. For cell-type refinement, keep candidate cell types biologically plausible and use curated marker sources when possible.

```bash
brick qr "pancreatic D cell, pancreatic PP cell" \
  --st Cell --tt Gene -r produce \
  --info-source CellMarker,PanglaoDB \
  -n 500 -o sub_adata_markers.csv

brick emb \
  -a adata_sub.h5ad \
  -i sub_adata_markers.csv \
  --method grace \
  --cell-neighbor \
  --grace-epochs 300 \
  --grace-hidden-dim 32 \
  --seed 20241204 \
  -o adata_grace.h5ad
```

Use `--cell-neighbor` only when `adata.obsp['connectivities']` has already been computed and is appropriate for the subset. Without it, cells remain connected to the graph through expressed query genes.

## Main Parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `--grace-epochs` | `300` | Contrastive training epochs |
| `--grace-hidden-dim` | `32` | Learned node embedding width |
| `--grace-layers` | `2` | GCN layers |
| `--grace-learning-rate` | `0.0001` | Adam learning rate |
| `--grace-edge-drop` | `0.3` | Edge-removal probability in each view |
| `--grace-feature-mask` | `0.3` | Feature-masking probability in each view |
| `--grace-temperature` | `0.2` | InfoNCE temperature |
| `--grace-feature-mode` | `projected_identity` | Scalable node identity projection plus type/degree features |
| `--grace-feature-dim` | `64` | Width of the projected identity component |
| `--grace-max-nodes` | `5000` | Safety cap for the all-node contrastive objective; `0` disables it |
| `--grace-device` | `auto` | Select `cuda`, `mps`, or `cpu` |
| `--distance-metric` | `cosine` | Distance used for cell-to-KG comparison under GRACE |
| `--seed` | `20241204` | Feature, augmentation, initialization, and UMAP seed |

`projected_identity` is the analysis default because the upstream dense identity feature matrix uses quadratic memory. Use `--grace-feature-mode identity` only for a small reproducibility comparison; `--grace-identity-max-nodes` guards accidental large allocations.

The local-to-local InfoNCE loss also forms all-node similarity matrices. Subset the biological question before raising `--grace-max-nodes`, and monitor host/GPU memory.

## Output Contract

The output `.h5ad` contains:

- `.X`: learned GRACE node representations.
- `.obsm['X_grace']`: the same explicit GRACE representation.
- `.obsm['X_umap']`: UMAP coordinates computed from `X_grace`.
- `.obsp['integrated_adjacency']`: the original integrated graph adjacency.
- `.uns['KG_Omics_distance']`: cell-by-KG distance table, cosine by default.
- `.uns['grace']`: hyperparameters, graph size, device, seed, per-epoch losses, final loss, and upstream source.
- `.uns['brick_embedding']`: compact method and distance metadata shared with other embedding modes.

## Analysis Checks

1. Confirm the retained query genes exist in `adata.var_names` and the query relation/source scope matches the question.
2. Plot `.uns['grace']['losses']`; investigate non-finite values, a flat curve, or an unstable late rise.
3. Repeat at least one important analysis with a second seed. Treat candidates that disappear as embedding-sensitive.
4. Compare GRACE with the legacy `--method umap` result and with expression-only annotations. Report major disagreements.
5. Inspect nearest KG candidates together with their marker genes and raw KG provenance. Distance alone is not causal or curated evidence.
6. For subtype assignment, examine the margin between the nearest and second-nearest KG candidates instead of using only `idxmin`.
