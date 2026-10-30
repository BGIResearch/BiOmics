# Cell Development Trajectory Analysis with brick CLI

## Task Overview

Cell development trajectory inference is an important part of single-cell analysis, used to reveal cell differentiation paths and developmental relationships. This tutorial demonstrates how to use brick CLI tools to complete an end-to-end trajectory analysis workflow: computing PAGA trajectories from single-cell data, querying brick knowledge graph for cell differentiation relationships, merging data layer and knowledge layer results, and generating distinguishable visualizations.

**Note: This workflow is for reference only. Please adjust according to user needs and real-time conditions.**

## Prerequisites

- BiOmics-BRICK package installed: `pip install BiOmics-BRICK`
- scanpy installed: `pip install scanpy`
- AnnData format single-cell data prepared (preprocessing, clustering, and dimensionality reduction completed)

## Analysis Workflow

### Step 1: Compute PAGA Trajectory

Use scanpy to compute PAGA (Partition-based Graph Abstraction) trajectory. To obtain a concise tree-like structure, use `connectivities_tree` instead of dense `connectivities`.

> **Note**: The following code is for demonstration only. Please ensure `sc.pp.neighbors` has been completed before PAGA computation.

```python
import scanpy as sc
import pandas as pd
import networkx as nx
import os
import subprocess

DATA_FILE = "./processed_data.h5ad"
OUTPUT_DIR = "./trajectory_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

adata = sc.read(DATA_FILE)
sc.tl.paga(adata, groups='cell_type')

# Save AnnData containing PAGA results
adata.write(os.path.join(OUTPUT_DIR, "adata_paga.h5ad"))
```

### Step 2: Convert PAGA Results to Graph File

Use `brick ga data2nxg` to extract PAGA results from AnnData into a NetworkX graph file for subsequent merging and visualization.

```bash
brick ga data2nxg \
    -i ./trajectory_output/adata_paga.h5ad \
    -o ./trajectory_output/paga.graphml \
    -t paga
```

**Parameter Description:**

| Parameter | Description |
|-----------|-------------|
| `data2nxg` | Convert data to NetworkX graph |
| `-i` | Input AnnData file (PAGA computed) |
| `-o` | Output graph file (`.graphml` format) |
| `-t paga` | Input type is PAGA |

### Step 3: Query brick KG Differentiation Relationships

Extract cell type nodes from the PAGA graph and use `brick qr` to query developmental relationships in the knowledge graph. Query **both** `develop_into` and `is_a` (server-side OR pattern `is_a|develop_into`; note the KG relation is spelled `develop_into` — the legacy typo `develope_into` is also accepted by the schema) and allow up to **10 hops** (`--hops 10`) so direct differentiation edges and ontological parent-child links can be connected through intermediate cell types.

This is the CLI equivalent of the legacy notebook workflow `brick.qr.query_shortest_path(..., relation="is_a|develop_into", multi_hop=(1,10))` followed by `filter_results`: pass the same cell type list as both source (`-i`) and target (`--target-file`), and use **`--shortest`** so each **unordered** cell pair gets exactly one server-side `shortestPath()` (one POST per source entity to `/api/v1/shortest_path`, targets passed as a list; already-queried unordered pairs are skipped, giving C(n,2) pairs total).

> **Why `--shortest` uses the server-side endpoint**: the generic `/api/v1/relation` endpoint returns all matching paths subject to an unordered server-side `LIMIT`. With a small `-n` the truncation hides entire source-target pairs; with a large `-n` the server enumerates all multi-hop paths and hits the 60s timeout (HTTP 504). `--shortest` instead delegates the per-pair `shortestPath()` computation to the database via `/api/v1/shortest_path`, which is fast and exact for deep hops.
>
> The `/api/v1/shortest_path` endpoint relies on two KG data conventions: path edges must carry a `relation` property, and the API response includes a `path.N.reversed` flag per edge so the CLI can restore the true edge direction (path column order is traversal order, which is arbitrary for undirected queries; `is_a` points child → parent regardless of traversal order).

```python
# Extract cell type list
g_paga = nx.read_graphml("./trajectory_output/paga.graphml")
cell_types = list(g_paga.nodes())

ct_file = os.path.join(OUTPUT_DIR, "cell_types.txt")
with open(ct_file, "w") as f:
    f.write("\n".join(cell_types))
```

```bash
brick qr -i ./trajectory_output/cell_types.txt \
    --target-file ./trajectory_output/cell_types.txt \
    --st Cell --tt Cell \
    -r "is_a|develop_into" \
    --hops 10 \
    --shortest \
    -o ./trajectory_output/kg_develope_into.csv
```

**Parameter Description:**

| Parameter | Description |
|-----------|-------------|
| `-i` | Input cell type list file (source entities) |
| `--target-file` | Restrict query endpoints to the same cell type list (pairwise source-target queries) |
| `--st Cell` | Source entity type is Cell |
| `--tt Cell` | Target entity type is Cell |
| `-r "is_a|develop_into"` | Query differentiation and ontological hierarchy relationships (`\|` = server-side OR) |
| `--hops 10` | Allow up to 10-hop paths between cell types (safe with `--shortest`: the DB computes one shortest path per pair) |
| `--shortest` | Server-side per-pair `shortestPath()` via `/api/v1/shortest_path`; unordered pairs only (C(n,2) paths). Exact legacy `query_shortest_path` semantics; no LIMIT truncation or 60s-timeout issues |
| `-o` | Output CSV file |

### Step 4: Merge PAGA Graph with KG Query Results

Use `brick ga complete_result` to merge PAGA original graph with KG query results into an annotated graph. This command automatically detects input format, filters node types, handles relationship direction, and highlights original edges not forming paths in query results.

```bash
brick ga complete_result \
    -g ./trajectory_output/paga.graphml \
    -q ./trajectory_output/kg_develope_into.csv \
    -o ./trajectory_output/merged.graphml \
    --node-attrs '{"color":"gray","shape":"hexagon","size":12}' \
    --edge-attrs '{"color":"green","weight":8,"arrows":"to"}' \
    --node-type-filter Cell \
    --reverse-relations is_a
```

**Parameter Description:**

| Parameter | Description |
|-----------|-------------|
| `-g` | Original graph file (PAGA) |
| `-q` | Query result file (CSV or graphml) |
| `-o` | Output merged graph file |
| `--node-attrs` | Default attributes for new KG nodes (JSON) |
| `--edge-attrs` | Default attributes for KG edges (JSON) |
| `--node-type-filter` | Filter query results by node type |
| `--reverse-relations` | List of relations with reversed direction (e.g., `is_a`) |

### Step 5: Post-processing —Distinguish Original Cells from KG Supplementary Cells

In the merged graph, original PAGA cells have scanpy auto-assigned colors, while KG supplementary cells are gray. To ensure correct color distinction under tree layout, execute the following post-processing: retain PAGA original colors, set KG nodes uniformly to gray, and distinguish by shape (PAGA as hexagon, KG as ellipse).

```python
G = nx.read_graphml("./trajectory_output/merged.graphml")
paga_nodes = set(g_paga.nodes())

for n in G.nodes:
    if n in paga_nodes:
        G.nodes[n]["shape"] = "hexagon"
        G.nodes[n]["size"] = 25
        # Retain original PAGA color, do not overwrite
    else:
        G.nodes[n]["color"] = "gray"
        G.nodes[n]["shape"] = "ellipse"
        G.nodes[n]["size"] = 18

# Change highlighted disconnected edges from red to light gray, avoiding being too prominent
for u, v, d in G.edges(data=True):
    if d.get("color") == "red":
        d["color"] = "#AAAAAA"
        d["weight"] = 1.5
        d["arrows"] = "none"

nx.write_graphml(G, "./trajectory_output/trajectory_processed.graphml")
```

### Step 6: Generate Visualization

Use `brick pl` to generate static and interactive HTML. Tree layout is suitable for displaying hierarchical differentiation relationships; spring layout is suitable for displaying overall network structure.

```bash
# Tree layout (recommended for trajectory inference)
brick pl static \
    -i ./trajectory_output/trajectory_processed.graphml \
    -o ./trajectory_output/trajectory_tree.png \
    --layout tree \
    --root "hematopoietic stem cell"

# Spring layout (display overall network)
brick pl static \
    -i ./trajectory_output/trajectory_processed.graphml \
    -o ./trajectory_output/trajectory_spring.png

# Interactive HTML (supports zoom, drag, click to view)
brick pl interact \
    -i ./trajectory_output/trajectory_processed.graphml \
    -o ./trajectory_output/trajectory.html \
    --width 1400 --height 1000
```

**Parameter Description:**

| Parameter | Description |
|-----------|-------------|
| `--layout tree` | Tree hierarchical layout (requires `--root`) |
| `--root` | Specify root node (e.g., `hematopoietic stem cell`) |
| `--width` / `--height` | Interactive canvas size (pixels) |

## Result Interpretation

- **paga.graphml**: PAGA tree structure (`connectivities_tree`), containing only original cell type nodes
- **kg_develope_into.csv**: brick KG query results, containing `develop_into`, `is_a`, and multi-hop developmental paths (up to 10 hops; with `--shortest`, exactly one shortest path per unordered cell pair)
- **merged.graphml**: Original graph after PAGA and KG merge
- **trajectory_processed.graphml**: Final graph after post-processing (colors/shapes distinguished)
- **trajectory_tree.png**: Tree layout static image, clear hierarchy
- **trajectory_spring.png**: Force-directed layout static image, overall structure at a glance
- **trajectory.html**: Interactive network, supports zooming to view node labels and edge attributes

### Legend

| Visual Element | Meaning |
|---------------|---------|
| Colored hexagon | PAGA cells from original data (retaining scanpy colors) |
| Gray ellipse | Cells supplemented by brick KG (`develope_into` / `is_a` targets from multi-hop query) |
| Green arrow edges | KG-verified developmental relationships (`develope_into` / `is_a`, `arrows='to'`) |
| Gray edges | PAGA data-inferred trajectory edges (`connectivities_tree`) |

## Notes

1. **Query both `develop_into` and `is_a` with up to 10 hops**: `develop_into` captures direct differentiation evidence; `is_a` adds ontological hierarchy. Use `-r "is_a|develop_into"` (server-side OR). `--hops 10` connects PAGA cell types through intermediate nodes; pair it with `--shortest` (server-side per-pair `shortestPath()` via `/api/v1/shortest_path`, the exact CLI equivalent of legacy `query_shortest_path` + `filter_results`, unordered C(n,2) pairs). Do not use the default `/api/v1/relation` route with a large `-n` for this workflow: it returns unordered LIMIT-truncated rows, and deep multi-hop enumeration may hit the server-side 60s timeout (HTTP 504). Use `--node-type-filter Cell` in Step 4 and Step 5 post-processing to keep the merged graph readable when parent nodes appear.
2. **Use `connectivities_tree`**: PAGA's `connectivities_tree` retains only tree backbone edges (usually 5—0), far fewer than dense `connectivities` (up to dozens), resulting in cleaner visualization.
3. **Root node selection**: Tree layout requires specifying a node with in-degree 0 as root (e.g., `hematopoietic stem cell`). If no in-degree 0 node exists in the graph, it will automatically attempt to generate a minimum spanning tree.
4. **Post-processing necessity**: `complete_result` retains original PAGA colors. Without Step 5 post-processing, some PAGA nodes may not display labels in spring layout due to color not equaling `default_node_color`.
5. **Node scale**: If there are too many cell types (>20), tree layout may experience label crowding. At this point, appropriately filter query results, or prioritize using interactive HTML for viewing.
