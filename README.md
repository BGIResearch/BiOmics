# BiOmics-KG-interpretation

An Agent Skill for querying and analyzing the **BiOmics Knowledge Graph (BiOmics-KG)**. Agents use it to invoke the BRICK toolchain for entity-relationship queries, result ranking/enrichment, graph conversion and merging, omics embedding, and visualization.

## Key Features

| Tool | Function | Typical Use |
|------|----------|-------------|
| `BRICK_qr` | Query | Cell markers, gene–disease associations, and other KG relationships |
| `BRICK_rk` | Rank / Enrich | Score and rank query results, or run enrichment analysis |
| `BRICK_ga` | Graph Adapter | Read/write graphs, format conversion, DataFrame/AnnData → graph, merge omics and KG graphs |
| `BRICK_emb` | Graph Embedding | Integrate single-cell data with KG results using PCA/UMAP or optional GRACE contrastive embeddings |
| `BRICK_pl` | Plot | Interactive/static network plots, radial enrichment plots, TF regulon plots |

Typical workflow: `BRICK_qr` → `BRICK_rk` → `BRICK_ga` / `BRICK_emb` → `BRICK_pl`

## Skill Architecture

```
BiOmics-KG-interpretation/
├── SKILL.md              # Main doc: tool reference, parameters, examples, troubleshooting
├── references/           # Domain workflow guides (loaded on demand)
│   ├── celltype_annotation_cli.md
│   ├── gene_enrichment.md
│   ├── drug_discovery_cli.md
│   └── ...
└── scripts/              # BRICK CLI implementation scripts
    ├── querygraph.py     # BRICK_qr
    ├── rankgraph.py      # BRICK_rk
    ├── graphadapter.py   # BRICK_ga
    ├── embedgraph.py     # BRICK_emb
    └── plotting/         # BRICK_pl
```

- **SKILL.md**: YAML frontmatter (`name` / `description`) plus the full BRICK tutorial; read first when the skill is triggered.
- **references/**: Scenario-specific workflows (cell annotation, DEG, drug repositioning, trajectory analysis, etc.) referenced by the agent as needed.
- **scripts/**: CLI entry points installed via `pip install BiOmics-BRICK`.

## Installation

### 1. Install BRICK CLI (required)

The CLI package lives in `BiOmics-KG-interpretation/scripts/` (where `pyproject.toml` is). **Run all install commands from that directory**, not the repo root.

#### Recommended: local editable install

**Windows (PowerShell)**

```powershell
cd BiOmics-KG-interpretation\scripts

# Use a dedicated venv (recommended; avoids Anaconda dependency conflicts)
python -m venv .venv
.\.venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
pip install -e .
```

For GRACE graph-contrastive embedding, install the optional runtime instead. The
legacy UMAP path remains available and does not require PyTorch or PyG:

```bash
pip install -e ".[grace]"
```

**Linux / macOS**

```bash
cd BiOmics-KG-interpretation/scripts

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
pip install -e .
```

`-e` installs in editable mode: code changes under `scripts/` take effect without reinstalling.

#### Verify

```bash
brick --help
brick agent --help
brick policy --help
brick qr --help
```

You should see subcommands including `qr`, `rk`, `ga`, `pl`, `emb`, `reason`, `plan`, `workflow`, `policy`, and `agent`.

#### Other install options

**Non-editable install**

```bash
pip install .
```

**Install with pipx** (isolated global CLI)

```bash
pipx install -e /path/to/BiOmics-KG-interpretation/scripts
```

**Run without installing**

```bash
cd BiOmics-KG-interpretation/scripts
python cli.py --help
```

#### If dependency install fails

On Anaconda or other shared environments, dependency conflicts (for example `llvmlite`) may block a full install. Try:

```bash
pip install -e . --no-deps
pip install pandas requests scipy networkx matplotlib pyvis pydot
pip install scanpy   # optional; needed for h5ad / celltype workflow
```

Core KG query, ranking, and reasoning commands work without `scanpy`; single-cell workflows (`brick workflow celltype`, `brick emb`) require it.

#### Requirements and notes

- Python >= 3.9
- Network access for `brick qr` (default API: `https://biomics.dcs.cloud/api`)
- Optional: [Graphviz](https://graphviz.org/download/) for `brick ga` SVG/PDF export
- If `brick` is not found after install, ensure the virtual environment is activated and open a new terminal if needed

### 2. Add to Your Agent

Copy or symlink the `BiOmics-KG-interpretation/` directory (including `SKILL.md` and subdirectories) into the skills path for your agent.

#### Cursor

| Scope | Path |
|-------|------|
| Personal (all projects) | `~/.cursor/skills/BiOmics-KG-interpretation/` |
| Project (shared via repo) | `.cursor/skills/BiOmics-KG-interpretation/` |

```bash
# Project-level example
mkdir -p .cursor/skills
cp -r BiOmics-KG-interpretation .cursor/skills/
```

#### Hermes

| Scope | Path |
|-------|------|
| Default | `~/.hermes/skills/BiOmics-KG-interpretation/` |
| Team-shared directory | Set `HERMES_OPTIONAL_SKILLS_DIR` to a directory containing this skill |

```bash
# Local directory
cp -r BiOmics-KG-interpretation ~/.hermes/skills/

# Or install from Hub / URL (single-file SKILL.md)
hermes skills install <hub-path-or-url>
```

In a session, search with `/skills search BiOmics` or invoke the skill by name.

#### Claude Code

| Scope | Path |
|-------|------|
| Personal | `~/.claude/skills/BiOmics-KG-interpretation/SKILL.md` |
| Project | `.claude/skills/BiOmics-KG-interpretation/SKILL.md` |

```bash
mkdir -p .claude/skills
cp -r BiOmics-KG-interpretation .claude/skills/
```

Claude Code picks up directory changes automatically — no restart needed.

#### Codex

| Scope | Path |
|-------|------|
| Project (recommended) | `.agents/skills/BiOmics-KG-interpretation/` |
| Global | `~/.agents/skills/BiOmics-KG-interpretation/` |

```bash
mkdir -p .agents/skills
cp -r BiOmics-KG-interpretation .agents/skills/
```

Codex scans `.agents/skills/` from the current working directory up to the repo root. Restart Codex if the skill does not appear. You can also use `$skill-installer` inside Codex to install published skill packages.

## GRACE Graph Embedding

GRACE is an optional `brick emb` method for analyses where the topology of the
integrated single-cell/KG graph should contribute to the representation. It
learns from two independently perturbed graph views using a shared GCN encoder
and a symmetric InfoNCE objective. The existing PCA/UMAP implementation remains
the default, so existing commands keep their previous behavior unless
`--method grace` is supplied.

Typical uses include cell-subtype refinement, prioritizing KG candidates near
omics cells, and comparing a topology-aware representation with the legacy
embedding. Embedding proximity is a prioritization signal, not proof of a
marker relationship, mechanism, or causal effect; always inspect the original
KG relations and provenance before making a biological claim.

### Inputs

| Argument | Required content |
|----------|------------------|
| `-a`, `--adata` | AnnData `.h5ad`; `var_names` must contain genes used by the query CSV |
| `-i`, `--input` | CSV produced by `brick qr` or `brick rk`, with `path.0.name`, `path.0.type`, `path.2.name`, and `path.2.type` |
| `-o`, `--output` | Output integrated graph in AnnData `.h5ad` format |

Rows whose Gene node is absent from `adata.var_names` are excluded. The command
stops with an explicit error if no query genes remain. Add `--cell-neighbor`
only when `adata.obsp['connectivities']` already contains a suitable cell-cell
neighbor graph.

### Install and run

Run the install command from `BiOmics-KG-interpretation/scripts/`:

```bash
pip install -e ".[grace]"

brick emb \
  -a adata_sub.h5ad \
  -i marker_hits.csv \
  --method grace \
  --cell-neighbor \
  --grace-epochs 300 \
  --grace-hidden-dim 32 \
  --seed 20241204 \
  -o adata_grace.h5ad
```

The command composes cell-gene expression edges, KG query edges, and optional
cell-cell edges. By default it uses a seeded projected node identity plus node
type and log-degree features, trains GRACE, computes UMAP from `X_grace`, and
calculates cosine distances between omics cells and KG entity nodes.

### Main parameters

| Parameter | Default | Purpose |
|-----------|--------:|---------|
| `--grace-epochs` | `300` | Number of contrastive training epochs |
| `--grace-hidden-dim` | `32` | Width of the learned node representation |
| `--grace-layers` | `2` | Number of GCN layers |
| `--grace-edge-drop` | `0.3` | Edge-removal probability in each augmented view |
| `--grace-feature-mask` | `0.3` | Feature-masking probability in each augmented view |
| `--grace-feature-mode` | `projected_identity` | Scalable default; use `identity` only for a small upstream-style comparison |
| `--grace-max-nodes` | `5000` | Safety cap for the all-node contrastive objective; `0` disables the cap |
| `--grace-device` | `auto` | Select CPU, CUDA, or Apple MPS automatically |
| `--distance-metric` | `cosine` | Metric for omics-cell to KG-entity distances under GRACE |
| `--seed` | `20241204` | Controls features, augmentation, initialization, and UMAP |

Run `brick emb --help` for the complete parameter list.

### Output fields

| AnnData field | Content |
|---------------|---------|
| `.X` | Learned GRACE node representations |
| `.obsm['X_grace']` | Explicit copy of the learned representation |
| `.obsm['X_umap']` | Two-dimensional UMAP computed from `X_grace` |
| `.obsp['integrated_adjacency']` | Original integrated graph adjacency matrix |
| `.uns['KG_Omics_distance']` | Cell-by-KG candidate distance table |
| `.uns['grace']` | Loss curve, parameters, graph size, device, seed, and upstream source |
| `.uns['brick_embedding']` | Compact method and distance metadata |

Example inspection:

```python
import scanpy as sc

graph = sc.read_h5ad("adata_grace.h5ad")
print(graph.uns["brick_embedding"])
print(graph.uns["grace"]["losses"])
print(graph.uns["KG_Omics_distance"].head())
```

### Validation and scaling notes

- Inspect the loss curve and confirm all embeddings and distances are finite.
- Repeat important analyses with a second seed and compare against the legacy
  `--method umap` result.
- The all-node InfoNCE objective uses O(N²) similarity matrices. Keep the graph
  biologically focused before raising `--grace-max-nodes`.
- `--grace-feature-mode identity` creates a dense N×N matrix and is guarded by
  `--grace-identity-max-nodes`; use the default `projected_identity` for routine
  analysis.
- Review nearest candidates together with marker expression, KG relation type,
  source database, and supporting evidence. Do not interpret distance alone as
  biological validation.

For more detail, see
[`references/grace_embedding.md`](BiOmics-KG-interpretation/references/grace_embedding.md).

## Quick Example

```bash
# Query T cell marker genes
BRICK_qr "T cell" --st Cell --tt Gene -r marker_of -o query.csv

# Rank and visualize
BRICK_rk -i query.csv --rank -o ranked.csv
BRICK_pl radial -i ranked.csv --top-k 10 -o result.png
```

For full parameters and troubleshooting, see [`BiOmics-KG-interpretation/SKILL.md`](BiOmics-KG-interpretation/SKILL.md).
