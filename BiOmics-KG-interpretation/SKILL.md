---
name: BiOmics-KG-interpretation
description: "Use when building modular BiOmics-BRICK analyses over BiOmics-KG: compose query, ranking, enrichment, graph adaptation, embedding, plotting, and reasoning bricks for bioinformatics tasks such as cell type annotation, trajectory interpretation, proteomics/PTM-aware PPI analysis, enrichment, drug repositioning, evidence tracing, and knowledge-graph hypothesis generation. Trigger when the user asks to use brick/BiOmics, query biomedical knowledge graph relationships, analyze proteome or post-translational-modification data with KG evidence, or design customizable brick workflows rather than a fixed pipeline."
---

# BiOmics-BRICK

BiOmics-BRICK is a set of composable command-line bricks for building biological interpretability models: omics inputs and BiOmics-KG relationships become evidence units, ranked hypotheses, enriched entities, merged graphs, plots, and traceable explanations.

Do not treat brick as a single fixed pipeline. Start from the biological question, choose the smallest useful bricks, then compose them into a transparent workflow the user can inspect and modify.

## Installation
**Linux / macOS**

```bash
cd BiOmics-KG-interpretation/scripts

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
pip install -e .
```

## On-Demand BRICK Extension

When the existing bricks and their compositions cannot perform a required analysis, the agent may extend the local BRICK source code instead of stopping at the missing capability or substituting an unsupported result. This includes adding task-required statistical methods, adapters, evidence transformations, validators, or other reusable operations.

Use the smallest safe extension workflow:

1. First check whether existing bricks, parameters, or compositions can solve the task. Modify the source only when the required capability is genuinely absent or insufficient.
2. Inspect neighboring implementations and preserve the repository's existing interfaces and conventions. Keep changes scoped to the task, avoid unrelated refactoring, and do not overwrite uncommitted user work.
3. Add or update the relevant module under `scripts/`; when introducing a user-facing command, register it in the CLI and document its inputs, outputs, parameters, assumptions, and limitations.
4. Install or refresh the local package in editable mode when needed so the extended `brick` command is available in the current environment:

   ```bash
   cd BiOmics-KG-interpretation/scripts
   python -m pip install -e .
   ```

5. Verify the extension before using its results: run focused tests, check the relevant `brick ... --help` entry, and execute a small end-to-end smoke test on representative data. For statistical extensions, check method assumptions and report effect sizes, uncertainty, multiple-testing correction, and fallback behavior where applicable.
6. Preserve the changed source, tests, installation record, and example invocation as reproducibility artifacts. In the final report, list the modified files and validation performed. Do not claim that an extension works if installation or validation failed.

Source extension is local and task-scoped by default. Do not commit, push, publish a release, or modify shared deployments unless the user explicitly requests it. If a durable repository change is not wanted, keep the extension in the current task workspace and clearly identify its source path.

## Operating Model

1. Identify the user's input: gene list, marker files, AnnData `.h5ad`, existing graph, brick query CSV, ranked CSV, or paper-derived entity set.
2. Identify the desired output: candidate table, enrichment table, augmented graph, figure, evidence summary, comparison against known references, or a reusable command recipe.
3. Choose source and target entity types before querying. Common types include `Gene`, `Cell`, `Disease`, `Chemical`, `Pathway`, `Function`, `Process`, `Phenotype`, `Tissue`, `Protein`, and `Cell_Component`. Match each query string to the base-ontology naming rules in **Query Entity Naming Contract** below.
4. Use standard English entity names in `brick qr` inputs. Do not pass informal abbreviations, acronyms, or colloquial aliases directly to `brick qr`. When the user question contains abbreviations, entity classes, disease subtypes, repeat motifs, vector wording, variant wording, or TF/driver wording, run `brick reason expand-query` first, then keep only base-ontology standard English names for the actual KG query.
5. Use relation constraints when they make the biological question sharper. Common KG relations include `produce`, `regulated`, `interact_with`, `participate_in`, `causally_related_to`, `disease_relation`, `is_a`, `develope_into`, and `pharmaco_effect_on`. Concept aliases such as `marker_of`, `regulate`, and `associated_with` are resolved before querying.
6. Select the smallest useful edge-level `info_source` scope before querying when the question clearly points to marker, drug, pathway, variant, disease, regulation, protein-interaction, ontology, or anatomy evidence. Nodes do not store database/source provenance; provenance is carried by relationship edges, returned as `path.1.info_source` and used during offline reasoning (it is not accepted as a `brick qr` query filter). For node-type primary databases, relation meanings, source-routing templates, and quality rules, see `references/kg_data_sources.md`.
7. Use one-hop KG querying by default. Enable optional multi-hop querying with `brick qr --hops 2` or `--hops 3` only when the question needs indirect mechanism, bridging entities, or path-style evidence; never exceed 3 hops. Exception: trajectory/shortest-path use cases (e.g., cell differentiation paths with `is_a|develop_into`) may go up to `--hops 10` when paired with `--shortest` (server-side per-pair `shortestPath()` via `/api/v1/shortest_path`, no LIMIT truncation or 60s-timeout risk).
8. Normalize raw query rows into evidence units when the result will feed a reasoning layer or model-like report.
9. Preserve intermediate files. brick should leave query, expansion, rank, graph, plot, plan, and reasoning artifacts that can be reviewed independently.
10. If the remote BiOmics API is unreachable, still build the workflow folder, commands, graph operations, plots, and offline reasoning from existing artifacts. Do not hide API failures inside a workflow.

## Entity Input Integrity

Preserve KG entity literals exactly. For SNP/HGVS values or any entity name containing spaces, commas, colons, or parentheses, prefer `brick qr -i entities.txt` or `--target-file targets.txt` with exactly one entity per non-empty line. These file modes remove only surrounding whitespace and preserve internal punctuation. The CLI rejoins adjacent source/target words for compatibility, but quote a single positional entity in the shell, for example `brick qr "Diabetes Mellitus, Type 2" ...`. Do not use legacy comma-separated positional input when a literal entity name itself contains a comma. When using `--dry-run-plan`, execute the emitted `command_sh`; it is shell-quoted to preserve each entity as one argument.

## Brick Map

| Brick | Role | Use When |
|---|---|---|
| `brick qr` | Query BiOmics-KG relationships | Need raw entity-relation-entity evidence; use `--hops 2/3` only for optional multi-hop paths, up to `--hops 10` with `--shortest` for trajectory/shortest-path evidence |
| `brick rk --rank` | Rank candidate targets | Need cell type, disease, chemical, pathway, or other target prioritization |
| `brick rk --enrich` | Enrichment analysis | Need over-represented biological entities from a gene/entity set |
| `brick ga data2nxg` | Convert tabular, AnnData, or PAGA data to graph | Need a graph artifact from omics or query data |
| `brick ga complete_result` | Merge an observed graph with KG query results | Need to annotate PAGA, trajectory, or other inferred graphs with KG evidence |
| `brick ga convert/read/info` | Inspect and convert graphs | Need interoperability across `.graphml`, `.gml`, `.pkl`, `.gpickle`, `.json`, `.svg`, or `.pdf` |
| `brick emb` | Embed omics plus KG associations with PCA/UMAP or GRACE | Need an integrated AnnData embedding for downstream exploration or graph-contrastive analysis |
| `brick pl` | Visualize graph, enrichment, or TF regulons | Need static, interactive, radial, or TF network plots |
| `brick reason expand-query` | Expand query terms before KG calls | Need aliases, acronyms, class members, motif variants, or target/relation hints |
| `brick reason filter-co-mentioned` | Drop unsupported co-mentioned query rows | Need to review `path.*.original_relation` and remove weak co-occurrence-only hits before rank/explain |
| `brick reason literature` | Inspect publications behind co-mentioned edges between two entities | Need the PMID → title/journal/date + source sentence for a specific entity pair |
| `brick reason` | Summarize evidence and hypotheses offline | Need explanation, evidence tracing, path search, or prediction/reference comparison |
| `brick plan` | Describe, inspect, and materialize block plans | Need machine-readable model blueprints, artifact contracts, or reusable block recipes |
| `brick workflow` | Generate editable workflow scaffolds | Need a fast starting recipe for `celltype`, `trajectory`, `enrichment`, or `drug` |
| `brick agent` | Check, dry-run, execute, and report deterministic brick agent plans | Need BioMaster-like plan execution while preserving KG source policy, bounded multi-hop, and evidence validation |
| `brick drug` | Stable-ID-first Chemical scope, GRACE consensus, freeze, and validation | Need leakage-resistant drug prioritization helpers |

Prefer direct brick composition for custom scientific questions. Use `brick workflow` as a scaffold generator, then edit the generated `commands.ps1` or `commands.sh`.

## Interpretability Model Layers

Use this layered contract when the goal is a biological interpretability foundation model:

| Layer | Artifact | brick Support |
|---|---|---|
| Query expansion layer | candidate query terms, aliases, class members, type/relation hints | `brick reason expand-query` |
| Query layer | raw KG edge rows | `brick qr` |
| Co-mentioned review layer | filtered query CSV with unsupported co-mentioned rows removed | `brick reason filter-co-mentioned` |
| Evidence layer | canonical `evidence_id, source, relation, target, evidence_text, provenance_file` table | `brick reason normalize` |
| Grounding layer | canonical entity names and IDs with grounding confidence | `brick reason ground` |
| Scoring layer | support score, polarity, negative/conflict evidence | `brick reason score`, `brick reason conflicts` |
| Hypothesis layer | ranked or enriched candidate table | `brick rk`, `brick reason hypotheses` |
| Graph layer | observed graph plus KG-completed graph | `brick ga`, `brick pl`, `brick reason paths` |
| Mechanism layer | scored explanatory paths between entities | `brick reason path-score` |
| Store layer | reusable local evidence database | `brick reason store` |
| Evaluation layer | reference-set metrics and evidence coverage | `brick reason eval` |
| Plan layer | JSON block DAG with inputs, outputs, commands, and API dependencies | `brick plan template/validate/materialize` |
| Agent orchestration layer | checked agent plan, run manifest, validation events, execution report | `brick agent check-plan/check/run/report` |
| Explanation layer | markdown or CSV report with sources, relations, metrics, snippets, and references | `brick reason evidence/explain/compare` |

The evidence layer is the core model substrate. A downstream model, agent, or report should reason over normalized evidence units rather than raw path-column CSVs whenever possible.

## Query Entity Naming Contract

BiOmics-KG entity matching is name-based at query time. Use the **standard English preferred label** that the KG base ontologies would recognize. The core base libraries are:

| `brick qr` entity type | Base library | Query naming rule |
|---|---|---|
| `Gene` | NCBI Gene | Official gene symbol for the intended species (e.g., human `TP53`, `HNF4A`). Do not use informal aliases, old symbols, or descriptive prose. |
| `Cell` | CL | Cell Ontology preferred label in English (e.g., `hepatocyte`, `CD4-positive, alpha-beta T cell`). Avoid cluster nicknames or author-specific labels unless mapped to CL first. |
| `Disease` | DOID | Disease Ontology standard English disease name (e.g., `colorectal cancer`, `Duchenne muscular dystrophy`). Do not query with acronyms such as `CRC`, `GBM`, or `AMD`. |
| `Phenotype` | HPO | Human Phenotype Ontology term label in English (e.g., `Seizure`, `Intellectual disability`). |
| `Tissue` | UBERON | UBERON anatomy term in English (e.g., `liver`, `cerebral cortex`). Prefer anatomy terms over vague region names. |
| `Pathway` | KEGG | Official KEGG pathway name in English (e.g., `Apoptosis`, `Cell cycle`). |
| `Chemical` | MeSH | MeSH preferred English descriptor/name (e.g., `Sildenafil`, `Bevacizumab`). Avoid class phrases such as `PDE5 inhibitors` in direct `brick qr` input. |
| `Function`, `Process`, `Cell_Component` | GO | Official Gene Ontology term name in English (e.g., `DNA binding`, `cell cycle`, `mitochondrion`). |
| `Stage` | Stage (`HSAPDV`, `MMUSDV`) | Standard developmental stage label in English for the intended species (e.g., human `embryonic stage`, mouse `Theiler stage 26`). |

### Required rules

1. **Standard English only.** Query strings must use the canonical English label from the relevant base ontology, not Chinese, mixed-language, or free-text prose.
2. **No informal abbreviations or acronyms in `brick qr`.** Examples to avoid in direct queries: `CRC`, `GBM`, `AMD`, `PDE5`, `DMD`, `CMS`, `ZIKV`. Resolve these with `brick reason expand-query`, ontology lookup, or `brick reason ground`, then query with the expanded standard names.
3. **Official gene symbols are allowed.** Symbols such as `TP53`, `IL2`, and `HNF4A` are not treated as informal abbreviations because they are the NCBI Gene standard for query matching. Keep species-appropriate symbols; do not mix human and mouse symbols in one query without explicit intent.
4. **Match entity type to ontology.** Do not query a CL cell label with `--st Gene`, or a GO process with `--tt Disease`. Wrong typing is a common cause of empty KG results.
5. **Prefer specific terms over overly broad parents.** Use `hepatocyte` rather than `cell`, and `colorectal cancer` rather than `cancer`, when the biological question allows it.
6. **Preserve official punctuation in valid names.** Hyphens and numerals that belong to the standard name should be kept (e.g., `COVID-19`, `IL-2`, `CD4-positive, alpha-beta T cell`).
7. **One clean entity token per query item.** Input files may use one entity per line or comma-separated lists. Remove prose, parenthetical notes, IDs, and ranking annotations before querying (e.g., use `hepatocyte`, not `hepatocyte (cluster 3 marker)`).
8. **Do not mix ontology IDs and labels unless validated.** Prefer standard names for `brick qr`. If an ID is used, confirm it resolves in the KG for the chosen entity type.
9. **Resolve omics-derived labels before querying.** Single-cell cluster names, DE table rows, GWAS locus descriptions, and paper shorthand should be mapped to the base ontology naming rules above before `brick qr`.

### Recommended workflow when naming is uncertain

```bash
# 1) Expand abbreviations / class phrases from the user question
brick reason expand-query -q "..." --entity "..." --entity-type Disease --format txt -o expanded_terms.txt

# 2) Query only with standard English names
brick qr -i expanded_terms.txt --st Disease --tt Gene -r disease_relation -n 100 -o query.csv

# 3) Optional post-query canonicalization for reports
brick reason ground -i evidence_units.csv --lexicon entity_lexicon.csv -o evidence_grounded.csv
```

For base-library priorities, cross-database mapping, and post-query source review, see `references/kg_data_sources.md`.

## Query Expansion Contract

Run query expansion before KG querying when a single literal term is likely to miss evidence, especially when the user input violates the **Query Entity Naming Contract** (abbreviations, class phrases, aliases, or non-standard labels).

Use `brick reason expand-query` for:

- Abbreviations and aliases: `AMD`, `DMD`, `GBM`, `CRC`, `CMS`, `PDE5`, `UCHL1`, `RFC1`, disease acronyms, gene aliases, and parenthetical forms.
- Entity classes: drug classes, inhibitor classes, TF families, cell-state labels, disease subtypes, and subtype labels.
- Motif and variant questions: exact motif plus comparator motif and disease context, while keeping motif identity separate.
- Vector or infection questions: pathogen plus likely organism/species target terms.
- BioASQ-style list questions where recall depends on class-member expansion.

Example:

```bash
brick reason expand-query -q "Please list the short half life Phosphodiesterase 5 inhibitors." --entity "PDE5 inhibitors" --entity-type Chemical --format json -o query_expansion.json
brick qr -i expanded_drugs.txt --st Chemical --tt Protein -r pharmaco_effect_on -n 50 -o kg_pde5_drugs.csv
```

Expansion terms are hypotheses for querying, not final answers. After querying, validate them with brick evidence, relation type, snippet relevance, and negative/conflict checks. For class/list questions, include likely negative or contrastive members, such as long-half-life drugs, so the reasoning layer can exclude them explicitly.

## Composition Patterns

### Query -> Rank -> Reason

Use for candidate prioritization.

```bash
brick qr -i genes.txt --st Gene --tt Cell -r produce -n 500 -o query_cell.csv
brick reason filter-co-mentioned -i query_cell.csv -o query_cell_filtered.csv --rejected-output query_cell_rejected.csv
brick reason normalize -i query_cell_filtered.csv -o evidence_units.csv
brick reason ground -i evidence_units.csv --lexicon entity_lexicon.csv -o evidence_grounded.csv
brick reason score -i evidence_grounded.csv --aggregate target -o target_support.csv
brick rk -i query_cell_filtered.csv --rank --metrics "match_count,info_source_count" --prune-parent -o ranked_cell.csv
brick reason hypotheses -i ranked_cell.csv --top-n 20 -o hypotheses.csv
brick reason explain -e evidence_grounded.csv --hypotheses ranked_cell.csv -o explanation.md
```

### Query -> Enrich -> Plot

Use for gene-set interpretation across multiple entity types.

```bash
brick qr -i genes.txt --st Gene --tt Pathway -n 1000 -o query_pathway.csv
brick rk -i query_pathway.csv --enrich -o enriched_pathway.csv
brick pl radial -i enriched_pathway.csv --top-k 15 -o enriched_pathway.png
```

### Omics -> Markers -> KG Candidates

Use when starting from `.h5ad`. Extract markers or DE genes first, then query and rank.

```bash
brick workflow celltype --adata data.h5ad --groupby leiden -o celltype_recipe
```

Then inspect and edit `celltype_recipe/commands.ps1` before running API-dependent steps.

### Omics + KG -> GRACE Embedding

Use GRACE when graph topology should shape the joint cell/gene/KG representation, especially for subtype refinement or candidate proximity analysis. Install the optional runtime, start from a biologically focused query CSV, and preserve the seed and training metadata stored in the output AnnData.

```bash
pip install -e ".[grace]"
brick emb -a adata_sub.h5ad -i marker_hits.csv --method grace --grace-epochs 300 --cell-neighbor -o adata_grace.h5ad
```

Use the default `projected_identity` node features for ordinary analysis. Use `--grace-feature-mode identity` only to reproduce the upstream dense one-hot construction on small graphs. Inspect `uns['grace']`, the loss curve, graph size, and `uns['KG_Omics_distance']` before interpreting proximity as biological support. See `references/grace_embedding.md` for parameters, outputs, scalability limits, and validation guidance.

### Observed Graph -> KG Completion -> Plot -> Paths

Use for trajectory, PAGA, lineage, or network interpretation.

```bash
brick ga complete_result -g paga.graphml -q kg_develope_into.csv -o merged.graphml --node-type-filter Cell
brick pl static -i merged.graphml --layout tree --root "hematopoietic stem cell" -o trajectory_tree.png
brick reason paths -g merged.graphml --source "stem cell" --target "T cell" -o paths.csv
```

### Prediction -> Reference Comparison

Use when evaluating candidate drugs, diseases, cell types, or pathways against known references.

```bash
brick reason compare -p predicted_drugs.csv -r known_drugs.csv -o validation.csv
```

## Task Recipes

Use these as starting shapes, not mandatory pipelines.

| Task | Minimal Brick Chain | Reference |
|---|---|---|
| Cell type annotation | marker genes -> `brick qr Gene->Cell produce` -> `brick rk --rank` -> `brick reason evidence` | `references/celltype_annotation_cli.md`, `references/celltype_refinement.md` |
| Trajectory interpretation | PAGA/trajectory graph -> `brick qr Cell->Cell develope_into` -> `brick ga complete_result` -> `brick pl static/interact` -> `brick reason paths` | `references/trajectory_analysis_cli.md` |
| Enrichment | gene set -> multi-target `brick qr` -> merge CSVs -> `brick rk --enrich` -> `brick pl radial` | `references/gene_enrichment.md` |
| Drug repositioning | disease/state DE genes -> independent Chemical ranking or fixed-seed GPU GRACE -> ranking freeze -> stable-ID post hoc validation | `references/drug_discovery_cli.md` |
| DEG-driven analysis | AnnData groups -> DE genes -> downstream brick query/rank/enrich | `references/DEG_analysis_cli.md` |
| GRN interpretation | TF/target sets -> KG validation/enrichment -> TF plot | `references/grn_analysis_cli.md` |
| GWAS/SNP reasoning | variants or candidate genes -> phenotype/disease/function KG evidence | `references/gwas_causal_snps_cli.md`, `references/gwas_phenotype_prediction_cli.md` |
| Proteome analysis | protein/gene set -> KG relationships and enrichment | `references/proteome_analysis_cli.md` |
| KG source selection | biological question -> entity type -> relation type -> post-query review of edge `info_source` | `references/kg_data_sources.md` |

Read the relevant reference only when that task is requested. For general brick use, rely on the brick map and command `--help` output.


## Controlled KG query policy

Before querying BiOmics-KG, choose a KG policy rather than blindly querying the whole graph. The optimized brick package exposes this as a first-class command:

```bash
brick policy route \
  -q "Which drugs target up-regulated macrophage genes?" \
  --source-type Gene --target-type Chemical \
  --kg-policy auto --command
```

The policy output records `use_kg`, inferred `intent`, relation whitelist, edge-level `info_source` whitelist, `max_hops`, required edge fields, and the `co_mentioned` handling rule. Use `--kg-policy off` when KG evidence is not needed, `direct` for explicit user-supplied filters only, `source-scoped` to require source filtering, and `multihop` for limited 2/3-hop mechanism/path questions.

`brick qr` accepts the same controls directly:

```bash
brick qr -i markers.txt --st Gene --tt Cell -r produce \
  --kg-policy auto \
  --info-source CellMarker,PanglaoDB,CellTaxonomy \
  --hops 1 -o marker_cell_query.csv
```

Broad queries with no target type, relation, explicit target, or `info_source` filter are blocked unless `--allow-broad-query` is supplied. Use `--dry-run-plan` to inspect the resolved policy and filters without calling the KG API.

To prevent policy drift between code and the reference JSON, export the built-in catalog when updating defaults and compare it with `references/kg_source_policy.json`:

```bash
brick policy export-default -o kg_source_policy.json
```

## KG Source Selection

Do not query the whole graph by default when the question maps cleanly to a specific evidence family. Choose a focused edge-level `relation_info_source` scope first, then broaden only if recall is too low. BiOmics-KG nodes do not record the original database/source; only edges carry provenance through `path.1.info_source`.

Use `brick qr --info-source` for command-line edge provenance filtering:

```bash
brick qr -i genes.txt --st Gene --tt Cell -r produce --info-source CellMarker,PanglaoDB -n 100 -o marker_hits.csv
brick qr "TP53" --st Gene --tt Disease -r causally_related_to --info-source ClinVar,OMIM,PHARMGKB -n 50 -o tp53_disease.csv
```

High-signal routes:

- Cell marker questions: `CellMarker`, `CellMarker_mm`, `PanglaoDB`, `CellTaxonomy`, then ontology sources such as `CL` and `CLO`.
- Drug, treatment, adverse event, or pharmacology questions: `DrugBank`, `DrugCentral`, `IUPHAR`, `PHARMGKB`, `STITCH`, `REPODB`, `SIDER`, `ClinicalTrials` (legacy raw source may appear as `ClinicalTrails`).
- Disease, phenotype, GWAS, and variant questions: `ClinVar`, `ClinVar_gene`, `GWAS_Catalog` (legacy raw source may appear as `GWAS_cataog`), `OMIM`, `HPO_anno`, `HPO_g2p_anno`, `MONDO`, `DOID`, `ORPHANET`, `MEDGEN`.
- Pathway and function questions: `GOA`, `GOA_mm`, `GO`, `KEGG_anno`, `KEGG_anno_mm`, `Reactome_anno`, `Reactome`, `SMPDB`, `METACYC`.
- Protein interaction and complex questions: `StringDB`, `StingDB`, `StingDB_mm`, `BioGRID_hs`, `BioGRID_mm`, `IntAct`, `CORUM`, `SIGNOR`.
- Regulation, TF, and signaling questions: `TRRUST`, `SCENIC`, `SCENIC_mm`, `PhosphoSitePlus`, `Reactome`, `CellChat`.
- Tissue, anatomy, developmental-stage, and trajectory questions: `UBERON`, `BTO`, `EMAPA`, `MA`, `CL`, `CLO`, `EFO`.

For the node primary database priorities, NER coverage summary, relation quick reference, query templates, quality-control rules, external reference links, and full source catalog, read `references/kg_data_sources.md`.

## KG schema-aware query policy

BiOmics-KG queries must be schema-aware. Do not invent relation names.

Allowed node types are `Cell`, `Cell_Component`, `Chemical`, `Disease`, `Function`, `Gene`, `Mutation`, `Other`, `Pathway`, `Phenotype`, `Process`, `Protein`, `Species`, `Stage`, and `Tissue`.

Allowed KG edge types are `capable_of`, `causally_related_to`, `chemical_relation`, `chemicals_relation`, `compare`, `converse_from_to`, `develope_into`, `disease_relation`, `equals_to`, `functionally_related_to`, `has_part`, `interact_with`, `is_a`, `no_relation`, `participate_in`, `pharmaco_effect_on`, `processly_related_to`, `produce`, `regulated`, `sequentially_related_to`, `spatially_related_to`, `species_relation`, `temporally_related_to`, and `transport_from_to`.

Important schema rules:

- `no_relation` may appear in raw KG rows, but must not be used as a query relation.
- `associated_with`, `associated`, `association`, `related`, `related_to`, `linked_to`, `relevant_to`, and `involved_in` are not valid hard query relations.
- If the question implies broad association, omit the relation constraint instead of using `associated_with` or `related_to`.
- `mentioned` and `co_mentioned` are weak literature-trace evidence. They can be inspected during soft probing, but should not be treated as curated mechanistic evidence.

Safe relation aliases:

- `marker_of` and `markers_of` -> `produce`
- `regulate`, `regulates`, and `regulated_by` -> `regulated`
- `interact`, `interacts`, and `interacts_with` -> `interact_with`
- `drug_effect`, `drug_response`, `treats`, and `treatment_for` -> `pharmaco_effect_on`
- `part_of` -> `has_part`

Drop these aliases instead of mapping them to KG relations: `associated_with`, `associated`, `association`, `related`, `related_to`, `linked_to`, `relevant_to`, and `involved_in`.

Never hard-filter by target unless the target is grounded to a canonical KG entity, database ID, or high-confidence synonym. For natural-language targets such as `consensus molecular subtype`, `respiratory function in DMD`, `new profiles in JASPAR 2024`, clinical trial names, challenge names, or study-specific outcomes, use soft target mode: keep `target_type` when useful, but use the target phrase only for post-query ranking.

Use relaxed-to-specific retrieval:

1. Stage 1: source + source_type + optional target_type; no hard target and no unsafe relation.
2. Stage 2: if Stage 1 returns evidence, narrow by schema-compatible relation.
3. Stage 3: only add target if target grounding is high-confidence.

Low-yield questions should still perform a small KG probe, not skip KG. Empty KG results are not negative evidence. For list questions, KG should verify, normalize, or supplement candidates; `co_mentioned`-only candidates are weak evidence and normally should not enter the final answer.

## Evidence Packet Boundary

Do not pass raw KG rows directly into final answer generation. Convert them into an evidence packet first:

```text
KG retrieval -> evidence packet -> reasoning validation -> answer rendering
```

Evidence packet groups:

- `strong_curated`: answer-eligible evidence from curated/schema relations such as `causally_related_to`, `disease_relation`, `pharmaco_effect_on`, `regulated`, `produce`, `participate_in`, `interact_with`, `functionally_related_to`, and `is_a`.
- `weak_trace`: `co_mentioned`, `mentioned`, and `no_relation`; these may support explanation or tracing, but must not by themselves create a final answer candidate.
- `background`: legal KG rows whose relation, answer type, or question scope does not match the requested answer boundary.

Use `brick reason packet` or the `evidence_packet` module when a runner, agent, or downstream model needs prompt-ready KG context. The default benchmark KG policy is `curated_preferred` plus `snippet_supported`: curated/schema rows are preferred, ordinary `co_mentioned` rows remain weak trace/background, and only snippet-supported `co_mentioned` rows may be considered as tightly capped factoid/list supplements.

Evidence tiers:

- `E3 strong_curated`: curated/schema relation evidence that can become an answer candidate when answer type, scope, candidate anchor, and stage constraints match.
- `E2 snippet_supported_weak`: `co_mentioned` or `mentioned` evidence whose `original_relation` snippet explicitly mentions the candidate/entity, question anchor or scope, and a relation cue without negation or contradiction. This is final-answer eligible only for factoid/list-style tasks and only under caps.
- `E1 weak_trace`: ordinary `co_mentioned` or `mentioned` evidence. It can support explanation, tracing, or diagnostics but cannot create a final answer candidate.
- `E0 drop/background`: `no_relation`, irrelevant rows, contradiction rows, or rows that fail answer-boundary checks.

Default P2 policy:

```bash
--kg-relation-policy curated_preferred
--weak-evidence-policy snippet_supported
--co-mentioned-max-share 0.15
--co-mentioned-max-rows-per-question 3
--allow-weak-final-types factoid,list
```

Use `curated_only` or `no_co_mentioned` only as ablation modes. They are useful for proving how much weak literature trace contributes, but they are not the default D-group strategy.

For multiple-choice questions, query each option as a hard anchor. Keep the option target in the KG probe; do not perform `fallback_drop_option`. Only `strong_curated` option evidence should enter the choice prompt. If every option has weak-only or no KG evidence, fall back to native LLM reasoning without KG rows.

For MCQ questions with multiple targets, judge scope at the option-packet level, not the single-row level. Multiple strong rows for the same option may collectively cover all required targets; rows that cover only part of the target set are `partial_strong` diagnostics and must not enter the choice prompt. Use `query_relation_with_fallback(stop_condition="candidate_eligible")` with packet-level validation so weak/background rows do not stop fallback before answer-eligible evidence has been checked.

For MCQ evidence override, candidate eligibility must also respect fallback stage specificity. `initial`, `fallback_no_relation`, and `fallback_relax_target_type` may produce answer-candidate evidence because the target constraint is retained. `fallback_drop_target` and `fallback_drop_target_no_relation` are broad diagnostic stages; rows from these stages should be background/diagnostic unless a later task explicitly opts into broad retrieval. When multiple options have strong evidence, rank options by target coverage, relation priority, source reliability, path specificity, and fallback-stage specificity. If the best option is not separated from the second-best option by the configured margin, fall back to native LLM reasoning rather than overriding with KG evidence.

For BioASQ/list-style answers, generate direct answer candidates first, then use KG to verify, standardize synonyms, and add only small curated/direct supplements. `co_mentioned`-only rows cannot add new final list answers; KG-added list candidates should be capped, for example `max_kg_added_answers=3`.

For BioASQ factoid/list P2 runs, the runner should log KG utility at three levels: `kg_raw_rows`, `kg_used_rows`, and `kg_final_answer_sources`, plus raw/used/final `co_mentioned` shares. Ordinary `co_mentioned` final sources should be zero; only `snippet_supported_weak` rows may appear in final-source counts, and only under the factoid/list caps.

For BioHopR `PathAnswerSet`, use the skill-owned `balanced` preset by default:

```bash
brick biohopr-eval --preset balanced --split-group single
```

BioHopR is an answer-set retrieval benchmark, so the default differs from BioASQ/list and MCQ policies. Use the strict question-only setting for main results: inference receives only the natural-language question, extracts the `PathAnswerSet` contract itself, queries KG paths, scores paths, merges `native_plus_kg`, and renders answers. Do not use `gold_answers`, `gold_hops`, official bridge entities, or official query/target types during inference; gold answers are only for the evaluator after prediction. The balanced preset uses `query_source=extracted`, `native_plus_kg`, `kg_limit=1000`, `kg_add_min_score_1hop=2.0`, `kg_add_min_score_2hop=4.0`, and allows weak/co_mentioned KG supplements only through BioHopR path scoring and diagnostics. `--preset strict_e2e_question_only` is an explicit alias for the same main-result policy. Source-restricted runs such as `--preset sourcefiltered` are ablations, not the default D_main strategy.

## Benchmark Evaluation Entrypoints

Use the skill-owned evaluation commands for D-group/new-skill benchmark reporting:

```bash
brick biohopr-eval --preset balanced --split-group single
brick biomix-eval --method deepseek_v4_pro_D
brick biomix-mcq-realkg --method hermes_brick_source_aware_mcq_20260707 --source-routing source_aware
brick biomix-mcq-consensus --output-method biomix_mcq_question_only_consensus_20260707
brick bioasq-eval --input "D:\code\Biomics\BioASQ Benchmark\results_old\deepseek_v4_pro_D.results.json"
brick benchmark-suite --datasets all
brick biomni-kg-overlap --biomni-repo "D:\code\Biomics\_tmp_biomni_official" --live-probe
```

`biohopr-eval` can regenerate and evaluate BioHopR predictions with the strict question-only `PathAnswerSet` policy. `biomix-eval` scores existing Biomix answer files with gold used only by the evaluator after prediction. `bioasq-eval` is the unified BioASQ reporting wrapper: it combines the strict prior-guided self-correction YesNo metrics with legacy-comparable Factoid/List metrics, writes one `bioasq_unified.metrics.json`, and records task-level protocol/status fields. `rejudge_yesno_deepseek.py` is the YesNo generation helper used to create the historical `prior_guarded` artifact; it is not a separate BioASQ dataset evaluator or an independent rerun. BioASQ YesNo is no-gold during generation but uses the previous raw answer as a prior and is not held-out without a separate split; BioASQ Factoid/List are `legacy_goldaware` because their standardization path uses gold exact answers.

For Biomix MCQ, use hard option anchoring and question-only disease surface alias expansion from `scripts/biomix_mcq.py`. The alias expansion is allowed because it is derived from the natural-language question text; it must not read `mcq_biomix.csv`, `mcq_biomix_gold.json`, or any answer file during inference. Details should record `target_alias_groups`, `query_targets`, fallback rate, selected KG evidence rate, and whether the final choice came from deterministic KG ranking or question-only fallback.

`biomix-mcq-realkg` is the strict independent MCQ route for Hermes + BRICK + real KG. It disables native fallback files, Biomni primary, question-only consensus, and gold-aware standardization. The default `--source-routing source_aware` preset parses the question/options, infers the option/source type, queries high-value edge `info_source` groups first for each option as a hypothesis, broadens only through lower-weight fallback/network evidence, writes an all-option evidence table, and asks Hermes to choose only from that evidence table plus the question. Gold labels are used only later by `biomix-eval`.

`biomix-mcq-gwas-verifier` is an external-source supplemental MCQ route for Biomix disease/trait association questions, not a strict D-main BiOmics-KG-only route. It reads only the question/options, queries the public GWAS Catalog REST API for each option against the two question traits, and uses deepseek-v4-pro only as a low-evidence fallback. It must not read `mcq_biomix.csv`, Biomni outputs, historical prediction files, or consensus outputs during generation. Report it as external source-aware GWAS evidence, or as D-main only after the GWAS Catalog records have been ingested into a frozen BiOmics-KG snapshot and inference queries only that KG.

`biomix-mcq-consensus` is a separate question-only ensemble utility for MCQ ablations. It reads the no-gold Biomix question file plus existing prediction files, starts from a primary source, and overrides only when the configured secondary sources have strong consensus on another option. It does not read gold labels during generation, but it is not a pure D-group BRICK-only run if `Biomni` is used as the primary or as a source. Use `biomix-eval` only after generation.

`biomni-kg-overlap` compares Biomni database query tools from a local Biomni checkout with the BiOmics-KG edge `info_source` catalog. It reports source-name/alias coverage without reading benchmark gold labels; `--live-probe` optionally performs small read-only KG API probes for mapped sources.

## Task Contract Architecture

Do not encode KG strategy as dataset-specific benchmark logic. Route user questions and benchmark items through:

```text
TaskContract -> schema-aware retrieval -> EvidencePacket -> ReasoningPolicy -> Renderer -> EvaluationAdapter
```

Use `scripts/task_contract.py` to map datasets into common task contracts:

- `ClosedChoice`: Biomix MCQ and other option-selection tasks.
- `BinaryClaim`: Biomix true/false and BioASQ yes/no.
- `ShortEntity`: BioASQ factoid and short entity answers.
- `EntityList`: BioASQ list.
- `EntityProfile`: Biomix OpenQA or entity-card answers.
- `MechanisticSummary`: BioASQ summary-style mechanism answers.
- `PathAnswerSet`: BioHopR 1-hop/2-hop answer/path sets.

Use `scripts/renderers.py` only for output format adaptation, such as BioASQ nested answer arrays, Biomix option text, true/false strings, or BioHopR path payloads. Use `scripts/evaluation_adapters.py` for contract-level normalization and diagnostic comparisons; official benchmark metrics may still live in dataset runners.

The runner should become a thin shell around dataset IO, cache/checkpoint, retries/timeouts, and metric collection. It should not redefine relation schema, weak evidence policy, answer eligibility, RankGraph scoring, or reasoning policy.

## Multi-Hop Querying

`brick qr` is direct one-hop by default. Multi-hop is optional and capped at 3 edges:

```bash
brick qr "TP53" --st Gene --tt Disease -n 50 -o tp53_disease_direct.csv
brick qr "TP53" --st Gene --tt Disease --hops 2 -n 50 -o tp53_disease_2hop.csv
brick qr "TP53" --st Gene --tt Disease --hops 3 --directed -n 50 -o tp53_disease_3hop_directed.csv
```

Use multi-hop when:

- A direct edge is absent or too sparse, but the question asks for a plausible mechanism or bridge.
- The desired answer may be mediated through pathway, process, protein interaction, phenotype, cell type, or disease ontology nodes.
- A plan needs both a high-precision one-hop query and a recall-oriented path query for comparison.

Avoid multi-hop when:

- The question asks for a direct marker, drug indication, variant annotation, or curated source-specific fact.
- `co_mentioned` paths dominate and snippets do not directly support the claim.
- The query is already broad; multi-hop can expand noise quickly.

When creating a reusable plan, let `brick plan` decide this explicitly:

```bash
brick plan template --task query-rank -i genes.txt --source-type Gene --target-type Disease --hops 1 -o direct_plan.json
brick plan template --task query-rank -i genes.txt --source-type Gene --target-type Disease --hops 3 -o multihop_plan.json
```

For every multi-hop result, inspect intermediate nodes (`path.2`, `path.4`, etc.) and every edge relation/source before using the path as biological evidence.

## Workflow Scaffolds

`brick workflow` creates reproducible folders with inputs, intermediate paths, and shell scripts. It should not replace scientific judgment.

```bash
brick workflow celltype --marker-dir markers -o wf_celltype
brick workflow enrichment --genes geneset.txt -o wf_enrichment
brick workflow drug --genes geneset.txt --real-drugs known_drugs.csv -o wf_drug
brick workflow trajectory -g paga.graphml --query-csv kg_results.csv -o wf_trajectory
```

Use `--run-query` only when the BiOmics API is reachable. Use `--run-analysis` for local steps that can run from existing files.

## Block Plans

Use `brick plan` when the user needs a reusable model blueprint rather than an immediate one-off command.

```bash
brick plan blocks
brick plan schema
brick plan inspect query_cell.csv ranked_cell.csv merged.graphml
brick plan template --task query-rank -i genes.txt --source-type Gene --target-type Cell -r produce --output-dir model_run -o model_run/brick_plan.json
brick plan template --task query-rank -i genes.txt --source-type Gene --target-type Cell -r produce --agent-plan --output-dir model_run -o model_run/agent_plan.json
brick plan validate model_run/brick_plan.json
brick plan materialize model_run/brick_plan.json -o model_run
```

Treat the JSON plan as the machine-readable definition of a brick model instance: blocks, commands, inputs, outputs, and which steps require the remote KG API.

## Agent Orchestration

Use `brick agent` as a deterministic BioMaster-like execution shell around the KG-aware brick core. It should help a long analysis run steadily, but it must not replace KG source policy, relation filtering, bounded multi-hop querying, co_mentioned review, or evidence-strength reasoning.

P0 agent behavior is intentionally conservative:

- `brick agent check-plan` validates the plan before execution.
- `brick agent run` is dry-run by default and writes `run_manifest.jsonl`.
- `brick agent run --execute` only allows safe `brick ...` commands.
- `brick agent check` validates generated artifacts.
- `brick agent report` creates a lightweight markdown summary from the manifest and optional evidence CSV.

Recommended sequence:

```bash
brick agent check-plan -p model_run/agent_plan.json
brick agent run -p model_run/agent_plan.json --workdir model_run
brick agent run -p model_run/agent_plan.json --workdir model_run --execute --max-retries 1
brick agent report --manifest model_run/run_manifest.jsonl --evidence model_run/evidence_units_cell.csv -o model_run/agent_report.md
```

Agent validators are grouped by purpose:

| Validator Group | Examples | Purpose |
|---|---|---|
| Plan | `depends_on_satisfied`, `no_output_collision`, `kg_policy_not_broad`, `relation_alias_resolved`, `max_hops_respected` | Reject unsafe or under-scoped plans before execution |
| File | `file_exists`, `non_empty_file`, `csv_columns_present`, `json_schema_valid` | Check artifact existence and basic schema |
| KG | `kg_edge_metadata_columns_present`, `kg_edge_metadata_required_for_co_mentioned`, `co_mentioned_weak_cap` | Ensure edge-level provenance is reviewable and weak co-mentioned evidence stays weak |
| Evidence | `evidence_strength_present`, `manual_review_notes_present`, `no_weak_only_strong_claim` | Prevent unsupported biological claims from weak evidence |

For workflow scaffolds, request an agent plan directly:

```bash
brick workflow drug --genes genes.txt --agent-plan -o wf_drug
brick agent check-plan -p wf_drug/agent_plan.json
```

## Reasoning Contract

When the user asks "why", "is this credible", "what evidence supports it", or "compare with known results", use `brick reason` before writing conclusions.

- `brick reason filter-co-mentioned` removes unsupported `co_mentioned` rows from raw `brick qr` CSVs by reviewing `path.*.original_relation` before rank/explain steps.
- `brick reason normalize` converts brick query rows into canonical evidence units for model-like reasoning.
- `brick reason ground` canonicalizes entity names/IDs using an optional lexicon CSV and fuzzy fallback.
- `brick reason score` scores evidence using confidence, information-source count, relation priors, grounding quality, and snippet availability.
- `brick reason conflicts` surfaces negative or contradictory evidence.
- `brick reason evidence` summarizes matched sources, relations, confidence-like fields, and snippets from `brick qr` CSVs.
- `brick reason hypotheses` converts ranked/enriched outputs into reviewable hypothesis tables.
- `brick reason explain` creates an interpretability report from evidence units plus optional ranked/enriched hypotheses.
- `brick reason paths` enumerates simple paths in a merged GraphML graph.
- `brick reason path-score` ranks mechanism paths by edge support and path length.
- `brick reason store` persists evidence units in a local SQLite evidence store.
- `brick reason eval` evaluates predictions against reference entities with precision/recall/hit/MRR/AP and optional evidence coverage.
- `brick reason compare` marks overlap between predictions and reference entities.

For deeper evidence search in raw query CSVs, see `references/data_processing_patterns.md`.

## Co-Mentioned Edge Handling

`co_mentioned` means two entities were mentioned together in the same literature sentence or text-mined snippet. It is **not** a structured biological relation by itself. The actual meaning must be read from the source sentence in `path.*.original_relation` (for one-hop rows, usually `path.1.original_relation`).

Never treat `co_mentioned` alone as proof of causality, treatment effect, marker status, regulation, clinical utility, or mechanism.

### What co-mentioned rows mean

- `path.1.relation == co_mentioned` indicates co-occurrence in source text, not a curated KG relation such as `produce`, `regulated`, or `pharmaco_effect_on`.
- `path.1.original_relation` is the primary evidence sentence. Read this column before using the row.
- If the sentence only shows that two entities appear in the same context but does **not** state a biologically useful relationship for the user's question, **drop the row** from downstream query tables.

### Required review steps

When a result row contains `co_mentioned`:

1. Read `path.1.original_relation` and decide whether the sentence supports the biological claim being queried.
2. Read `path.1.info_source`, `path.1.info_source_length`, `path.1.condition`, and `path.1.relation_confidence` for context, but do not let source count override snippet relevance.
3. Check whether the answer entity appears in the snippet, not only as `path.2.name`. Some co-mentioned rows may have an incidental target label while the useful clue is inside the edge text.
4. If the snippet only shows broad co-occurrence, treat it as weak background support or exclude it.
5. If the snippet directly contains the requested fact, keep the row and cite the relevant source field.
6. If the snippet does not support the requested fact, **exclude the row** before ranking, hypothesis generation, or explanation.

### Automated filtering command

Use `brick reason filter-co-mentioned` immediately after `brick qr` when co-mentioned rows may appear:

```bash
brick qr -i genes.txt --st Gene --tt Disease -n 200 -o query.csv
brick reason filter-co-mentioned \
  -i query.csv \
  -o query_filtered.csv \
  --rejected-output query_comentioned_rejected.csv \
  --audit-output query_comentioned_audit.csv
brick rk -i query_filtered.csv --rank -o ranked.csv
```

Default behavior in `strict` mode:

- Keep all non-`co_mentioned` structured relation rows unchanged.
- For `co_mentioned` rows, require a non-empty `path.*.original_relation` snippet.
- Require both endpoints to appear in the snippet unless `--allow-missing-endpoint` is set.
- Require a surface relationship cue in the sentence (for example marker, regulate, associate, inhibit, treat, express, interact) unless `--keep-bare-co-occurrence` is set.
- Drop unresolved or metadata-only co-mentioned rows.
- Write dropped rows to `--rejected-output` when provided.

Use `--mode permissive` only when you want to drop clearly unusable co-mentioned rows but keep bare co-occurrence sentences for manual review.

`brick reason normalize`, `brick reason score`, and `brick reason evidence` still classify remaining co-mentioned rows from edge metadata and emit `edge_support_class`, `edge_detail_score`, and `edge_review_note`. Filtering happens **before** those steps when you use `filter-co-mentioned`.

Recommended review pattern for BioASQ-style questions:

1. Query with the most specific grounded entity.
2. Run `brick reason filter-co-mentioned` on the raw query CSV.
3. If only filtered `co_mentioned` rows remain, inspect snippets again by answer terms, synonyms, and question terms.
4. Rank candidate rows by snippet relevance first, then `info_source_length`, then evidence score.
5. Synthesize an answer only after the edge text has been inspected.

## Validation Rules

Before reporting a result as biologically meaningful:

- Check that source and target entity types match the question.
- Check that relation filters did not remove essential evidence or introduce overly broad parent terms.
- Inspect at least one raw `brick qr` CSV and, when available, `path.1.original_relation`.
- For multi-hop rows, inspect every edge and intermediate node in the path. Treat the whole path as weak unless each step is biologically meaningful and source-supported.
- For `co_mentioned` rows, read `path.1.original_relation`, exclude unsupported rows with `brick reason filter-co-mentioned`, and only then rank or explain the filtered query table.
- For ranking, report the ranking metric and whether lower or higher is better.
- For enrichment, report p-value/FDR columns if present.
- For graph outputs, inspect node and edge counts with `brick ga info` before plotting.
- For drug or benchmark tasks, freeze the primary ranking before reference access; compare by stable ID and label reused references as post hoc exploratory evidence.

## Environment Notes

On Windows PowerShell, prefer the generated `commands.ps1` files. If JSON-like graph attributes cause quoting trouble, first run the simpler command without `--node-attrs` or `--edge-attrs`, then post-process the GraphML.

The BiOmics API can be affected by proxy settings. By default, `brick qr` should avoid inheriting environment proxies unless `BRICK_TRUST_ENV_PROXY=1` is set. The legacy alias `BRICK_API_TRUST_ENV_PROXY=1` is also accepted. If API calls time out, verify network access separately and continue with offline bricks where possible.

The standard test suite avoids live KG calls. To run the optional compact API schema smoke test, set:

```bash
RUN_BRICK_API_TESTS=1 pytest tests/test_live_api_schema.py
```

## Quick Help

```bash
brick qr --help
brick rk --help
brick ga --help
brick emb --help
brick drug --help
brick-gpu drug grace-rank --help
brick pl --help
brick workflow --help
brick reason --help
brick plan --help
```

## Evidence strength rules added in optimized fix

- `co_mentioned`, co-occurrence, text-mined, literature-mined, or abstract-only edges are normalized as `evidence_strength=weak_co_mentioned_review_required`.
- Weak evidence is capped during scoring/ranking and is marked `review_required=True`.
- Reports must surface `edge_support_class`, `evidence_strength`, `source_quality`, `condition`, `relation_confidence`, and `edge_review_note` when present.
- Multi-hop query outputs are normalized one edge per hop, preserving `path_id`, `path_hops`, `edge_index`, and `edge_position`.
