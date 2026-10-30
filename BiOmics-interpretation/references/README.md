# BiOmics-BRICK References

This folder contains reference documentation for KG source selection and
task-specific BRICK workflows.

## Key Files for Benchmarks

- `kg_data_sources.md`: KG source catalog, canonical source names, common
  misspellings, and source-routing guidance.
- `kg_source_policy.json`: machine-readable source routing and policy defaults.
- `kg_schema.json`: node types, relation types, and schema constraints.
- `D_main_policy.yaml`: D-main benchmark policy notes.

## Task Recipes

- `celltype_annotation_cli.md`
- `trajectory_analysis_cli.md`
- `gene_enrichment.md`
- `drug_discovery_cli.md`
- `gwas_causal_snps_cli.md`
- `gwas_phenotype_prediction_cli.md`
- `proteome_analysis_cli.md`
- `DEG_analysis_cli.md`
- `grn_analysis_cli.md`

## Source-Name Hygiene

Use `kg_data_sources.md` as the human-readable source of truth for names such
as `GWAS_Catalog` / `GWAS Catalog`, `StringDB`, `ClinVar`, `GOA`, `Reactome`,
and other KG edge `info_source` values. Benchmark code should canonicalize
source names before source-aware routing and should not silently introduce new
spellings without updating the policy/reference files.
