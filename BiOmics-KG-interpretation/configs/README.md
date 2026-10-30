# BiOmics-BRICK Benchmark Configs

This folder stores configuration files for benchmark-specific BRICK policies.

## Files

- `biomix_mcq_realkg_policy.yaml`: source-aware option-verification policy for
  Biomix MCQ experiments.

## Usage

The YAML policy is documentation plus a reproducible preset for MCQ runs that
treat each option as a hypothesis and score it with KG evidence. The main
uploaded D result uses the script-level arguments documented in:

```text
evaluate_pipeline/skill_benchmark_eval/D_main_hopaware/D_GROUP_NEW_RESULTS_20260707.md
```

The supplemental GWAS run is not configured here because it calls the external
GWAS Catalog REST API through:

```text
skills/biomics_brick_fix/scripts/biomix_mcq_gwas_verifier.py
```

## Protocol Note

Config files may route KG source families, weak-evidence policy, and answer
selection behavior. They must not contain gold answers, official BioHopR hops,
official bridges, or benchmark labels used during generation.
