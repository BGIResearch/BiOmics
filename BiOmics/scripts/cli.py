#!/usr/bin/env python3
"""Unified CLI entry point for BiOmics brick tools."""

import importlib
import sys


COMMANDS = {
    "qr": ("querygraph", "_cli_main", "Query entity relations from the knowledge graph"),
    "rk": ("rankgraph", "_cli_main", "Rank or enrich brick query results"),
    "ga": ("graphadapter", "_cli_main", "Graph adapter: read, convert, save, and combine"),
    "emb": ("embedgraph", "_cli_main", "Embed an integrated omics/KG graph with UMAP or GRACE"),
    "pl": ("plotting.cli", "main", "Visualization tools for brick graphs"),
    "reason": ("reasoning", "_cli_main", "Evidence normalization, scoring, and reasoning helpers"),
    "plan": ("brickplan", "_cli_main", "Inspect, validate, and materialize brick block plans"),
    "workflow": ("workflow", "_cli_main", "Generate curated workflow scaffolds"),
    "policy": ("kgpolicy", "_cli_main", "Infer controlled KG source/hop query policies"),
    "agent": ("agent", "_cli_main", "Check, dry-run, execute, and report brick agent plans"),
    "drug": ("drugdiscovery", "_cli_main", "Stable-ID-first drug ranking, GRACE, and validation helpers"),
    "biohopr-eval": ("biohopr_eval", "_cli_main", "Run BioHopR benchmark evaluation with skill-owned PathAnswerSet presets"),
    "bioasq-eval": ("bioasq_eval", "_cli_main", "Run unified BioASQ YesNo/Factoid/List evaluation with task protocol labels"),
    "biomix-eval": ("biomix_eval", "_cli_main", "Evaluate existing Biomix D-group/new-skill outputs"),
    "biomix-mcq-realkg": ("biomix_mcq_realkg", "_cli_main", "Generate Biomix MCQ strict Hermes+BRICK real-KG predictions"),
    "biomix-mcq-gwas-verifier": ("biomix_mcq_gwas_verifier", "_cli_main", "Generate strict Biomix MCQ GWAS Catalog option-verifier predictions"),
    "benchmark-suite": ("benchmark_suite_eval", "_cli_main", "Evaluate BioHopR, BioASQ, and Biomix D-group outputs"),
    "biomni-kg-overlap": ("biomni_kg_overlap", "_cli_main", "Compare Biomni database tools with BiOmics-KG info_source catalog"),
}


def _print_help():
    print("Usage: brick <command> [options]")
    print()
    print("Commands:")
    for name, (_, _, description) in COMMANDS.items():
        print(f"  {name:<10} {description}")
    print()
    print("Run 'brick <command> --help' for command-specific help.")


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        _print_help()
        sys.exit(0)

    command = sys.argv[1]
    if command not in COMMANDS:
        print(f"Error: unknown command '{command}'", file=sys.stderr)
        print(file=sys.stderr)
        _print_help()
        sys.exit(1)

    module_name, entry_name, _ = COMMANDS[command]
    module = importlib.import_module(module_name)
    entry = getattr(module, entry_name)
    sys.argv = [f"brick {command}"] + sys.argv[2:]
    entry()


if __name__ == "__main__":
    main()
