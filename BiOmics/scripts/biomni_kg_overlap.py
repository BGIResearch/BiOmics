#!/usr/bin/env python3
"""Compare Biomni database tools with BiOmics-KG edge info_source catalog."""

from __future__ import annotations

import argparse
import ast
import csv
import json
import re
from pathlib import Path
from typing import Iterable


STRUCTURE_OR_EXTERNAL = {
    "alphafold",
    "emdb",
    "iucn",
    "paleobiology",
    "pdb",
    "pride",
    "worms",
}


BIOMNI_TO_KG_ALIASES = {
    "cbioportal": ["cBioPortal"],
    "chembl": ["ChEMBL"],
    "clinicaltrials": ["ClinicalTrials", "ClinicalTrails", "ClinicalTrials.gov"],
    "clinvar": ["ClinVar", "ClinVar_gene"],
    "dailymed": ["DailyMed"],
    "dbsnp": ["dbSNP"],
    "encode": ["ENCODE"],
    "ensembl": ["ENSEMBL"],
    "geo": ["GEO"],
    "gnomad": ["gnomAD"],
    "gtopdb": ["IUPHAR"],
    "gwas_catalog": ["GWAS_Catalog", "GWAS Catalog", "GWAS_cataog"],
    "interpro": ["InterPro"],
    "jaspar": ["JASPAR"],
    "kegg": ["KEGG", "KEGG_anno", "KEGG_anno_mm"],
    "monarch": ["MONDO", "HPO", "DOID"],
    "mpd": ["MPD"],
    "openfda": ["OpenFDA"],
    "opentarget": ["OpenTargets"],
    "pubchem": ["PUBCHEMCID"],
    "quickgo": ["GO", "GOA", "GOA_mm", "GO-PLUS"],
    "reactome": ["Reactome", "Reactome_anno"],
    "regulomedb": ["RegulomeDB"],
    "remap": ["ReMap"],
    "stringdb": ["StringDB", "StingDB", "StingDB_mm"],
    "synapse": ["Synapse"],
    "ucsc": ["UCSC"],
    "unichem": ["UniChem"],
    "uniprot": ["UNIPROT_human", "UNIPROT_rodents", "UNIPROT_mapping"],
}


LIVE_PROBES = [
    ("clinvar", "ClinVar", "TP53", "Gene", "Disease", "causally_related_to"),
    ("clinvar_gene", "ClinVar_gene", "TP53", "Gene", "Disease", "causally_related_to"),
    ("clinicaltrials", "ClinicalTrails", "imatinib", "Chemical", "Disease", None),
    ("ensembl", "ENSEMBL", "TP53", "Gene", "Protein", "sequentially_related_to"),
    ("gtopdb", "IUPHAR", "sildenafil", "Chemical", "Protein", None),
    ("kegg", "KEGG", "TP53", "Gene", "Pathway", "participate_in"),
    ("kegg_anno", "KEGG_anno", "TP53", "Gene", "Pathway", "participate_in"),
    ("monarch_mondo", "MONDO", "breast cancer", "Disease", "Disease", "is_a"),
    ("monarch_doid", "DOID", "breast cancer", "Disease", "Disease", "is_a"),
    ("pubchem", "PUBCHEMCID", "sildenafil", "Chemical", None, None),
    ("quickgo_goa", "GOA", "TP53", "Gene", "Process", "participate_in"),
    ("quickgo_go", "GO", "apoptotic process", "Process", "Process", "is_a"),
    ("reactome", "Reactome", "TP53", "Gene", "Pathway", "participate_in"),
    ("reactome_anno", "Reactome_anno", "TP53", "Gene", "Pathway", "participate_in"),
    ("stringdb", "StringDB", "TP53", "Protein", "Protein", "interact_with"),
    ("uniprot_human", "UNIPROT_human", "TP53", "Protein", None, None),
    ("uniprot_mapping", "UNIPROT_mapping", "TP53", "Gene", "Protein", "sequentially_related_to"),
]


def _norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def _add_source(sources: set[str], value: str) -> None:
    value = (value or "").strip()
    if value and re.match(r"^[A-Za-z0-9_+\-. /]+$", value):
        sources.add(value)


def _walk_policy_sources(node, sources: set[str]) -> None:
    if isinstance(node, dict):
        for key, value in node.items():
            if key in {"info_sources", "primary_info_sources", "fallback_info_sources"} and isinstance(value, list):
                for item in value:
                    _add_source(sources, str(item))
            else:
                _walk_policy_sources(value, sources)
    elif isinstance(node, list):
        for item in node:
            _walk_policy_sources(item, sources)


def load_kg_sources(skill_dir: Path) -> set[str]:
    sources: set[str] = set()
    policy_path = skill_dir / "references" / "kg_source_policy.json"
    if policy_path.exists():
        _walk_policy_sources(json.loads(policy_path.read_text(encoding="utf-8")), sources)

    md_path = skill_dir / "references" / "kg_data_sources.md"
    if md_path.exists():
        md = md_path.read_text(encoding="utf-8", errors="ignore")
        for line in md.splitlines():
            match = re.match(r"^\|\s*([A-Za-z][A-Za-z0-9_+\-. /]+?)\s*\|", line)
            if not match:
                continue
            value = match.group(1).strip()
            if value not in {"Source", "Node type", "Entity type", "Resource", "Edge type"}:
                _add_source(sources, value)
    return sources


def load_biomni_databases(biomni_repo: Path) -> list[str]:
    names: set[str] = set()
    schema_dir = biomni_repo / "biomni" / "tool" / "schema_db"
    if schema_dir.exists():
        names.update(path.stem for path in schema_dir.glob("*.pkl"))

    database_py = biomni_repo / "biomni" / "tool" / "database.py"
    if database_py.exists():
        module = ast.parse(database_py.read_text(encoding="utf-8"))
        for node in module.body:
            if isinstance(node, ast.FunctionDef) and node.name.startswith("query_"):
                name = node.name.removeprefix("query_")
                if name != "pdb_identifiers":
                    names.add(name)
    return sorted(names)


def build_overlap_rows(biomni_names: Iterable[str], kg_sources: set[str]) -> list[dict[str, str]]:
    norm_to_source = {_norm(source): source for source in kg_sources}
    rows = []
    for name in sorted(biomni_names):
        mapped = BIOMNI_TO_KG_ALIASES.get(name, [])
        present = []
        for candidate in mapped:
            if candidate in kg_sources:
                present.append(candidate)
                continue
            normalized = norm_to_source.get(_norm(candidate))
            if normalized:
                present.append(normalized)
        if present:
            status = "covered"
        elif name in STRUCTURE_OR_EXTERNAL:
            status = "not_applicable_structure_or_external"
        else:
            status = "not_in_catalog"
        rows.append(
            {
                "biomni_database": name,
                "status": status,
                "mapped_kg_names": ";".join(mapped),
                "present_kg_names": ";".join(present),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else [])
        writer.writeheader()
        writer.writerows(rows)


def run_live_probe() -> list[dict[str, str]]:
    from querygraph import query_relation

    rows = []
    for name, info_source, entity, source_type, target_type, relation in LIVE_PROBES:
        try:
            df = query_relation(
                entity,
                relation=relation,
                source_entity_type=source_type,
                target_entity_type=target_type,
                relation_info_source=[info_source],
                limit_number=3,
                return_type="dataframe",
            )
            first_relation = ""
            if not df.empty and "path.1.relation" in df.columns:
                first_relation = str(df.iloc[0].get("path.1.relation") or "")
            rows.append(
                {
                    "probe": name,
                    "info_source": info_source,
                    "rows": str(len(df)),
                    "status": "ok",
                    "first_relation": first_relation,
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "probe": name,
                    "info_source": info_source,
                    "rows": "0",
                    "status": f"error:{type(exc).__name__}:{str(exc)[:120]}",
                    "first_relation": "",
                }
            )
    return rows


def _print_rows(rows: list[dict[str, str]], max_rows: int | None = None) -> None:
    shown = rows if max_rows is None else rows[:max_rows]
    for row in shown:
        print(json.dumps(row, ensure_ascii=False))


def _cli_main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--biomni-repo", required=True, help="Path to a local Biomni checkout.")
    parser.add_argument(
        "--brick-skill-dir",
        default=str(Path(__file__).resolve().parents[1]),
        help="Path to BiOmics_BRICK_SKILL_fix.",
    )
    parser.add_argument("--output", help="Optional CSV path for Biomni-vs-KG source overlap.")
    parser.add_argument("--live-probe", action="store_true", help="Run small BRICK KG API probes for mapped sources.")
    parser.add_argument("--probe-output", help="Optional CSV path for live probe rows.")
    args = parser.parse_args()

    biomni_repo = Path(args.biomni_repo).resolve()
    skill_dir = Path(args.brick_skill_dir).resolve()
    biomni_names = load_biomni_databases(biomni_repo)
    kg_sources = load_kg_sources(skill_dir)
    rows = build_overlap_rows(biomni_names, kg_sources)

    summary = {
        "biomni_database_count": len(biomni_names),
        "kg_catalog_source_count": len(kg_sources),
        "covered_or_alias": sum(row["status"] == "covered" for row in rows),
        "not_applicable_structure_or_external": sum(row["status"] == "not_applicable_structure_or_external" for row in rows),
        "not_in_catalog": sum(row["status"] == "not_in_catalog" for row in rows),
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    _print_rows(rows)

    if args.output:
        write_csv(Path(args.output), rows)
        print(f"overlap_csv={args.output}")

    if args.live_probe:
        probe_rows = run_live_probe()
        probe_summary = {
            "probe_count": len(probe_rows),
            "probe_hit_count": sum(int(row["rows"]) > 0 for row in probe_rows if row["rows"].isdigit()),
        }
        print(json.dumps(probe_summary, ensure_ascii=False, indent=2))
        _print_rows(probe_rows)
        if args.probe_output:
            write_csv(Path(args.probe_output), probe_rows)
            print(f"probe_csv={args.probe_output}")


if __name__ == "__main__":
    _cli_main()
