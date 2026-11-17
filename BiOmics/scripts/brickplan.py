"""brick plan - inspect, describe, and compose brick analysis blocks."""

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

try:
    from kgpolicy import KG_POLICY_CHOICES, query_options_from_policy, select_kg_policy
except Exception:  # pragma: no cover
    KG_POLICY_CHOICES = ("manual", "auto", "off", "direct", "source-scoped", "multihop")
    query_options_from_policy = None
    select_kg_policy = None


ENTITY_TYPES = [
    "Cell",
    "Gene",
    "Protein",
    "Disease",
    "Phenotype",
    "Chemical",
    "Pathway",
    "Function",
    "Process",
    "Cell_Component",
    "Tissue",
    "Species",
    "Stage",
    "Mutation",
    "Publication",
]

DEFAULT_DRUG_GRACE_TARGET_TYPES = (
    "Disease,Phenotype,Process,Cell_Component,Pathway,Function,"
    "Tissue,Cell,Chemical,Stage"
)


RELATIONS = [
    {"name": "produce", "typical_source": "Gene", "typical_target": "Cell", "use": "cell marker and annotation"},
    {"name": "causally_related_to", "typical_source": "Gene", "typical_target": "Disease", "use": "disease association lookup"},
    {"name": "interact_with", "typical_source": "Gene/Protein", "typical_target": "Gene/Protein", "use": "interaction evidence"},
    {"name": "regulated", "typical_source": "Gene/Protein", "typical_target": "Gene", "use": "regulatory evidence"},
    {"name": "participate_in", "typical_source": "Gene/Protein", "typical_target": "Pathway/Process", "use": "function and pathway context"},
    {"name": "is_a", "typical_source": "Any", "typical_target": "Any", "use": "ontology parent-child structure"},
    {"name": "develope_into", "typical_source": "Cell", "typical_target": "Cell", "use": "cell trajectory and differentiation"},
    {"name": "pharmaco_effect_on", "typical_source": "Chemical", "typical_target": "Disease/Gene", "use": "known drug effect evidence"},
    {"name": "causally_related_to", "typical_source": "Gene", "typical_target": "Disease/Phenotype", "use": "causal hypothesis evidence"},
    {"name": "functionally_related_to", "typical_source": "Gene/Protein", "typical_target": "Function", "use": "functional interpretation"},
]


BLOCKS = [
    {
        "name": "query",
        "command": "brick qr",
        "category": "query",
        "requires_api": True,
        "inputs": ["entity list or source entity"],
        "outputs": ["query CSV"],
        "purpose": "retrieve raw entity-relation-entity evidence from BiOmics-KG, direct by default with optional <=3-hop path search",
    },
    {
        "name": "rank",
        "command": "brick rk --rank",
        "category": "rank",
        "requires_api": False,
        "inputs": ["brick qr CSV"],
        "outputs": ["ranked CSV"],
        "purpose": "prioritize target entities by match and evidence metrics",
    },
    {
        "name": "enrich",
        "command": "brick rk --enrich",
        "category": "rank",
        "requires_api": False,
        "inputs": ["brick qr CSV"],
        "outputs": ["enrichment CSV"],
        "purpose": "test enriched target entities for a source set",
    },
    {
        "name": "data-to-graph",
        "command": "brick ga data2nxg",
        "category": "graph",
        "requires_api": False,
        "inputs": ["CSV, AnnData, or PAGA AnnData"],
        "outputs": ["graph file"],
        "purpose": "turn tabular or omics structure into a graph artifact",
    },
    {
        "name": "complete-graph",
        "command": "brick ga complete_result",
        "category": "graph",
        "requires_api": False,
        "inputs": ["observed graph", "query CSV or graph"],
        "outputs": ["merged graph"],
        "purpose": "augment an observed graph with knowledge-graph evidence",
    },
    {
        "name": "plot",
        "command": "brick pl",
        "category": "visualize",
        "requires_api": False,
        "inputs": ["graph or analysis CSV"],
        "outputs": ["HTML, PNG, PDF"],
        "purpose": "visualize networks, enrichment results, or TF regulons",
    },
    {
        "name": "reason",
        "command": "brick reason",
        "category": "reason",
        "requires_api": False,
        "inputs": ["query CSV, ranked CSV, enrichment CSV, graph, reference set"],
        "outputs": ["evidence summaries, hypotheses, paths, comparisons"],
        "purpose": "convert brick artifacts into reviewable evidence and hypotheses",
    },
    {
        "name": "normalize-evidence",
        "command": "brick reason normalize",
        "category": "reason",
        "requires_api": False,
        "inputs": ["brick qr CSV"],
        "outputs": ["canonical evidence-unit CSV or JSONL"],
        "purpose": "standardize raw query rows into the evidence substrate for interpretability models",
    },
    {
        "name": "ground-entities",
        "command": "brick reason ground",
        "category": "reason",
        "requires_api": False,
        "inputs": ["evidence-unit table", "optional lexicon CSV"],
        "outputs": ["evidence-unit table with canonical IDs/names"],
        "purpose": "standardize entity mentions before scoring, storage, and evaluation",
    },
    {
        "name": "score-evidence",
        "command": "brick reason score",
        "category": "reason",
        "requires_api": False,
        "inputs": ["evidence-unit table"],
        "outputs": ["scored evidence or aggregated support table"],
        "purpose": "score support strength from confidence, source count, relation prior, grounding, and snippets",
    },
    {
        "name": "detect-conflicts",
        "command": "brick reason conflicts",
        "category": "reason",
        "requires_api": False,
        "inputs": ["scored or raw evidence-unit table"],
        "outputs": ["conflict and negative-evidence table"],
        "purpose": "surface contradictory or negative evidence instead of only positive support",
    },
    {
        "name": "score-paths",
        "command": "brick reason path-score",
        "category": "reason",
        "requires_api": False,
        "inputs": ["GraphML graph", "optional evidence-unit table"],
        "outputs": ["mechanism path score table"],
        "purpose": "rank explanatory graph paths by edge evidence strength and path length",
    },
    {
        "name": "evidence-store",
        "command": "brick reason store",
        "category": "reason",
        "requires_api": False,
        "inputs": ["evidence-unit table"],
        "outputs": ["SQLite evidence store"],
        "purpose": "persist local evidence units for reusable offline interpretability",
    },
    {
        "name": "evaluate",
        "command": "brick reason eval",
        "category": "reason",
        "requires_api": False,
        "inputs": ["prediction table", "reference table", "optional evidence-unit table"],
        "outputs": ["precision, recall, hit, MRR, AP, and evidence-coverage metrics"],
        "purpose": "benchmark biological hypotheses against known references",
    },
    {
        "name": "explain",
        "command": "brick reason explain",
        "category": "reason",
        "requires_api": False,
        "inputs": ["evidence-unit table", "optional ranked or enriched hypotheses"],
        "outputs": ["interpretability report"],
        "purpose": "trace hypotheses back to sources, relations, metrics, snippets, and provenance",
    },
    {
        "name": "merge-csv",
        "command": "brick plan merge-csv",
        "category": "utility",
        "requires_api": False,
        "inputs": ["one or more CSV files"],
        "outputs": ["merged CSV"],
        "purpose": "join query outputs before ranking or enrichment",
    },
]


QUERY_COLUMNS = {"path.0.name", "path.2.name"}
RANK_COLUMNS = {"path.2.rank_voting"}
ENRICH_COLUMNS = {"path.2.enrich_pvalue"}
RELATION_COLUMNS = {"path.1.relation", "path.1.name", "relation"}


def _utc_now():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _write_json(data, output=None):
    text = json.dumps(data, indent=2, ensure_ascii=False)
    if output:
        Path(output).write_text(text + "\n", encoding="utf-8")
        print(f"Wrote {output}")
    else:
        print(text)


def _print_table(rows, columns):
    widths = {col: len(col) for col in columns}
    for row in rows:
        for col in columns:
            widths[col] = max(widths[col], len(str(row.get(col, ""))))
    header = "  ".join(col.ljust(widths[col]) for col in columns)
    print(header)
    print("  ".join("-" * widths[col] for col in columns))
    for row in rows:
        print("  ".join(str(row.get(col, "")).ljust(widths[col]) for col in columns))


def _command_to_ps(command):
    return subprocess.list2cmdline([str(x) for x in command])


def _command_to_sh(command):
    return " ".join(shlex.quote(str(x)) for x in command)


def _block(block_id, brick, command, inputs=None, outputs=None, requires_api=False, notes=None):
    return {
        "id": block_id,
        "brick": brick,
        "requires_api": bool(requires_api),
        "command": [str(x) for x in command],
        "inputs": [str(x) for x in (inputs or [])],
        "outputs": [str(x) for x in (outputs or [])],
        "notes": notes or "",
    }


def _step_type_for_block(block):
    brick = block.get("brick", "")
    if brick == "query":
        return "kg_query"
    if brick == "normalize-evidence":
        return "evidence_normalization"
    if brick == "rank":
        return "ranking"
    if brick == "enrich":
        return "enrichment"
    if brick in {"complete-graph", "data-to-graph", "graph-info"}:
        return "graph"
    if brick == "plot":
        return "visualization"
    if brick == "explain":
        return "report"
    if brick == "reason":
        return "reasoning"
    return "utility"


def _default_agent_validators(block):
    brick = block.get("brick", "")
    common = ["file_exists", "non_empty_file"]
    if brick == "query":
        return [
            "kg_policy_not_broad",
            "max_hops_respected",
            *common,
            "kg_edge_metadata_columns_present",
            "kg_edge_metadata_required_for_co_mentioned",
        ]
    if brick == "normalize-evidence":
        return [
            *common,
            "evidence_strength_present",
            "manual_review_notes_present",
            "co_mentioned_weak_cap",
            "max_hops_respected",
        ]
    if brick == "explain":
        return [*common, "no_weak_only_strong_claim"]
    return common


def _infer_agent_dependencies(blocks):
    output_owner = {}
    dependencies = {}
    for block in blocks:
        block_id = block.get("id")
        deps = set()
        for input_path in block.get("inputs", []) or []:
            owner = output_owner.get(str(input_path))
            if owner and owner != block_id:
                deps.add(owner)
        dependencies[block_id] = sorted(deps)
        for output_path in block.get("outputs", []) or []:
            output_owner[str(output_path)] = block_id
    return dependencies


def to_agent_plan(plan):
    """Convert a block plan into the P0 deterministic agent plan contract."""
    dependencies = _infer_agent_dependencies(plan.get("blocks", []))
    steps = []
    for block in plan.get("blocks", []):
        step_type = _step_type_for_block(block)
        validators = _default_agent_validators(block)
        retry_count = 1 if block.get("requires_api") else 0
        steps.append(
            {
                "id": block.get("id"),
                "type": step_type,
                "command": block.get("command", []),
                "inputs": block.get("inputs", []),
                "outputs": block.get("outputs", []),
                "depends_on": dependencies.get(block.get("id"), []),
                "validators": validators,
                "validator_args": {"max_hops_respected": {"max_hops": 3}},
                "retry_policy": {
                    "max_retries": retry_count,
                    "retry_on": ["command_failed", "validator_failed"],
                },
                "checkpoint": bool(block.get("requires_api") or step_type in {"kg_query", "report", "reasoning"}),
                "allow_failure": False,
                "notes": block.get("notes", ""),
            }
        )
    return {
        "schema_version": "1.0",
        "task": plan.get("title", "brick agent plan"),
        "description": plan.get("description", ""),
        "question": plan.get("question", ""),
        "workdir": plan.get("output_dir", "."),
        "execution": {
            "mode": "dry_run",
            "executor": "brick",
            "allow_shell": False,
            "max_retries": 1,
            "timeout_seconds": 600,
        },
        "steps": steps,
    }


def _cmd(command, *args):
    return ["brick", command, *args]


def _base_plan(title, description, out_dir):
    return {
        "version": 1,
        "title": title,
        "description": description,
        "created_at": _utc_now(),
        "output_dir": str(out_dir),
        "blocks": [],
    }


def _split_csv_arg(value):
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def _ensure_out_dir(args):
    out_dir = Path(args.output_dir or "brick_plan_output")
    return out_dir


def _query_hops(args) -> int:
    return max(1, min(int(getattr(args, "hops", 1) or 1), 3))


def _resolve_policy(args, task=None, source_type=None, target_type=None, relation=None):
    if select_kg_policy is None:
        return {}
    return select_kg_policy(
        question=getattr(args, "question", None),
        source_type=source_type or getattr(args, "source_type", None),
        target_type=target_type or getattr(args, "target_type", None),
        relation=relation or getattr(args, "relation", None),
        info_source=getattr(args, "info_source", None),
        hops=_query_hops(args),
        task=task or getattr(args, "task", None),
        kg_policy=getattr(args, "kg_policy", "auto"),
        policy_file=getattr(args, "policy_file", None),
    )


def _policy_relation(args, task=None, source_type=None, target_type=None, relation=None):
    """Return explicit/hardcoded relation or a policy-derived relation string."""
    if relation:
        return relation
    if getattr(args, "relation", None):
        return args.relation
    policy = _resolve_policy(args, task=task, source_type=source_type, target_type=target_type, relation=relation)
    rels = policy.get("relations") or []
    return ",".join(rels) if rels else None


def _query_options(args, task=None, source_type=None, target_type=None, relation=None):
    options = []
    policy = _resolve_policy(args, task=task, source_type=source_type, target_type=target_type, relation=relation)
    explicit_info_source = getattr(args, "info_source", None)
    if explicit_info_source:
        options.extend(["--info-source", explicit_info_source])
    elif policy.get("primary_info_sources"):
        options.extend(["--info-source", ",".join(policy["primary_info_sources"])])
    hops = int(policy.get("hops") or _query_hops(args))
    if hops > 1:
        options.extend(["--hops", str(hops)])
    if getattr(args, "directed", False):
        options.append("--directed")
    if getattr(args, "kg_policy", "auto") != "auto":
        options.extend(["--kg-policy", getattr(args, "kg_policy")])
    if getattr(args, "question", None):
        options.extend(["--question", getattr(args, "question")])
    if getattr(args, "policy_file", None):
        options.extend(["--policy-file", getattr(args, "policy_file")])
    return options


def _query_notes(args, task=None, source_type=None, target_type=None, relation=None):
    policy = _resolve_policy(args, task=task, source_type=source_type, target_type=target_type, relation=relation)
    notes = []
    if policy:
        notes.append(
            f"KG policy intent={policy.get('intent')} policy={policy.get('kg_policy')} reason={policy.get('reason')}"
        )
        if policy.get("primary_info_sources"):
            notes.append("Use edge info_source filter: " + ",".join(policy["primary_info_sources"]))
    hops = int(policy.get("hops") or _query_hops(args))
    if hops > 1:
        direction = "directed" if getattr(args, "directed", False) else "undirected"
        notes.append(
            f"Optional multi-hop KG query enabled: up to {hops} hops ({direction}). Inspect intermediate nodes and every edge before treating the path as biological evidence."
        )
    notes.append("co_mentioned/text-mined edges must remain weak evidence unless original_relation/info_source/condition/confidence support a concrete claim.")
    return " ".join(notes)


def _append_relation_if_needed(command, relation):
    if relation and "-r" not in command and "--relation" not in command:
        command.extend(["-r", relation])
    return command

def _query_rank_plan(args):
    out_dir = _ensure_out_dir(args)
    source = args.input or args.genes
    if not source:
        raise ValueError("--input or --genes is required for query-rank.")
    target_type = args.target_type or "Cell"
    query_csv = out_dir / f"query_{target_type.lower()}.csv"
    evidence_csv = out_dir / f"evidence_units_{target_type.lower()}.csv"
    ranked_csv = out_dir / f"ranked_{target_type.lower()}.csv"
    explanation_md = out_dir / f"explanation_{target_type.lower()}.md"

    plan = _base_plan(
        "brick query-rank-reason plan",
        "Query BiOmics-KG, rank target candidates, and summarize evidence.",
        out_dir,
    )
    query_cmd = _cmd("qr", "-i", source, "--st", args.source_type, "--tt", target_type, "-n", args.query_limit, "-o", query_csv)
    _append_relation_if_needed(query_cmd, _policy_relation(args, task="query-rank", source_type=args.source_type, target_type=target_type))
    query_cmd.extend(_query_options(args, task="query-rank", source_type=args.source_type, target_type=target_type))
    plan["blocks"].append(_block("query", "query", query_cmd, [source], [query_csv], True, notes=_query_notes(args, task="query-rank", source_type=args.source_type, target_type=target_type)))
    plan["blocks"].append(_block("normalize_evidence", "normalize-evidence", _cmd("reason", "normalize", "-i", query_csv, "--deduplicate", "-o", evidence_csv), [query_csv], [evidence_csv]))
    plan["blocks"].append(_block("rank", "rank", _cmd("rk", "-i", query_csv, "--rank", "-o", ranked_csv), [query_csv], [ranked_csv]))
    plan["blocks"].append(_block("explain", "explain", _cmd("reason", "explain", "-e", evidence_csv, "--hypotheses", ranked_csv, "-o", explanation_md), [evidence_csv, ranked_csv], [explanation_md]))
    return plan


def _query_enrich_plan(args):
    out_dir = _ensure_out_dir(args)
    source = args.input or args.genes
    if not source:
        raise ValueError("--input or --genes is required for query-enrich.")
    target_types = _split_csv_arg(args.target_types) or ["Pathway", "Process", "Function", "Disease", "Cell"]
    query_files = []
    plan = _base_plan(
        "brick query-enrich-plot plan",
        "Query one source set against multiple target types, merge results, enrich, and plot.",
        out_dir,
    )
    for target_type in target_types:
        query_csv = out_dir / f"query_{target_type.lower()}.csv"
        query_files.append(query_csv)
        query_cmd = _cmd("qr", "-i", source, "--st", args.source_type, "--tt", target_type, "-n", args.query_limit, "-o", query_csv)
        _append_relation_if_needed(query_cmd, _policy_relation(args, task="query-enrich", source_type=args.source_type, target_type=target_type))
        query_cmd.extend(_query_options(args, task="query-enrich", source_type=args.source_type, target_type=target_type))
        plan["blocks"].append(_block(f"query_{target_type.lower()}", "query", query_cmd, [source], [query_csv], True, notes=_query_notes(args, task="query-enrich", source_type=args.source_type, target_type=target_type)))

    merged_csv = out_dir / "combined_query.csv"
    evidence_csv = out_dir / "evidence_units.csv"
    enriched_csv = out_dir / "enriched.csv"
    explanation_md = out_dir / "enrichment_explanation.md"
    radial_png = out_dir / "enriched_radial.png"
    plan["blocks"].append(_block("merge_queries", "merge-csv", _cmd("plan", "merge-csv", "-o", merged_csv, *query_files), query_files, [merged_csv]))
    plan["blocks"].append(_block("normalize_evidence", "normalize-evidence", _cmd("reason", "normalize", "-i", merged_csv, "--deduplicate", "-o", evidence_csv), [merged_csv], [evidence_csv]))
    plan["blocks"].append(_block("enrich", "enrich", _cmd("rk", "-i", merged_csv, "--enrich", "-o", enriched_csv), [merged_csv], [enriched_csv]))
    plan["blocks"].append(_block("explain", "explain", _cmd("reason", "explain", "-e", evidence_csv, "--hypotheses", enriched_csv, "-o", explanation_md), [evidence_csv, enriched_csv], [explanation_md]))
    plan["blocks"].append(_block("radial_plot", "plot", _cmd("pl", "radial", "-i", enriched_csv, "--top-k", args.top_k, "-o", radial_png), [enriched_csv], [radial_png]))
    return plan


def _graph_complete_plan(args):
    out_dir = _ensure_out_dir(args)
    if not args.graph or not args.query:
        raise ValueError("--graph and --query are required for graph-complete.")
    merged = out_dir / "merged.graphml"
    info = out_dir / "graph_info.txt"
    static_png = out_dir / "graph.png"
    html = out_dir / "graph.html"
    plan = _base_plan(
        "brick graph-completion plan",
        "Merge an observed graph with KG query results, inspect it, and create plots.",
        out_dir,
    )
    cmd = _cmd("ga", "complete_result", "-g", args.graph, "-q", args.query, "-o", merged)
    if args.node_type_filter:
        cmd.extend(["--node-type-filter", args.node_type_filter])
    plan["blocks"].append(_block("complete_graph", "complete-graph", cmd, [args.graph, args.query], [merged]))
    plan["blocks"].append(_block("graph_info", "graph-info", _cmd("ga", "info", "-i", merged), [merged], [info], notes="Redirect stdout to this file if desired."))
    plot_cmd = _cmd("pl", "static", "-i", merged, "-o", static_png)
    if args.root:
        plot_cmd.extend(["--layout", "tree", "--root", args.root])
    plan["blocks"].append(_block("static_plot", "plot", plot_cmd, [merged], [static_png]))
    plan["blocks"].append(_block("interactive_plot", "plot", _cmd("pl", "interact", "-i", merged, "-o", html), [merged], [html]))
    return plan


def _celltype_plan(args):
    out_dir = _ensure_out_dir(args)
    if not args.marker_dir:
        raise ValueError("--marker-dir is required for celltype plans. Use brick workflow celltype for AnnData marker extraction.")
    marker_dir = Path(args.marker_dir)
    marker_files = sorted(marker_dir.glob("*.txt"))
    if not marker_files:
        raise ValueError(f"No marker .txt files found in {marker_dir}")
    plan = _base_plan(
        "brick cell type annotation block plan",
        "For each marker file, query Gene->Cell marker evidence, rank candidates, and summarize evidence.",
        out_dir,
    )
    for marker_file in marker_files:
        stem = marker_file.stem
        query_csv = out_dir / f"{stem}_query_cell.csv"
        evidence_csv = out_dir / f"{stem}_evidence_units.csv"
        ranked_csv = out_dir / f"{stem}_ranked_cell.csv"
        explanation_md = out_dir / f"{stem}_explanation.md"
        query_cmd = _cmd("qr", "-i", marker_file, "--st", "Gene", "--tt", "Cell", "-n", args.query_limit, "-o", query_csv)
        _append_relation_if_needed(query_cmd, _policy_relation(args, task="celltype", source_type="Gene", target_type="Cell"))
        query_cmd.extend(_query_options(args, task="celltype", source_type="Gene", target_type="Cell"))
        plan["blocks"].append(_block(f"{stem}_query", "query", query_cmd, [marker_file], [query_csv], True, notes=_query_notes(args, task="celltype", source_type="Gene", target_type="Cell")))
        plan["blocks"].append(_block(f"{stem}_normalize", "normalize-evidence", _cmd("reason", "normalize", "-i", query_csv, "--deduplicate", "-o", evidence_csv), [query_csv], [evidence_csv]))
        plan["blocks"].append(_block(f"{stem}_rank", "rank", _cmd("rk", "-i", query_csv, "--rank", "--metrics", "match_count,info_source_count", "--prune-parent", "-o", ranked_csv), [query_csv], [ranked_csv]))
        plan["blocks"].append(_block(f"{stem}_explain", "explain", _cmd("reason", "explain", "-e", evidence_csv, "--hypotheses", ranked_csv, "-o", explanation_md), [evidence_csv, ranked_csv], [explanation_md]))
    return plan


def _drug_plan(args):
    out_dir = _ensure_out_dir(args)
    genes = args.genes or args.input
    if not genes:
        raise ValueError("--genes or --input is required for drug plans.")
    if (args.reference or args.disease) and not args.evaluation_only:
        raise ValueError("Drug references require --evaluation-only and are appended after ranking freeze.")
    if args.method == "grace" and not args.adata:
        raise ValueError("--method grace requires --adata.")
    plan = _base_plan(
        "brick drug repositioning block plan",
        "Build a stable-ID-first ranking, freeze it, then optionally evaluate against a post-freeze reference.",
        out_dir,
    )
    target_types = (
        [value.strip() for value in (args.target_types or DEFAULT_DRUG_GRACE_TARGET_TYPES).split(",") if value.strip()]
        if args.method == "grace"
        else ["Chemical", "Disease"]
    )
    query_files = []
    for target_type in target_types:
        query_csv = out_dir / f"query_{target_type.lower()}.csv"
        query_files.append(query_csv)
        query_cmd = _cmd("qr", "-i", genes, "--st", "Gene", "--tt", target_type, "-n", args.query_limit, "-o", query_csv)
        if args.method == "rank":
            relation = "causally_related_to,disease_relation,participate_in" if target_type == "Disease" else None
            task = "disease" if target_type == "Disease" else "drug"
            _append_relation_if_needed(query_cmd, _policy_relation(args, task=task, source_type="Gene", target_type=target_type, relation=relation))
            query_cmd.extend(_query_options(args, task=task, source_type="Gene", target_type=target_type, relation=relation))
        plan["blocks"].append(_block(f"query_{target_type.lower()}", "query", query_cmd, [genes], [query_csv], True))

    q_chem = out_dir / "query_chemical.csv"
    q_dis = out_dir / "query_disease.csv"
    merged = out_dir / "combined_query.csv"
    evidence_csv = out_dir / "evidence_units.csv"
    primary_ranking = out_dir / "primary_chemical_ranking.csv"
    predicted = out_dir / "predicted_drugs.csv"
    scope_csv = out_dir / "chemical_scope.csv"
    freeze = out_dir / "ranking_frozen.sha256"
    explanation = out_dir / "drug_explanation.md"
    plan["blocks"].append(_block("merge_queries", "merge-csv", _cmd("plan", "merge-csv", "-o", merged, *query_files), query_files, [merged]))
    plan["blocks"].append(_block("normalize_evidence", "normalize-evidence", _cmd("reason", "normalize", "-i", merged, "--deduplicate", "-o", evidence_csv), [merged], [evidence_csv]))
    if args.method == "rank":
        ranked_chemicals = out_dir / "ranked_chemicals.csv"
        ranked_diseases = out_dir / "ranked_diseases.csv"
        plan["blocks"].append(_block("rank_chemicals", "rank", _cmd("rk", "-i", q_chem, "--rank", "-o", ranked_chemicals), [q_chem], [ranked_chemicals]))
        plan["blocks"].append(_block("rank_diseases_context", "rank", _cmd("rk", "-i", q_dis, "--rank", "-o", ranked_diseases), [q_dis], [ranked_diseases]))
        plan["blocks"].append(_block("scope_chemicals", "reason", _cmd("drug", "scope", "-i", q_chem, "-o", scope_csv), [q_chem], [scope_csv], notes="Diagnostic tag only; does not filter the formal candidate set."))
        plan["blocks"].append(_block("extract_primary_ranking", "reason", _cmd("drug", "extract", "-i", ranked_chemicals, "--scope", scope_csv, "--top-n", 0, "-o", primary_ranking), [ranked_chemicals, scope_csv], [primary_ranking]))
        plan["blocks"].append(_block("extract_top_candidates", "reason", _cmd("drug", "extract", "-i", ranked_chemicals, "--scope", scope_csv, "--top-n", args.top_n_drugs, "-o", predicted), [ranked_chemicals, scope_csv], [predicted]))
        plan["blocks"].append(_block("freeze_primary_ranking", "reason", _cmd("drug", "freeze", "-i", primary_ranking, "-o", freeze), [primary_ranking], [freeze]))
    else:
        enrichment_all = out_dir / "native_enrichment_all.csv"
        enrichment_strict = out_dir / "native_enrichment_p1e-4.csv"
        graph = out_dir / "primary_joint_graph.graphml"
        graph_nodes = out_dir / "primary_joint_graph_nodes.csv"
        graph_audit = out_dir / "primary_joint_graph_audit.json"
        plan["blocks"].append(_block("enrich_kg_wide", "enrich", _cmd("rk", "-i", merged, "--enrich", "--background-method", "auto", "--pvalue-cutoff", args.pvalue_cutoff, "--relation-batch-size", 10, "--relation-workers", 4, "--relation-retries", 2, "-o", enrichment_all), [merged], [enrichment_all]))
        plan["blocks"].append(_block("filter_strict_enrichment", "reason", _cmd("drug", "filter-enrichment", "-i", enrichment_all, "--pvalue-cutoff", args.pvalue_cutoff, "-o", enrichment_strict), [enrichment_all], [enrichment_strict]))
        plan["blocks"].append(_block("scope_chemicals", "reason", _cmd("drug", "scope", "-i", merged, "-o", scope_csv), [merged], [scope_csv], notes="Diagnostic tag only; does not filter the formal candidate set."))
        plan["blocks"].append(_block("build_joint_graph", "graph", _cmd("drug", "build-graph", "-a", args.adata, "--genes", genes, "-i", enrichment_strict, "--disease", args.anchor_disease, "--pvalue-cutoff", args.pvalue_cutoff, "--identity-output", graph_nodes, "--audit-output", graph_audit, "-o", graph), [args.adata, genes, enrichment_strict], [graph, graph_nodes, graph_audit]))
        grace_cmd = ["brick-gpu", "drug", "grace-rank", "-g", str(graph), "--identity", str(graph_nodes), "--disease", str(args.anchor_disease), "--seeds", str(args.grace_seeds), "--epochs", str(args.grace_epochs), "--device", "cuda", "-o", str(out_dir)]
        plan["blocks"].append(_block("grace_consensus_and_freeze", "embedding", grace_cmd, [graph, graph_nodes], [primary_ranking, freeze]))
        plan["blocks"].append(_block("extract_top_candidates", "reason", _cmd("drug", "extract", "-i", primary_ranking, "--scope", scope_csv, "--top-n", args.top_n_drugs, "-o", predicted), [primary_ranking, scope_csv], [predicted]))
    plan["blocks"].append(_block("explain", "explain", _cmd("reason", "explain", "-e", evidence_csv, "--hypotheses", primary_ranking, "-o", explanation), [evidence_csv, primary_ranking], [explanation]))

    reference_path = Path(args.reference) if args.reference else out_dir / "known_drugs_for_disease.csv"
    if args.disease:
        query_known_cmd = _cmd("qr", args.disease, "--st", "Disease", "--tt", "Chemical", "-r", "pharmaco_effect_on", "-n", args.query_limit, "-o", reference_path)
        plan["blocks"].append(_block("query_known_drugs_post_freeze", "query", query_known_cmd, [freeze], [reference_path], True, notes="Evaluation only; must not alter the frozen ranking."))
    if args.reference or args.disease:
        validation = out_dir / "ranking_validation.csv"
        identity_audit = out_dir / "chemical_identity_audit.csv"
        plan["blocks"].append(_block("validate_frozen_ranking", "evaluation", _cmd("drug", "validate", "--ranking", primary_ranking, "--freeze", freeze, "--reference", reference_path, "--identity-audit", identity_audit, "-o", validation), [primary_ranking, freeze, reference_path], [validation, identity_audit]))
    return plan


def _trajectory_plan(args):
    if not args.graph or not args.query:
        raise ValueError("--graph and --query are required for trajectory plans.")
    args.node_type_filter = args.node_type_filter or "Cell"
    return _graph_complete_plan(args)


def build_template(args):
    builders = {
        "query-rank": _query_rank_plan,
        "query-enrich": _query_enrich_plan,
        "graph-complete": _graph_complete_plan,
        "celltype": _celltype_plan,
        "enrichment": _query_enrich_plan,
        "drug": _drug_plan,
        "trajectory": _trajectory_plan,
    }
    plan = builders[args.task](args)
    return plan


def detect_artifact(path):
    path = Path(path)
    result = {
        "path": str(path),
        "exists": path.exists(),
        "kind": "missing",
        "rows": None,
        "columns": [],
        "notes": [],
    }
    if not path.exists():
        return result

    suffix = path.suffix.lower()
    result["size_bytes"] = path.stat().st_size
    if suffix == ".csv":
        try:
            df = pd.read_csv(path, nrows=200)
        except pd.errors.EmptyDataError:
            result["kind"] = "empty_csv"
            return result
        result["rows"] = _count_csv_rows(path)
        result["columns"] = list(df.columns)
        cols = set(df.columns)
        if RANK_COLUMNS.issubset(cols):
            result["kind"] = "brick_ranked_csv"
        elif ENRICH_COLUMNS.issubset(cols):
            result["kind"] = "brick_enrichment_csv"
        elif QUERY_COLUMNS.issubset(cols) and cols.intersection(RELATION_COLUMNS):
            result["kind"] = "brick_query_csv"
        else:
            result["kind"] = "csv"
        if "path.1.original_relation" in cols:
            result["notes"].append("contains evidence snippets")
        return result

    if suffix in {".txt", ".tsv"}:
        text = path.read_text(encoding="utf-8", errors="ignore")
        non_empty = [line for line in text.splitlines() if line.strip()]
        result["kind"] = "entity_list" if suffix == ".txt" else "tsv"
        result["rows"] = len(non_empty)
        return result

    if suffix == ".h5ad":
        result["kind"] = "anndata"
        return result

    if suffix in {".graphml", ".gml", ".gpickle", ".pkl", ".json"}:
        result["kind"] = "graph"
        _fill_graph_counts(path, result)
        return result

    result["kind"] = suffix.lstrip(".") or "file"
    return result


def _count_csv_rows(path):
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            return max(sum(1 for _ in csv.reader(handle)) - 1, 0)
    except UnicodeDecodeError:
        with open(path, "r", encoding="latin-1", newline="") as handle:
            return max(sum(1 for _ in csv.reader(handle)) - 1, 0)


def _fill_graph_counts(path, result):
    try:
        import networkx as nx
        suffix = path.suffix.lower()
        if suffix == ".graphml":
            graph = nx.read_graphml(path)
        elif suffix == ".gml":
            graph = nx.read_gml(path)
        else:
            result["notes"].append("graph counts skipped for this format")
            return
        result["nodes"] = graph.number_of_nodes()
        result["edges"] = graph.number_of_edges()
    except Exception as exc:
        result["notes"].append(f"graph inspection failed: {exc}")


def cmd_blocks(args):
    if args.format == "json":
        _write_json(BLOCKS)
        return
    rows = [
        {
            "name": item["name"],
            "category": item["category"],
            "api": "yes" if item["requires_api"] else "no",
            "command": item["command"],
            "purpose": item["purpose"],
        }
        for item in BLOCKS
    ]
    _print_table(rows, ["name", "category", "api", "command", "purpose"])


def cmd_schema(args):
    data = {"entity_types": ENTITY_TYPES, "relations": RELATIONS}
    if args.format == "json":
        _write_json(data)
        return
    print("Entity types:")
    print(", ".join(ENTITY_TYPES))
    print("")
    _print_table(RELATIONS, ["name", "typical_source", "typical_target", "use"])


def cmd_inspect(args):
    results = [detect_artifact(path) for path in args.inputs]
    if args.format == "json":
        _write_json(results, args.output)
        return
    rows = []
    for item in results:
        rows.append({
            "path": item["path"],
            "kind": item["kind"],
            "rows": item.get("rows", ""),
            "nodes": item.get("nodes", ""),
            "edges": item.get("edges", ""),
            "notes": "; ".join(item.get("notes", [])),
        })
    _print_table(rows, ["path", "kind", "rows", "nodes", "edges", "notes"])


def cmd_template(args):
    plan = build_template(args)
    if getattr(args, "question", None):
        plan["question"] = args.question
    if getattr(args, "agent_plan", False):
        plan = to_agent_plan(plan)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    _write_json(plan, args.output)


def load_plan(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def validate_plan(plan):
    issues = []
    for block in plan.get("blocks", []):
        command = block.get("command") or []
        if not command:
            issues.append({"block": block.get("id"), "level": "error", "message": "missing command"})
            continue
        executable = str(command[0])
        if shutil.which(executable) is None:
            issues.append({"block": block.get("id"), "level": "warning", "message": f"command not found on PATH: {executable}"})
        for input_path in block.get("inputs", []):
            if input_path and "{" not in input_path and not Path(input_path).exists():
                issues.append({"block": block.get("id"), "level": "warning", "message": f"input does not exist yet: {input_path}"})
        for output_path in block.get("outputs", []):
            parent = Path(output_path).parent
            if str(parent) not in ("", ".") and not parent.exists():
                issues.append({"block": block.get("id"), "level": "info", "message": f"output directory will be created: {parent}"})
    return issues


def cmd_validate(args):
    plan = load_plan(args.plan)
    issues = validate_plan(plan)
    if args.format == "json":
        _write_json({"valid": not any(x["level"] == "error" for x in issues), "issues": issues}, args.output)
        return
    if not issues:
        print("Plan is valid.")
        return
    rows = [{"level": x["level"], "block": x["block"], "message": x["message"]} for x in issues]
    _print_table(rows, ["level", "block", "message"])


def cmd_materialize(args):
    plan = load_plan(args.plan)
    out_dir = Path(args.output_dir or plan.get("output_dir") or ".")
    out_dir.mkdir(parents=True, exist_ok=True)
    is_agent_plan = "steps" in plan and "blocks" not in plan
    plan_copy = out_dir / ("agent_plan.json" if is_agent_plan else "brick_plan.json")
    plan_copy.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if is_agent_plan:
        commands = [step["command"] for step in plan.get("steps", []) if step.get("command")]
    else:
        commands = [block["command"] for block in plan.get("blocks", []) if block.get("command")]
    ps1 = out_dir / "commands.ps1"
    sh = out_dir / "commands.sh"
    ps1.write_text("$ErrorActionPreference = 'Stop'\n\n" + "\n".join(_command_to_ps(command) for command in commands) + "\n", encoding="utf-8")
    sh.write_text("#!/usr/bin/env bash\nset -euo pipefail\n\n" + "\n".join(_command_to_sh(command) for command in commands) + "\n", encoding="utf-8")
    print(f"Wrote {plan_copy}")
    print(f"Wrote {ps1}")
    print(f"Wrote {sh}")


def cmd_merge_csv(args):
    frames = []
    for input_path in args.inputs:
        path = Path(input_path)
        if not path.exists() or path.stat().st_size == 0:
            continue
        try:
            df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            continue
        if not df.empty:
            frames.append(df)
    if not frames:
        raise ValueError("No non-empty CSV inputs to merge.")
    merged = pd.concat(frames, ignore_index=True, sort=False)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, index=False)
    print(f"Wrote {args.output} ({len(merged)} rows)")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="brick plan",
        description="Describe, inspect, and compose brick graph-knowledge blocks.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("blocks", help="List available brick block contracts.")
    p.add_argument("--format", choices=["table", "json"], default="table")
    p.set_defaults(func=cmd_blocks)

    p = sub.add_parser("schema", help="List common entity types and relations.")
    p.add_argument("--format", choices=["table", "json"], default="table")
    p.set_defaults(func=cmd_schema)

    p = sub.add_parser("inspect", help="Inspect brick artifacts and infer their type.")
    p.add_argument("inputs", nargs="+")
    p.add_argument("--format", choices=["table", "json"], default="table")
    p.add_argument("-o", "--output")
    p.set_defaults(func=cmd_inspect)

    p = sub.add_parser("template", help="Create a JSON block plan for a common composition.")
    p.add_argument("--task", required=True, choices=["query-rank", "query-enrich", "graph-complete", "celltype", "enrichment", "drug", "trajectory"])
    p.add_argument("-i", "--input")
    p.add_argument("--genes")
    p.add_argument("--adata")
    p.add_argument("--marker-dir")
    p.add_argument("--graph")
    p.add_argument("--query")
    p.add_argument("--reference")
    p.add_argument("--disease")
    p.add_argument("--source-type", default="Gene")
    p.add_argument("--target-type")
    p.add_argument("--target-types", help="Comma-separated target types for query-enrich/enrichment or GRACE drug discovery.")
    p.add_argument("--method", choices=["rank", "grace"], default="rank", help="Drug task method.")
    p.add_argument("--anchor-disease", default="COVID-19")
    p.add_argument("--pvalue-cutoff", default="0.0001")
    p.add_argument("--grace-seeds", default="20241204,20241205,20241206")
    p.add_argument("--grace-epochs", default="1000")
    p.add_argument("--top-n-drugs", default="50")
    p.add_argument("--evaluation-only", action="store_true")
    p.add_argument("-r", "--relation")
    p.add_argument("--hops", "--multi-hop", dest="hops", type=int, choices=[1, 2, 3], default=1, help="Optional brick qr path hops for query blocks: 1 direct only, 2/3 multi-hop.")
    p.add_argument("--directed", action="store_true", help="Add --directed to generated brick qr query blocks.")
    p.add_argument("--question", help="Natural-language question/context for automatic KG source routing.")
    p.add_argument("--kg-policy", choices=KG_POLICY_CHOICES, default="auto", help="KG query policy for generated query blocks.")
    p.add_argument("--info-source", help="Override edge info_source filter for generated query blocks.")
    p.add_argument("--policy-file", help="Optional JSON KG source policy catalog.")
    p.add_argument("--query-limit", default="500")
    p.add_argument("--top-k", default="15")
    p.add_argument("--root")
    p.add_argument("--node-type-filter")
    p.add_argument("--output-dir", default="brick_plan_output")
    p.add_argument("--agent-plan", action="store_true", help="Emit a deterministic brick agent plan instead of the legacy block plan.")
    p.add_argument("-o", "--output")
    p.set_defaults(func=cmd_template)

    p = sub.add_parser("validate", help="Validate a brick JSON plan.")
    p.add_argument("plan")
    p.add_argument("--format", choices=["table", "json"], default="table")
    p.add_argument("-o", "--output")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("materialize", help="Write commands.ps1 and commands.sh from a brick JSON plan.")
    p.add_argument("plan")
    p.add_argument("-o", "--output-dir")
    p.set_defaults(func=cmd_materialize)

    p = sub.add_parser("merge-csv", help="Merge non-empty CSV files.")
    p.add_argument("inputs", nargs="+")
    p.add_argument("-o", "--output", required=True)
    p.set_defaults(func=cmd_merge_csv)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


def _cli_main():
    return main()


if __name__ == "__main__":
    sys.exit(main())
