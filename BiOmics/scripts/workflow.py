"""brick workflow - curated brick workflows for common BiOmics tasks."""

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

try:
    from kgpolicy import KG_POLICY_CHOICES, select_kg_policy
except Exception:  # pragma: no cover
    KG_POLICY_CHOICES = ("manual", "auto", "off", "direct", "source-scoped", "multihop")
    select_kg_policy = None


DEFAULT_TARGET_TYPES = [
    "Cell",
    "Tissue",
    "Disease",
    "Phenotype",
    "Chemical",
    "Process",
    "Function",
    "Cell_Component",
    "Pathway",
]

DEFAULT_ENRICH_RELATIONS = (
    "produce,participate_in,causally_related_to,capable_of,"
    "spatially_related_to,processly_related_to,functionally_related_to,"
    "disease_relation,pharmaco_effect_on"
)

DEFAULT_DRUG_GRACE_TARGET_TYPES = (
    "Disease,Phenotype,Process,Cell_Component,Pathway,Function,"
    "Tissue,Cell,Chemical,Stage"
)


def _ensure_dir(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write_text(path, text):
    Path(path).write_text(text, encoding="utf-8")


def _run(cmd, cwd=None):
    print("+ " + " ".join(str(x) for x in cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def _resolve_policy(args, task, source_type, target_type, relation=None):
    if select_kg_policy is None:
        return {}
    return select_kg_policy(
        question=getattr(args, "question", None),
        source_type=source_type,
        target_type=target_type,
        relation=relation,
        info_source=getattr(args, "info_source", None),
        hops=getattr(args, "hops", 1),
        task=task,
        kg_policy=getattr(args, "kg_policy", "auto"),
        policy_file=getattr(args, "policy_file", None),
    )


def _qr_command(args, source_or_file, output, source_type, target_type, relation=None, task=None, positional=False):
    """Build a policy-scoped brick qr command list for workflow scripts."""
    policy = _resolve_policy(args, task or "workflow", source_type, target_type, relation)
    cmd = ["brick", "qr"]
    if positional:
        cmd.append(str(source_or_file))
    else:
        cmd.extend(["-i", str(source_or_file)])
    cmd.extend(["--st", source_type, "--tt", target_type])
    rel = relation or ",".join(policy.get("relations") or [])
    if rel:
        cmd.extend(["-r", rel])
    info_source = getattr(args, "info_source", None)
    if info_source:
        cmd.extend(["--info-source", info_source])
    elif policy.get("primary_info_sources"):
        cmd.extend(["--info-source", ",".join(policy["primary_info_sources"])])
    hops = int(policy.get("hops") or getattr(args, "hops", 1) or 1)
    if hops > 1:
        cmd.extend(["--hops", str(hops)])
    if getattr(args, "kg_policy", "auto") != "auto":
        cmd.extend(["--kg-policy", getattr(args, "kg_policy")])
    if getattr(args, "question", None):
        cmd.extend(["--question", getattr(args, "question")])
    if getattr(args, "policy_file", None):
        cmd.extend(["--policy-file", getattr(args, "policy_file")])
    cmd.extend(["-n", str(args.query_limit), "-o", str(output)])
    return [str(x) for x in cmd]


def _cmd_to_sh(cmd):
    return " ".join(shlex.quote(str(x)) for x in cmd)


def _reasoning_commands(query_csv, evidence_csv, explanation_md=None, hypotheses_csv=None):
    commands = [["brick", "reason", "normalize", "-i", str(query_csv), "--deduplicate", "-o", str(evidence_csv)]]
    if explanation_md:
        cmd = ["brick", "reason", "explain", "-e", str(evidence_csv), "-o", str(explanation_md)]
        if hypotheses_csv:
            cmd.extend(["--hypotheses", str(hypotheses_csv)])
        commands.append(cmd)
    return commands


def _copy_or_write_genes(args, out_dir, default_name="geneset.txt"):
    genes_path = out_dir / default_name
    if getattr(args, "genes", None):
        shutil.copyfile(args.genes, genes_path)
        return genes_path

    if not getattr(args, "adata", None):
        raise ValueError("Provide either --genes or --adata.")

    import scanpy as sc

    adata = sc.read(args.adata)
    if getattr(args, "layer", None):
        adata.X = adata.layers[args.layer].copy()
    groupby = args.groupby
    group = getattr(args, "group", None)
    top_n = getattr(args, "top_n_genes", 100)

    if groupby not in adata.obs:
        raise ValueError(f"groupby column not found in adata.obs: {groupby}")

    if "rank_genes_groups" not in adata.uns:
        sc.tl.rank_genes_groups(adata, groupby=groupby, method=args.de_method, pts=True)

    names = adata.uns["rank_genes_groups"]["names"]
    groups = [group] if group else list(names.dtype.names or [])
    if not groups:
        raise ValueError("No rank_genes_groups entries found.")

    genes = []
    for group_name in groups:
        if group_name not in names.dtype.names:
            raise ValueError(f"group not found in rank_genes_groups: {group_name}")
        genes.extend([str(x) for x in list(names[group_name])[:top_n]])

    unique = sorted(set(g for g in genes if g and g != "nan"))
    _write_text(genes_path, "\n".join(unique) + "\n")
    return genes_path


def _merge_csvs(paths, output):
    frames = []
    for path in paths:
        path = Path(path)
        if path.exists() and path.stat().st_size > 0:
            df = pd.read_csv(path)
            if not df.empty:
                frames.append(df)
    if not frames:
        return None
    merged = pd.concat(frames, ignore_index=True)
    merged.to_csv(output, index=False)
    return output


def _commands_files(out_dir, commands):
    sh = out_dir / "commands.sh"
    ps1 = out_dir / "commands.ps1"
    _write_text(sh, "#!/usr/bin/env bash\nset -euo pipefail\n\n" + "\n".join(commands) + "\n")
    _write_text(ps1, "$ErrorActionPreference = 'Stop'\n\n" + "\n".join(commands) + "\n")
    return sh, ps1


def _extract_output_from_command(command_text):
    matches = re.findall(r"(?:^|\s)-o\s+(?:\"([^\"]+)\"|'([^']+)'|([^\s]+))", command_text)
    if not matches:
        return []
    outputs = []
    for quoted, single, bare in matches:
        outputs.append(quoted or single or bare)
    return outputs


def _agent_type_and_validators(command_text):
    text = command_text.lower()
    common = ["file_exists", "non_empty_file"]
    if text.startswith("brick qr"):
        return "kg_query", [
            "kg_policy_not_broad",
            "max_hops_respected",
            *common,
            "kg_edge_metadata_columns_present",
            "kg_edge_metadata_required_for_co_mentioned",
        ]
    if text.startswith("brick reason normalize"):
        return "evidence_normalization", [
            *common,
            "evidence_strength_present",
            "manual_review_notes_present",
            "co_mentioned_weak_cap",
            "max_hops_respected",
        ]
    if text.startswith("brick rk") and "--enrich" in text:
        return "enrichment", common
    if text.startswith("brick rk"):
        return "ranking", common
    if text.startswith("brick ga"):
        return "graph", common
    if text.startswith("brick plan"):
        return "utility", common
    if text.startswith("brick pl"):
        return "visualization", common
    if text.startswith("brick-gpu drug grace-rank") or text.startswith("brick drug grace-rank"):
        return "embedding", [*common, "cuda_required"]
    if text.startswith("brick drug validate"):
        return "evaluation", [*common, "stable_id_primary_identity", "ranking_frozen_before_validation"]
    if text.startswith("brick drug"):
        return "drug_analysis", common
    if text.startswith("brick reason explain"):
        return "report", [*common, "no_weak_only_strong_claim"]
    if text.startswith("brick reason"):
        return "reasoning", common
    return "utility", common


def _write_agent_plan(out_dir, task, commands, question=None):
    steps = []
    previous_step = None
    for index, command_text in enumerate(commands, start=1):
        step_type, validators = _agent_type_and_validators(command_text)
        step_id = f"step_{index:03d}_{step_type}"
        steps.append(
            {
                "id": step_id,
                "type": step_type,
                "command": command_text,
                "inputs": [],
                "outputs": _extract_output_from_command(command_text),
                "depends_on": [previous_step] if previous_step else [],
                "validators": validators,
                "validator_args": {"max_hops_respected": {"max_hops": 3}},
                "retry_policy": {
                    "max_retries": 1 if step_type == "kg_query" else 0,
                    "retry_on": ["command_failed", "validator_failed"],
                },
                "checkpoint": step_type in {"kg_query", "reasoning", "report"},
                "allow_failure": False,
            }
        )
        previous_step = step_id
    plan = {
        "schema_version": "1.0",
        "task": task,
        "question": question or "",
        "workdir": str(out_dir),
        "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "execution": {
            "mode": "dry_run",
            "executor": "brick",
            "allow_shell": False,
            "max_retries": 1,
            "timeout_seconds": 600,
        },
        "steps": steps,
    }
    output = out_dir / "agent_plan.json"
    output.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return output


def _extract_markers(adata_path, groupby, out_dir, top_n=10, layer=None, de_method="wilcoxon"):
    import scanpy as sc

    adata = sc.read(adata_path)
    if layer:
        adata.X = adata.layers[layer].copy()
    if groupby not in adata.obs:
        raise ValueError(f"groupby column not found in adata.obs: {groupby}")
    if "rank_genes_groups" not in adata.uns:
        sc.tl.rank_genes_groups(adata, groupby=groupby, method=de_method, pts=True)

    names = adata.uns["rank_genes_groups"]["names"]
    groups = list(names.dtype.names or [])
    marker_dir = _ensure_dir(out_dir / "markers")
    marker_files = []
    for group in groups:
        genes = [str(x) for x in list(names[group])[:top_n] if str(x) != "nan"]
        path = marker_dir / f"cluster_{group}_markers.txt"
        _write_text(path, "\n".join(genes) + "\n")
        marker_files.append((group, path))

    adata.write(out_dir / "adata_with_rank_genes.h5ad")
    return marker_files


def cmd_celltype(args):
    out_dir = _ensure_dir(args.output_dir)
    marker_files = []
    if args.adata:
        marker_files = _extract_markers(
            args.adata,
            args.groupby,
            out_dir,
            top_n=args.top_n_markers,
            layer=args.layer,
            de_method=args.de_method,
        )

    if args.marker_dir:
        marker_dir = Path(args.marker_dir)
        marker_files = [(p.stem, p) for p in sorted(marker_dir.glob("*.txt"))]

    if not marker_files:
        raise ValueError("Provide --adata or --marker-dir for celltype workflow.")

    commands = []
    for group, marker_file in marker_files:
        query_csv = out_dir / f"cluster_{group}_query.csv"
        ranked_csv = out_dir / f"cluster_{group}_ranked.csv"
        evidence_csv = out_dir / f"cluster_{group}_evidence_units.csv"
        explanation_md = out_dir / f"cluster_{group}_explanation.md"
        qr_cmd = _qr_command(args, marker_file, query_csv, "Gene", "Cell", task="celltype")
        commands.append(_cmd_to_sh(qr_cmd))
        for reason_cmd in _reasoning_commands(query_csv, evidence_csv):
            commands.append(_cmd_to_sh(reason_cmd))
        rk_cmd = ["brick", "rk", "-i", str(query_csv), "--rank", "--metrics", args.rank_metrics]
        if args.prune_parent:
            rk_cmd.append("--prune-parent")
        rk_cmd.extend(["-o", str(ranked_csv)])
        commands.append(_cmd_to_sh(rk_cmd))
        commands.append(_cmd_to_sh(["brick", "reason", "explain", "-e", str(evidence_csv), "--hypotheses", str(ranked_csv), "-o", str(explanation_md)]))

        if args.run_query:
            _run(qr_cmd)
            for reason_cmd in _reasoning_commands(query_csv, evidence_csv):
                _run(reason_cmd)
            _run(rk_cmd)
            _run(["brick", "reason", "explain", "-e", str(evidence_csv), "--hypotheses", str(ranked_csv), "-o", str(explanation_md)])

    _commands_files(out_dir, commands)
    if getattr(args, "agent_plan", False):
        agent_plan = _write_agent_plan(out_dir, "celltype", commands, question=getattr(args, "question", None))
        print(f"Wrote {agent_plan}")
    _summarize_celltype(out_dir)
    print(f"Cell type workflow prepared in {out_dir}")


def _summarize_celltype(out_dir):
    rows = []
    for ranked in sorted(Path(out_dir).glob("cluster_*_ranked.csv")):
        df = pd.read_csv(ranked)
        if df.empty:
            continue
        top = df.sort_values("path.2.rank_voting").iloc[0]
        cluster = ranked.name.replace("cluster_", "").replace("_ranked.csv", "")
        rows.append({
            "cluster": cluster,
            "celltype": top.get("path.2.name"),
            "rank_voting": top.get("path.2.rank_voting"),
            "match_count": top.get("path.2.match_count"),
        })
    if rows:
        pd.DataFrame(rows).to_csv(Path(out_dir) / "celltype_annotations.csv", index=False)


def cmd_enrichment(args):
    out_dir = _ensure_dir(args.output_dir)
    genes_path = _copy_or_write_genes(args, out_dir)
    target_types = [x.strip() for x in args.target_types.split(",") if x.strip()]
    commands = []
    query_paths = []
    for target_type in target_types:
        query_csv = out_dir / f"query_{target_type.lower()}.csv"
        query_paths.append(query_csv)
        qr_cmd = _qr_command(args, genes_path, query_csv, "Gene", target_type, relation=args.relations or None, task="enrichment")
        commands.append(_cmd_to_sh(qr_cmd))
        if args.run_query:
            _run(qr_cmd)

    combined = _merge_csvs(query_paths, out_dir / "combined_query.csv")
    evidence_csv = out_dir / "evidence_units.csv"
    explanation_md = out_dir / "enrichment_explanation.md"
    commands.append(_cmd_to_sh(["brick", "plan", "merge-csv", "-o", str(out_dir / "combined_query.csv"), *[str(path) for path in query_paths]]))
    commands.append(_cmd_to_sh(["brick", "reason", "normalize", "-i", str(out_dir / "combined_query.csv"), "--deduplicate", "-o", str(evidence_csv)]))
    commands.append(f'brick rk -i "{out_dir / "combined_query.csv"}" --enrich -o "{out_dir / "enrichment_results.csv"}"')
    commands.append(_cmd_to_sh(["brick", "reason", "explain", "-e", str(evidence_csv), "--hypotheses", str(out_dir / "enrichment_results.csv"), "-o", str(explanation_md)]))
    commands.append(f'brick pl radial -i "{out_dir / "enrichment_results.csv"}" -o "{out_dir / "radial_enrichment.png"}" --top-k {args.top_k}')

    if args.run_analysis and combined:
        _run(["brick", "reason", "normalize", "-i", str(combined), "--deduplicate", "-o", str(evidence_csv)])
        _run(["brick", "rk", "-i", str(combined), "--enrich", "-o", str(out_dir / "enrichment_results.csv")])
        _run(["brick", "reason", "explain", "-e", str(evidence_csv), "--hypotheses", str(out_dir / "enrichment_results.csv"), "-o", str(explanation_md)])
        _run([
            "brick",
            "pl",
            "radial",
            "-i",
            str(out_dir / "enrichment_results.csv"),
            "-o",
            str(out_dir / "radial_enrichment.png"),
            "--top-k",
            str(args.top_k),
        ])

    _commands_files(out_dir, commands)
    if getattr(args, "agent_plan", False):
        agent_plan = _write_agent_plan(out_dir, "enrichment", commands, question=getattr(args, "question", None))
        print(f"Wrote {agent_plan}")
    print(f"Enrichment workflow prepared in {out_dir}")


def cmd_drug(args):
    """Prepare or run leakage-resistant rank or GRACE drug discovery."""
    if (args.disease or args.real_drugs or args.run_validation) and not args.evaluation_only:
        raise ValueError(
            "Known-drug/reference access requires --evaluation-only and is scheduled only "
            "after ranking_frozen.sha256 exists."
        )
    if args.run_validation and not (args.disease or args.real_drugs):
        raise ValueError("--run-validation requires --disease or --real-drugs")
    if args.method == "grace" and not args.adata:
        raise ValueError("--method grace requires --adata for the native joint graph")

    out_dir = _ensure_dir(args.output_dir)
    genes_path = _copy_or_write_genes(args, out_dir)
    combined = out_dir / "combined_query.csv"
    evidence_csv = out_dir / "evidence_units.csv"
    explanation_md = out_dir / "drug_explanation.md"
    scope_csv = out_dir / "chemical_scope.csv"
    primary_ranking = out_dir / "primary_chemical_ranking.csv"
    predicted = out_dir / "predicted_drugs.csv"
    freeze_path = out_dir / "ranking_frozen.sha256"

    if args.method == "grace":
        target_types = [value.strip() for value in args.target_types.split(",") if value.strip()]
    else:
        target_types = ["Chemical", "Disease"]
    if "Chemical" not in target_types:
        raise ValueError("Drug discovery target types must include Chemical")
    if args.method == "grace" and "Disease" not in target_types:
        raise ValueError("GRACE drug discovery target types must include Disease")

    query_paths = []
    query_commands = []
    for target_type in target_types:
        query_path = out_dir / f"query_{target_type.lower()}.csv"
        query_paths.append(query_path)
        if args.method == "grace":
            # The notebook-parity graph uses the complete one-hop target union;
            # no relation/source label is selected from validation outcomes.
            command = [
                "brick", "qr", "-i", str(genes_path), "--st", "Gene", "--tt", target_type,
                "-n", str(args.query_limit), "-o", str(query_path),
            ]
        elif target_type == "Chemical":
            command = _qr_command(args, genes_path, query_path, "Gene", "Chemical", task="drug")
        else:
            command = _qr_command(
                args,
                genes_path,
                query_path,
                "Gene",
                "Disease",
                relation="causally_related_to,disease_relation,participate_in",
                task="disease",
            )
        query_commands.append(command)

    merge_cmd = ["brick", "plan", "merge-csv", "-o", str(combined), *[str(path) for path in query_paths]]
    normalize_cmd = ["brick", "reason", "normalize", "-i", str(combined), "--deduplicate", "-o", str(evidence_csv)]
    analysis_commands = [normalize_cmd]

    if args.method == "rank":
        query_chemical = out_dir / "query_chemical.csv"
        query_disease = out_dir / "query_disease.csv"
        ranked_chemicals = out_dir / "ranked_chemicals.csv"
        ranked_diseases = out_dir / "ranked_diseases.csv"
        analysis_commands.extend(
            [
                ["brick", "rk", "-i", str(query_chemical), "--rank", "-o", str(ranked_chemicals)],
                ["brick", "rk", "-i", str(query_disease), "--rank", "-o", str(ranked_diseases)],
                ["brick", "drug", "scope", "-i", str(query_chemical), "-o", str(scope_csv)],
                ["brick", "drug", "extract", "-i", str(ranked_chemicals), "--scope", str(scope_csv), "--top-n", "0", "-o", str(primary_ranking)],
                ["brick", "drug", "extract", "-i", str(ranked_chemicals), "--scope", str(scope_csv), "--top-n", str(args.top_n_drugs), "-o", str(predicted)],
                ["brick", "drug", "freeze", "-i", str(primary_ranking), "-o", str(freeze_path)],
                ["brick", "reason", "explain", "-e", str(evidence_csv), "--hypotheses", str(primary_ranking), "-o", str(explanation_md)],
            ]
        )
    else:
        enrichment_all = out_dir / "native_enrichment_all.csv"
        enrichment_strict = out_dir / "native_enrichment_p1e-4.csv"
        graph_path = out_dir / "primary_joint_graph.graphml"
        graph_nodes = out_dir / "primary_joint_graph_nodes.csv"
        graph_audit = out_dir / "primary_joint_graph_audit.json"
        analysis_commands.extend(
            [
                [
                    "brick", "rk", "-i", str(combined), "--enrich",
                    "--background-method", "auto", "--pvalue-cutoff", str(args.pvalue_cutoff),
                    "--relation-batch-size", "10", "--relation-workers", "4", "--relation-retries", "2",
                    "-o", str(enrichment_all),
                ],
                ["brick", "drug", "filter-enrichment", "-i", str(enrichment_all), "--pvalue-cutoff", str(args.pvalue_cutoff), "-o", str(enrichment_strict)],
                ["brick", "drug", "scope", "-i", str(combined), "-o", str(scope_csv)],
                [
                    "brick", "drug", "build-graph", "-a", str(args.adata), "--genes", str(genes_path),
                    "-i", str(enrichment_strict), "--disease", args.anchor_disease,
                    "--pvalue-cutoff", str(args.pvalue_cutoff), "--identity-output", str(graph_nodes),
                    "--audit-output", str(graph_audit), "-o", str(graph_path),
                ],
                [
                    "brick-gpu", "drug", "grace-rank", "-g", str(graph_path), "--identity", str(graph_nodes),
                    "--disease", args.anchor_disease, "--seeds", args.grace_seeds,
                    "--epochs", str(args.grace_epochs), "--device", "cuda", "-o", str(out_dir),
                ],
                ["brick", "drug", "extract", "-i", str(primary_ranking), "--scope", str(scope_csv), "--top-n", str(args.top_n_drugs), "-o", str(predicted)],
                ["brick", "reason", "explain", "-e", str(evidence_csv), "--hypotheses", str(primary_ranking), "-o", str(explanation_md)],
            ]
        )

    validation_commands = []
    reference_path = Path(args.real_drugs) if args.real_drugs else out_dir / "covid_chemical_reference.csv"
    if args.disease:
        known_cmd = [
            "brick", "qr", args.disease, "--st", "Disease", "--tt", "Chemical",
            "-r", "pharmaco_effect_on", "-n", str(args.query_limit), "-o", str(reference_path),
        ]
        validation_commands.append(known_cmd)
    if args.disease or args.real_drugs:
        validation_commands.append(
            [
                "brick", "drug", "validate", "--ranking", str(primary_ranking),
                "--freeze", str(freeze_path), "--reference", str(reference_path),
                "--identity-audit", str(out_dir / "chemical_identity_audit.csv"),
                "-o", str(out_dir / "ranking_validation.csv"),
            ]
        )

    command_lists = [*query_commands, merge_cmd, *analysis_commands, *validation_commands]
    commands = [_cmd_to_sh(command) for command in command_lists]

    if args.run_query:
        for command in query_commands:
            _run(command)

    merged = _merge_csvs(query_paths, combined)
    if args.run_analysis:
        if not merged:
            raise ValueError("No query files are available for analysis")
        for command in analysis_commands:
            _run(command)

    if args.run_validation:
        if not freeze_path.exists():
            raise ValueError("Validation is blocked until ranking_frozen.sha256 exists")
        for command in validation_commands:
            _run(command)

    _commands_files(out_dir, commands)
    if getattr(args, "agent_plan", False):
        agent_plan = _write_agent_plan(out_dir, "drug", commands, question=getattr(args, "question", None))
        print(f"Wrote {agent_plan}")
    print(f"Drug workflow prepared in {out_dir} (method={args.method})")


def cmd_trajectory(args):
    out_dir = _ensure_dir(args.output_dir)
    paga_graph = out_dir / "paga.graphml"
    cell_types = out_dir / "cell_types.txt"
    kg_query = Path(args.query_csv) if args.query_csv else out_dir / "kg_develope_into.csv"
    merged = out_dir / "merged.graphml"
    commands = []

    if args.graph:
        shutil.copyfile(args.graph, paga_graph)

    if args.adata:
        import scanpy as sc
        from graphadapter import data2nxg, save_graph

        adata = sc.read(args.adata)
        if "paga" not in adata.uns:
            if args.groupby not in adata.obs:
                raise ValueError(f"groupby column not found in adata.obs: {args.groupby}")
            if "neighbors" not in adata.uns:
                sc.pp.neighbors(adata)
            sc.tl.paga(adata, groups=args.groupby)
        save_graph(data2nxg(adata, data_type="paga"), str(paga_graph))

    if paga_graph.exists():
        import networkx as nx

        graph = nx.read_graphml(paga_graph)
        _write_text(cell_types, "\n".join(str(x) for x in graph.nodes()) + "\n")

    if not args.query_csv:
        traj_cmd = _qr_command(args, cell_types, kg_query, "Cell", "Cell", relation="develope_into", task="trajectory")
        commands.append(_cmd_to_sh(traj_cmd))
    commands.append(
        f'brick ga complete_result -g "{paga_graph}" -q "{kg_query}" -o "{merged}"'
    )
    if args.root:
        commands.append(f'brick pl static -i "{merged}" -o "{out_dir / "trajectory_tree.png"}" --layout tree --root "{args.root}"')
    commands.append(f'brick pl interact -i "{merged}" -o "{out_dir / "trajectory.html"}" --width 1400 --height 1000')

    if args.run_query and not args.query_csv:
        _run(_qr_command(args, cell_types, kg_query, "Cell", "Cell", relation="develope_into", task="trajectory"))

    if args.run_analysis and paga_graph.exists() and kg_query.exists():
        _run(["brick", "ga", "complete_result", "-g", str(paga_graph), "-q", str(kg_query), "-o", str(merged)])

    _commands_files(out_dir, commands)
    if getattr(args, "agent_plan", False):
        agent_plan = _write_agent_plan(out_dir, "trajectory", commands, question=getattr(args, "question", None))
        print(f"Wrote {agent_plan}")
    print(f"Trajectory workflow prepared in {out_dir}")


def build_parser():
    parser = argparse.ArgumentParser(prog="brick workflow", description="Curated brick workflows")
    subparsers = parser.add_subparsers(dest="command")

    cell = subparsers.add_parser("celltype", help="Reference-free cell type annotation workflow")
    cell.add_argument("-a", "--adata", help="Input AnnData file")
    cell.add_argument("--marker-dir", help="Directory with marker gene txt files")
    cell.add_argument("--groupby", default="leiden", help="obs column for clusters")
    cell.add_argument("--layer", help="Layer to copy into adata.X before DEG")
    cell.add_argument("--de-method", default="wilcoxon")
    cell.add_argument("--top-n-markers", type=int, default=10)
    cell.add_argument("--query-limit", type=int, default=500)
    cell.add_argument("--hops", type=int, choices=[1, 2, 3], default=1, help="KG query hops for workflow query steps")
    cell.add_argument("--kg-policy", choices=KG_POLICY_CHOICES, default="auto", help="KG source/hop policy")
    cell.add_argument("--info-source", help="Override edge info_source filter")
    cell.add_argument("--question", help="Question/context for KG source routing")
    cell.add_argument("--policy-file", help="Optional JSON KG source policy catalog")
    cell.add_argument("--rank-metrics", default="match_count,info_source_count")
    cell.add_argument("--prune-parent", dest="prune_parent", action="store_true", default=True)
    cell.add_argument("--no-prune-parent", dest="prune_parent", action="store_false")
    cell.add_argument("--run-query", action="store_true")
    cell.add_argument("--agent-plan", action="store_true", help="Also write agent_plan.json for brick agent check/run.")
    cell.add_argument("-o", "--output-dir", default="brick_celltype_workflow")
    cell.set_defaults(func=cmd_celltype)

    enrich = subparsers.add_parser("enrichment", help="Gene set enrichment workflow")
    enrich.add_argument("--genes", help="Input gene list")
    enrich.add_argument("-a", "--adata", help="Input AnnData file")
    enrich.add_argument("--groupby", default="leiden")
    enrich.add_argument("--group", help="Group to extract genes from")
    enrich.add_argument("--layer", help="Layer to use for DEG")
    enrich.add_argument("--de-method", default="wilcoxon")
    enrich.add_argument("--top-n-genes", type=int, default=100)
    enrich.add_argument("--target-types", default=",".join(DEFAULT_TARGET_TYPES))
    enrich.add_argument("--relations", default=None, help="Override query relations. By default, policy selects target-specific KG relations.")
    enrich.add_argument("--query-limit", type=int, default=1000)
    enrich.add_argument("--hops", type=int, choices=[1, 2, 3], default=1, help="KG query hops for workflow query steps")
    enrich.add_argument("--kg-policy", choices=KG_POLICY_CHOICES, default="auto", help="KG source/hop policy")
    enrich.add_argument("--info-source", help="Override edge info_source filter")
    enrich.add_argument("--question", help="Question/context for KG source routing")
    enrich.add_argument("--policy-file", help="Optional JSON KG source policy catalog")
    enrich.add_argument("--top-k", type=int, default=15)
    enrich.add_argument("--run-query", action="store_true")
    enrich.add_argument("--run-analysis", action="store_true")
    enrich.add_argument("--agent-plan", action="store_true", help="Also write agent_plan.json for brick agent check/run.")
    enrich.add_argument("-o", "--output-dir", default="brick_enrichment_workflow")
    enrich.set_defaults(func=cmd_enrichment)

    drug = subparsers.add_parser("drug", help="Single-cell drug repositioning workflow")
    drug.add_argument("--genes", help="Input gene list")
    drug.add_argument("-a", "--adata", help="Input AnnData file")
    drug.add_argument("--groupby", default="leiden")
    drug.add_argument("--group", help="Group to extract genes from")
    drug.add_argument("--layer", help="Layer to use for DEG")
    drug.add_argument("--de-method", default="wilcoxon")
    drug.add_argument("--top-n-genes", type=int, default=100)
    drug.add_argument("--query-limit", type=int, default=500)
    drug.add_argument("--method", choices=["rank", "grace"], default="rank", help="Independent RankGraph or GPU GRACE consensus")
    drug.add_argument("--target-types", default=DEFAULT_DRUG_GRACE_TARGET_TYPES, help="Comma-separated one-hop target union for --method grace")
    drug.add_argument("--anchor-disease", default="COVID-19", help="Exact Disease node used only as the GRACE similarity anchor")
    drug.add_argument("--pvalue-cutoff", type=float, default=0.0001, help="Predeclared raw-p threshold for the strict GRACE graph")
    drug.add_argument("--grace-seeds", default="20241204,20241205,20241206", help="Fixed comma-separated GRACE seeds")
    drug.add_argument("--grace-epochs", type=int, default=1000)
    drug.add_argument("--hops", type=int, choices=[1, 2, 3], default=1, help="KG query hops for workflow query steps")
    drug.add_argument("--kg-policy", choices=KG_POLICY_CHOICES, default="auto", help="KG source/hop policy")
    drug.add_argument("--info-source", help="Override edge info_source filter")
    drug.add_argument("--question", help="Question/context for KG source routing")
    drug.add_argument("--policy-file", help="Optional JSON KG source policy catalog")
    drug.add_argument("--top-n-drugs", type=int, default=50)
    drug.add_argument("--disease", help="Disease name for known-drug lookup command")
    drug.add_argument("--real-drugs", help="Existing real-drugs CSV for validation")
    drug.add_argument("--evaluation-only", action="store_true", help="Allow a reference query only after the primary ranking is frozen")
    drug.add_argument("--run-query", action="store_true")
    drug.add_argument("--run-analysis", action="store_true")
    drug.add_argument("--run-validation", action="store_true", help="Run the post-freeze reference query and stable-ID evaluation")
    drug.add_argument("--agent-plan", action="store_true", help="Also write agent_plan.json for brick agent check/run.")
    drug.add_argument("-o", "--output-dir", default="brick_drug_workflow")
    drug.set_defaults(func=cmd_drug)

    traj = subparsers.add_parser("trajectory", help="PAGA plus KG trajectory workflow")
    traj.add_argument("-a", "--adata", help="Input AnnData file with or without PAGA")
    traj.add_argument("-g", "--graph", help="Existing PAGA/trajectory GraphML file")
    traj.add_argument("--groupby", default="cell_type")
    traj.add_argument("--query-csv", help="Existing KG query CSV")
    traj.add_argument("--root", help="Root node for tree plot")
    traj.add_argument("--query-limit", type=int, default=500)
    traj.add_argument("--hops", type=int, choices=[1, 2, 3], default=1, help="KG query hops for workflow query steps")
    traj.add_argument("--kg-policy", choices=KG_POLICY_CHOICES, default="auto", help="KG source/hop policy")
    traj.add_argument("--info-source", help="Override edge info_source filter")
    traj.add_argument("--question", help="Question/context for KG source routing")
    traj.add_argument("--policy-file", help="Optional JSON KG source policy catalog")
    traj.add_argument("--run-query", action="store_true")
    traj.add_argument("--run-analysis", action="store_true")
    traj.add_argument("--agent-plan", action="store_true", help="Also write agent_plan.json for brick agent check/run.")
    traj.add_argument("-o", "--output-dir", default="brick_trajectory_workflow")
    traj.set_defaults(func=cmd_trajectory)

    return parser


def _cli_main():
    parser = build_parser()
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit(1)
    try:
        args.func(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _cli_main()
