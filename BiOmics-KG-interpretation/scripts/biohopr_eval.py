#!/usr/bin/env python3
"""Run BioHopR benchmark evaluation from the BiOmics-BRICK skill layer."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import biohopr_paths


SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_DIR = SCRIPT_DIR.parent


def _default_benchmark_root() -> Path:
    candidates = [
        Path.cwd() / "biomics-benchmark",
        Path.cwd(),
        Path(r"D:\code\Biomics\biomics-benchmark"),
    ]
    for candidate in candidates:
        if (candidate / "run_agent_pipeline" / "brickkg" / "run_biohopr_graph_answers.py").exists():
            return candidate.resolve()
    return candidates[-1]


def _add_value(cmd: list[str], flag: str, value: Any) -> None:
    if value is not None and value != "":
        cmd.extend([flag, str(value)])


def _default_method_for_preset(preset: str) -> str:
    if preset in {"balanced", "strict_e2e_question_only"}:
        return "strict_e2e_question_only"
    return f"deepseek_v4_pro_D_main_v2_biohopr_{preset}_skill"


def build_runner_command(args: argparse.Namespace) -> tuple[list[str], dict[str, Any]]:
    benchmark_root = args.benchmark_root.resolve()
    runner = benchmark_root / "run_agent_pipeline" / "brickkg" / "run_biohopr_graph_answers.py"
    if not runner.exists():
        raise FileNotFoundError(f"BioHopR benchmark runner not found: {runner}")

    policy = biohopr_paths.biohopr_benchmark_policy(args.preset)
    method = args.method or _default_method_for_preset(args.preset)
    native_answer_dir = args.native_answer_dir or (
        benchmark_root / "agent_answer" / "BioHopR" / "deepseek_v4_pro_plain"
    )
    gold_dir = args.gold_dir or (benchmark_root / "evaluate_pipeline" / "BioHopR" / "gold")

    cmd = [
        sys.executable,
        str(runner),
        "--method",
        method,
        "--workers",
        str(args.workers),
        "--kg-request-timeout",
        str(args.kg_request_timeout),
        "--kg-limit",
        str(policy["kg_limit"]),
        "--max-answers",
        str(policy["max_answers"]),
        "--kg-answer-strength",
        str(policy["kg_answer_strength"]),
        "--path-co-mentioned-policy",
        str(policy["path_co_mentioned_policy"]),
        "--co-mentioned-fallback-min-score",
        str(policy["co_mentioned_fallback_min_score"]),
        "--answer-selection-mode",
        str(policy["answer_selection_mode"]),
        "--score-drop",
        str(policy["score_drop"]),
        "--min-answers",
        str(policy["min_answers"]),
        "--answer-policy",
        str(policy["answer_policy"]),
        "--kg-add-max-1hop",
        str(policy["kg_add_max_1hop"]),
        "--kg-add-max-2hop",
        str(policy["kg_add_max_2hop"]),
        "--kg-add-min-score-1hop",
        str(policy["kg_add_min_score_1hop"]),
        "--kg-add-min-score-2hop",
        str(policy["kg_add_min_score_2hop"]),
        "--kg-add-max-weak-1hop",
        str(policy["kg_add_max_weak_1hop"]),
        "--kg-add-max-weak-2hop",
        str(policy["kg_add_max_weak_2hop"]),
        "--biohopr-source-policy",
        str(policy["biohopr_source_policy"]),
        "--weak-source-rescue-max-1hop",
        str(policy["weak_source_rescue_max_1hop"]),
        "--weak-source-rescue-max-2hop",
        str(policy["weak_source_rescue_max_2hop"]),
        "--brick-skill-dir",
        str(SKILL_DIR),
        "--query-source",
        str(policy["query_source"]),
        "--gold-dir",
        str(gold_dir),
        "--entity-linking",
        str(policy.get("entity_linking", "off")),
        "--max-query-entities",
        str(policy.get("max_query_entities", 1)),
        "--max-bridge-entities",
        str(policy.get("max_bridge_entities", 25)),
    ]
    _add_value(cmd, "--entity-linking-1hop", policy.get("entity_linking_1hop"))
    _add_value(cmd, "--entity-linking-2hop", policy.get("entity_linking_2hop"))
    _add_value(cmd, "--max-query-entities-1hop", policy.get("max_query_entities_1hop"))
    _add_value(cmd, "--max-query-entities-2hop", policy.get("max_query_entities_2hop"))

    if args.splits:
        cmd.extend(["--splits", *args.splits])
    else:
        cmd.extend(["--split-group", args.split_group])
    if args.limit is not None:
        cmd.extend(["--limit", str(args.limit)])
    if args.overwrite:
        cmd.append("--overwrite")
    if args.dry_run:
        cmd.append("--dry-run")
    if args.no_evaluate:
        cmd.append("--no-evaluate")
    if policy.get("kg_add_allow_weak"):
        cmd.append("--kg-add-allow-weak")
    if policy.get("type_fallback"):
        cmd.append("--type-fallback")
    if policy.get("discover_bridge_fallback"):
        cmd.append("--discover-bridge-fallback")
    if policy.get("type_fallback_1hop"):
        cmd.append("--type-fallback-1hop")
    if policy.get("type_fallback_2hop"):
        cmd.append("--type-fallback-2hop")
    if policy.get("discover_bridge_fallback_1hop"):
        cmd.append("--discover-bridge-fallback-1hop")
    if policy.get("discover_bridge_fallback_2hop"):
        cmd.append("--discover-bridge-fallback-2hop")
    if policy.get("relax_source_type_fallback"):
        cmd.append("--relax-source-type-fallback")
    if policy.get("relax_target_type_fallback"):
        cmd.append("--relax-target-type-fallback")
    _add_value(cmd, "--native-answer-dir", native_answer_dir)
    _add_value(cmd, "--min-path-score", policy.get("min_path_score"))
    _add_value(cmd, "--weak-source-rescue-min-score-1hop", policy.get("weak_source_rescue_min_score_1hop"))
    _add_value(cmd, "--weak-source-rescue-min-score-2hop", policy.get("weak_source_rescue_min_score_2hop"))
    _add_value(cmd, "--kg-info-source", policy.get("kg_info_source"))

    run_manifest = {
        "preset": args.preset,
        "policy": policy,
        "benchmark_root": str(benchmark_root),
        "method": method,
        "runner": str(runner),
        "native_answer_dir": str(native_answer_dir),
        "gold_dir": str(gold_dir),
        "command": cmd,
    }
    return cmd, run_manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run BioHopR evaluation using BiOmics-BRICK skill presets."
    )
    parser.add_argument(
        "--preset",
        choices=sorted(biohopr_paths.BIOHOPR_BENCHMARK_POLICIES),
        default="balanced",
        help="Skill-owned BioHopR benchmark policy preset.",
    )
    parser.add_argument("--benchmark-root", type=Path, default=_default_benchmark_root())
    parser.add_argument("--method", default=None)
    parser.add_argument("--split-group", choices=["single", "full", "all"], default="single")
    parser.add_argument(
        "--splits",
        nargs="+",
        choices=["1hop_single", "2hop_single", "1hop_full", "2hop_full"],
        default=None,
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--kg-request-timeout", type=float, default=30.0)
    parser.add_argument("--native-answer-dir", type=Path, default=None)
    parser.add_argument("--gold-dir", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-evaluate", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--print-policy",
        action="store_true",
        help="Print the resolved skill preset before running.",
    )
    return parser.parse_args()


def _cli_main() -> int:
    args = parse_args()
    cmd, manifest = build_runner_command(args)
    if args.print_policy or args.dry_run:
        print(json.dumps(manifest, ensure_ascii=False, indent=2), flush=True)
    print("RUN " + " ".join(str(part) for part in cmd), flush=True)
    completed = subprocess.run(cmd, cwd=str(args.benchmark_root.resolve()))
    return int(completed.returncode)


if __name__ == "__main__":
    raise SystemExit(_cli_main())
