#!/usr/bin/env python3
"""Run strict Biomix MCQ source-aware KG-first option verification.

This wrapper launches the benchmark BRICK runner in an independent
Hermes+BRICK+real-KG mode:

question/options -> parse entities -> source-routed option-wise real KG
retrieval -> all-option evidence table -> Hermes final choice.

It does not read MCQ gold labels, Biomni predictions, consensus outputs, or
native fallback answer files during generation.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


DEFAULT_METHOD = "hermes_brick_source_aware_mcq_20260707"


def _default_benchmark_root() -> Path:
    candidates = [
        Path.cwd() / "biomics-benchmark",
        Path.cwd(),
        Path(r"D:\code\Biomics\biomics-benchmark"),
    ]
    for candidate in candidates:
        if (candidate / "run_agent_pipeline" / "brickkg" / "run_biomix_batch_glm51_brickkg.py").exists():
            return candidate.resolve()
    return candidates[-1].resolve()


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def build_command(args: argparse.Namespace) -> list[str]:
    runner = args.benchmark_root / "run_agent_pipeline" / "brickkg" / "run_biomix_batch_glm51_brickkg.py"
    if not runner.exists():
        raise FileNotFoundError(f"Biomix BRICK runner not found: {runner}")
    command = [
        sys.executable,
        str(runner),
        "--tasks",
        "mcq",
        "--method",
        args.method,
        "--workers",
        str(args.workers),
        "--model",
        args.model,
        "--kg-limit",
        str(args.kg_limit),
        "--request-timeout",
        str(args.request_timeout),
        "--kg-request-timeout",
        str(args.kg_request_timeout),
        "--max-tokens",
        str(args.max_tokens),
        "--mcq-kg-mode",
        "direct_only",
        "--mcq-decision-mode",
        "realkg_verifier",
        "--mcq-query-hops",
        str(args.max_hops),
        "--mcq-fallback-stop-condition",
        args.fallback_stop_condition,
        "--mcq-no-native-fallback",
        "--mcq-reverse-query",
        "--native-mcq-from",
        "",
        "--mcq-strong-choice-policy",
        "llm",
        "--mcq-score-policy",
        "dual_support",
        "--mcq-verifier-evidence-view",
        args.verifier_evidence_view,
        "--mcq-source-routing",
        args.source_routing,
        "--mcq-source-stage-mode",
        args.source_stage_mode,
        "--kg-relation-policy",
        args.kg_relation_policy,
        "--weak-evidence-policy",
        args.weak_evidence_policy,
        "--mcq-scope-aggregation",
        "option",
        "--mcq-candidate-stage-policy",
        "target_retained_only",
        "--mcq-evidence-ranking",
        "on",
        "--brick-skill-dir",
        str(args.brick_skill_dir),
    ]
    if args.limit is not None:
        command.extend(["--limit", str(args.limit)])
    if args.source_broad_fallback:
        command.append("--mcq-source-broad-fallback")
    if args.native_if_no_kg_support:
        command.append("--mcq-native-if-no-kg-support")
    if args.native_if_low_kg_confidence:
        command.append("--mcq-native-if-low-kg-confidence")
    if args.overwrite:
        command.append("--overwrite")
    return command


def run_realkg(args: argparse.Namespace) -> int:
    command = build_command(args)
    manifest_dir = args.benchmark_root / "run_agent_pipeline" / "brickkg" / "summaries" / args.method
    manifest = {
        "dataset": "Biomix",
        "task": "mcq_5_choice",
        "method": args.method,
        "solver": "hermes_brick_source_aware_option_verifier"
        if args.source_routing == "source_aware"
        else "hermes_brick_realkg_option_verifier",
        "prediction_source": "independent",
        "gold_access_during_generation": "none",
        "forbidden_sources": [
            "Biomni primary",
            "question-only consensus",
            "native fallback answer files",
            "gold-aware standardizer",
            "mcq_biomix.csv correct_answer",
        ],
        "pipeline": [
            "question/options parse",
            "question-type/source-aware edge info_source routing",
            "option-wise BRICK real KG query",
            "evidence table per option",
            "Hermes final choice from evidence table",
        ],
        "policy": {
            "source_routing": args.source_routing,
            "source_stage_mode": args.source_stage_mode,
            "source_broad_fallback": args.source_broad_fallback,
            "max_hops": args.max_hops,
            "fallback_stop_condition": args.fallback_stop_condition,
            "kg_relation_policy": args.kg_relation_policy,
            "weak_evidence_policy": args.weak_evidence_policy,
            "verifier_evidence_view": args.verifier_evidence_view,
            "native_if_no_kg_support": args.native_if_no_kg_support,
            "native_if_low_kg_confidence": args.native_if_low_kg_confidence,
            "native_fallback": False,
            "biomni": False,
            "consensus": False,
        },
        "outputs": {
            "answers": str(args.benchmark_root / "agent_answer" / "biomix" / args.method / "mcq_biomix_pure_question.json"),
            "details": str(manifest_dir / "mcq_biomix_pure_question.details.json"),
            "kg_summary": str(manifest_dir / "kg_summary.json"),
        },
        "command": command,
    }
    if args.dry_run:
        print("DRY RUN " + " ".join(command), flush=True)
        write_json(manifest_dir / "realkg_manifest.json", manifest)
        return 0

    completed = subprocess.run(command, cwd=str(args.benchmark_root), text=True, check=False)
    write_json(manifest_dir / "realkg_manifest.json", manifest)
    return int(completed.returncode)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Biomix MCQ predictions with strict Hermes+BRICK real-KG option verification.")
    parser.add_argument("--benchmark-root", type=Path, default=_default_benchmark_root())
    parser.add_argument("--brick-skill-dir", type=Path, default=Path(r"D:\code\Biomics\BiOmics_BRICK_SKILL_fix"))
    parser.add_argument("--method", default=DEFAULT_METHOD)
    parser.add_argument("--model", default="deepseek-v4-pro")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--kg-limit", type=int, default=30)
    parser.add_argument("--max-hops", type=int, default=1)
    parser.add_argument("--fallback-stop-condition", choices=["nonempty", "candidate_eligible", "all_stages"], default="candidate_eligible")
    parser.add_argument("--source-routing", choices=["off", "source_aware"], default="source_aware")
    parser.add_argument("--source-stage-mode", choices=["merged", "layered"], default="merged")
    parser.add_argument("--source-broad-fallback", action="store_true")
    parser.add_argument("--verifier-evidence-view", choices=["all", "answer_eligible"], default="all")
    parser.add_argument("--native-if-no-kg-support", action="store_true")
    parser.add_argument("--native-if-low-kg-confidence", action="store_true")
    parser.add_argument("--kg-relation-policy", default="curated_preferred")
    parser.add_argument("--weak-evidence-policy", default="snippet_supported")
    parser.add_argument("--request-timeout", type=float, default=120)
    parser.add_argument("--kg-request-timeout", type=float, default=30)
    parser.add_argument("--max-tokens", type=int, default=256)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    args.benchmark_root = args.benchmark_root.resolve()
    args.brick_skill_dir = args.brick_skill_dir.resolve()
    return args


def _cli_main() -> int:
    try:
        return run_realkg(parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(_cli_main())
