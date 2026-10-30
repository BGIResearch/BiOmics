#!/usr/bin/env python3
"""Run the unified BioASQ evaluation bundle from the BiOmics-BRICK skill layer.

The unified BioASQ bundle keeps one public dataset/evaluator entrypoint while
recording task-level protocols. YesNo is read from the strict no-gold
prior-guided self-correction metrics, and Factoid/List are reproduced through
the legacy gold-aware standardizer plus the official BioASQ Java evaluator.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


METRIC_NAMES = [
    "YesNo_Acc",
    "Factoid_Strict_Acc",
    "Factoid_Lenient_Acc",
    "Factoid_MRR",
    "List_Precision",
    "List_Recall",
    "List_F1",
    "YesNo_macroF1",
    "YesNo_F1_yes",
    "YesNo_F1_no",
]

REPORT_METRIC_NAMES = [
    "Factoid_Strict_Acc",
    "Factoid_Lenient_Acc",
    "Factoid_MRR",
    "List_Precision",
    "List_Recall",
    "List_F1",
]

STRICT_YESNO_METRIC_NAMES = [
    "YesNo_Acc",
    "YesNo_macroF1",
    "YesNo_F1_yes",
    "YesNo_F1_no",
]

YESNO_VARIANT = "prior_guided_self_correction"
YESNO_PROTOCOL = "strict_nogold_prior_guided_self_correction"
YESNO_STATUS = "strict_nogold_prior_guided_self_correction_not_heldout"


SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_ROOT = SCRIPT_DIR.parents[2]


def _default_bioasq_root() -> Path:
    candidates = [
        BENCHMARK_ROOT / "evaluate_pipeline" / "BioASQ",
        Path.cwd() / "BioASQ Benchmark",
        Path.cwd(),
        Path(r"D:\code\Biomics\BioASQ Benchmark"),
    ]
    for candidate in candidates:
        if (candidate / "scripts" / "extract_and_standardize.py").exists():
            return candidate.resolve()
    return candidates[-1].resolve()


def _default_gold(root: Path) -> Path:
    candidates = [
        BENCHMARK_ROOT / "evaluate_pipeline" / "BioASQ" / "gold" / "answer.json",
        Path.cwd() / "evaluate_pipeline" / "BioASQ" / "gold" / "answer.json",
        Path.cwd() / "biomics-benchmark" / "evaluate_pipeline" / "BioASQ" / "gold" / "answer.json",
        root / "data" / "answer.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def _default_input(root: Path) -> Path | None:
    candidates = [
        BENCHMARK_ROOT / "agent_answer" / "BioASQ" / "deepseek_v4_pro_D" / "deepseek_v4_pro_D.results.json",
        root / "results_old" / "deepseek_v4_pro_D.results.json",
        root / "results_old" / "deepseek_v4_pro_D_biomics_githubinstall_latestkg_20260706.results.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return None


def _prefix_from_input(input_path: Path, output_dir: Path, tag: str) -> Path:
    name = input_path.name
    if name.endswith(".json"):
        name = name[:-5]
    if name.endswith(".results"):
        name = name[:-8]
    return output_dir / f"{name}.{tag}"


def _path_from_prefix(prefix: Path, suffix: str) -> Path:
    return Path(str(prefix) + suffix)


def _existing_legacy_outputs(input_path: Path) -> tuple[Path | None, Path | None]:
    name = input_path.name
    if name.endswith(".json"):
        name = name[:-5]
    if name.endswith(".results"):
        name = name[:-8]
    standardized = input_path.parent / f"{name}.legacy_goldaware_standardized.json"
    eval_output = input_path.parent / f"{name}.legacy_goldaware_standardized.eval.json"
    if standardized.exists() and eval_output.exists():
        return standardized.resolve(), eval_output.resolve()
    return None, None


def _resolve_args(args: argparse.Namespace) -> argparse.Namespace:
    root = args.bioasq_root.resolve()
    args.bioasq_root = root
    if args.input is None:
        default_input = _default_input(root)
        if default_input is None:
            raise FileNotFoundError(
                "No default new-skill BioASQ result was found. Pass --input explicitly."
            )
        args.input = default_input
    else:
        args.input = args.input.resolve()

    if args.output_dir is not None:
        args.output_dir = args.output_dir.resolve()
    elif args.metrics_json is not None:
        args.output_dir = args.metrics_json.parent.resolve()
    else:
        args.output_dir = args.input.parent.resolve()
    tag = "legacy_goldaware_standardized"
    args.output_prefix = (args.output_prefix or _prefix_from_input(args.input, args.output_dir, tag)).resolve()
    existing_standardized, existing_eval = _existing_legacy_outputs(args.input)
    reuse_existing_legacy = args.skip_standardizer or args.skip_standardizer_if_exists
    if args.standardized_output is None and existing_standardized is not None and reuse_existing_legacy:
        args.standardized_output = existing_standardized
    else:
        args.standardized_output = (args.standardized_output or _path_from_prefix(args.output_prefix, ".json")).resolve()
    if args.eval_output is None and existing_eval is not None and reuse_existing_legacy:
        args.eval_output = existing_eval
    else:
        args.eval_output = (args.eval_output or _path_from_prefix(args.output_prefix, ".eval.json")).resolve()
    args.metrics_json = (args.metrics_json or _path_from_prefix(args.output_prefix, ".metrics.json")).resolve()
    args.metrics_csv = (args.metrics_csv or _path_from_prefix(args.output_prefix, ".metrics.csv")).resolve()
    args.manifest = (args.manifest or _path_from_prefix(args.output_prefix, ".manifest.json")).resolve()
    args.yesno_metrics_json = (
        args.yesno_metrics_json
        or (args.output_dir / "bioasq_yesno_prior_guarded.official_metrics.json")
    ).resolve()
    args.yesno_eval_json = (
        args.yesno_eval_json
        or (args.output_dir / "bioasq_yesno_prior_guarded.eval.json")
    ).resolve()
    args.unified_eval_output = (
        args.unified_eval_output
        or (args.output_dir / "bioasq_unified.eval.json")
    ).resolve()
    args.unified_metrics_json = (
        args.unified_metrics_json
        or (args.output_dir / "bioasq_unified.metrics.json")
    ).resolve()
    args.unified_metrics_csv = (
        args.unified_metrics_csv
        or (args.output_dir / "bioasq_unified.metrics.csv")
    ).resolve()
    args.unified_manifest = (
        args.unified_manifest
        or (args.output_dir / "bioasq_unified.manifest.json")
    ).resolve()
    args.unified_official_metrics_json = (
        args.unified_official_metrics_json
        or (args.output_dir / "bioasq_unified.official_metrics.json")
    ).resolve()
    args.gold = (args.gold or _default_gold(root)).resolve()
    args.extract_script = (args.extract_script or (root / "scripts" / "extract_and_standardize.py")).resolve()
    args.evaluator_jar = (
        args.evaluator_jar
        or (root / "Evaluation-Measures" / "flat" / "BioASQEvaluation" / "dist" / "BioASQEvaluation.jar")
    ).resolve()
    return args


def _check_paths(args: argparse.Namespace) -> None:
    required = {
        "input": args.input,
        "gold": args.gold,
        "extract_script": args.extract_script,
        "evaluator_jar": args.evaluator_jar,
    }
    if not args.allow_missing_yesno:
        required["yesno_metrics_json"] = args.yesno_metrics_json
        required["yesno_eval_json"] = args.yesno_eval_json
    missing = [f"{name}: {path}" for name, path in required.items() if not path.exists()]
    if missing:
        raise FileNotFoundError("Required BioASQ evaluation files are missing:\n" + "\n".join(missing))


def _run(cmd: list[str], cwd: Path, *, capture: bool = False) -> subprocess.CompletedProcess[str]:
    printable = " ".join(str(part) for part in cmd)
    print("RUN " + printable, flush=True)
    return subprocess.run(
        cmd,
        cwd=str(cwd),
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=capture,
        check=False,
    )


def run_legacy_goldaware_standardizer(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.skip_standardizer or (
        args.skip_standardizer_if_exists and args.standardized_output.exists() and args.eval_output.exists()
    ):
        if not args.eval_output.exists():
            raise FileNotFoundError(f"Cannot skip BioASQ standardizer because eval output is missing: {args.eval_output}")
        print(f"Reusing BioASQ standardized/eval output: {args.eval_output}", flush=True)
        return
    cmd = [
        sys.executable,
        str(args.extract_script),
        "--input",
        str(args.input),
        "--output",
        str(args.standardized_output),
        "--eval-output",
        str(args.eval_output),
        "--gold",
        str(args.gold),
        "--no-extract",
        "--redo-standardize",
        "--workers",
        str(args.workers),
        "--request-timeout",
        str(args.request_timeout),
        "--max-tokens",
        str(args.max_tokens),
        "--model",
        str(args.standardizer_model),
        "--base-url",
        str(args.base_url),
        "--api-key-env",
        str(args.api_key_env),
        "--gold-aware-standardizer",
        str(args.gold_aware_standardizer),
        "--gold-aware-list-unmatched",
        str(args.gold_aware_list_unmatched),
    ]
    if args.enable_thinking:
        cmd.append("--enable-thinking")
    if args.dry_run:
        print("DRY-RUN standardizer command:")
        print(" ".join(str(part) for part in cmd))
        return
    completed = _run(cmd, args.bioasq_root)
    if completed.returncode != 0:
        raise RuntimeError(f"BioASQ standardizer failed with exit code {completed.returncode}")


def parse_evaluator_output(text: str) -> dict[str, float]:
    numbers: list[float] = []
    for token in re.split(r"\s+", text.strip()):
        if not token:
            continue
        try:
            numbers.append(float(token))
        except ValueError:
            continue
    if len(numbers) < len(METRIC_NAMES):
        raise ValueError(f"Could not parse BioASQ evaluator output: {text!r}")
    return dict(zip(METRIC_NAMES, numbers[: len(METRIC_NAMES)]))


def run_evaluator(args: argparse.Namespace, eval_path: Path | None = None) -> dict[str, dict[str, float]]:
    results: dict[str, dict[str, float]] = {}
    if args.dry_run:
        return results
    eval_path = eval_path or args.eval_output
    for edition in args.editions:
        cmd = [
            "java",
            "-Xmx10G",
            "-cp",
            str(args.evaluator_jar),
            "evaluation.EvaluatorTask1b",
            "-phaseB",
            "-e",
            str(edition),
            str(args.gold),
            str(eval_path),
        ]
        completed = _run(cmd, args.bioasq_root, capture=True)
        if completed.stdout:
            print(completed.stdout.strip(), flush=True)
        if completed.stderr:
            print(completed.stderr.strip(), file=sys.stderr, flush=True)
        if completed.returncode != 0:
            raise RuntimeError(f"BioASQ evaluator edition {edition} failed with exit code {completed.returncode}")
        results[str(edition)] = parse_evaluator_output(completed.stdout)
    return results


def write_unified_eval_input(args: argparse.Namespace) -> None:
    legacy_payload = read_json(args.eval_output)
    yesno_payload = read_json(args.yesno_eval_json)
    yesno_by_id = {
        question.get("id"): question
        for question in yesno_payload.get("questions", [])
        if question.get("type") == "yesno"
    }
    replaced = 0
    unified_questions = []
    for question in legacy_payload.get("questions", []):
        if question.get("type") == "yesno" and question.get("id") in yesno_by_id:
            unified_questions.append(yesno_by_id[question["id"]])
            replaced += 1
        else:
            unified_questions.append(question)
    if replaced != len(yesno_by_id):
        raise ValueError(
            f"Unified BioASQ eval input replaced {replaced} YesNo rows, "
            f"but {len(yesno_by_id)} prior-guided self-correction YesNo rows were available."
        )
    payload = dict(legacy_payload)
    payload["questions"] = unified_questions
    write_json(args.unified_eval_output, payload)


def read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as f:
        return json.load(f)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(BENCHMARK_ROOT.resolve())).replace("\\", "/")
    except ValueError:
        return str(path)


def _task_from_metric(metric: str) -> str:
    return metric.split("_", 1)[0]


def _write_rows_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    preferred = [
        "dataset",
        "task",
        "metric",
        "method",
        "variant",
        "edition",
        "value",
        "status",
        "protocol",
        "path",
    ]
    fields = sorted({key for row in rows for key in row})
    fields = [field for field in preferred if field in fields] + [field for field in fields if field not in preferred]
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_metrics(args: argparse.Namespace, metrics: dict[str, dict[str, float]]) -> None:
    report_metrics = {
        edition: {name: values[name] for name in REPORT_METRIC_NAMES if name in values}
        for edition, values in metrics.items()
    }
    payload = {
        "mode": "legacy_goldaware",
        "warning": (
            "Legacy-comparable evaluation uses gold exact answers during standardization. "
            "Use only for Factoid/List comparison with prior BioASQ tables, not strict no-gold inference. "
            "Legacy YesNo metrics are intentionally omitted; use the strict prior-guided self-correction "
            "YesNo result instead."
        ),
        "input": str(args.input),
        "standardized_output": str(args.standardized_output),
        "eval_output": str(args.eval_output),
        "metrics": report_metrics,
        "omitted_metrics": ["YesNo_Acc", "YesNo_macroF1", "YesNo_F1_yes", "YesNo_F1_no"],
    }
    write_json(args.metrics_json, payload)
    with args.metrics_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["edition", *REPORT_METRIC_NAMES])
        for edition, values in metrics.items():
            writer.writerow([edition, *[values[name] for name in REPORT_METRIC_NAMES]])


def write_manifest(args: argparse.Namespace, metrics: dict[str, dict[str, float]]) -> None:
    report_metrics = {
        edition: {name: values[name] for name in REPORT_METRIC_NAMES if name in values}
        for edition, values in metrics.items()
    }
    manifest = {
        "mode": "legacy_goldaware",
        "scope": "BioASQ Factoid/List only; legacy YesNo metrics are intentionally omitted.",
        "bioasq_root": str(args.bioasq_root),
        "input": str(args.input),
        "gold": str(args.gold),
        "extract_script": str(args.extract_script),
        "evaluator_jar": str(args.evaluator_jar),
        "standardizer_model": args.standardizer_model,
        "gold_aware_standardizer": args.gold_aware_standardizer,
        "gold_aware_list_unmatched": args.gold_aware_list_unmatched,
        "editions": args.editions,
        "outputs": {
            "standardized": str(args.standardized_output),
            "eval": str(args.eval_output),
            "metrics_json": str(args.metrics_json),
            "metrics_csv": str(args.metrics_csv),
        },
        "metrics": report_metrics,
        "omitted_metrics": ["YesNo_Acc", "YesNo_macroF1", "YesNo_F1_yes", "YesNo_F1_no"],
    }
    write_json(args.manifest, manifest)


def write_unified_official_metrics(
    args: argparse.Namespace, unified_full_metrics: dict[str, dict[str, float]]
) -> None:
    payload = {
        "mode": "bioasq_unified_official_java",
        "method": "deepseek_v4_pro_D",
        "prediction_file": display_path(args.unified_eval_output),
        "gold_file": display_path(args.gold),
        "java_output_order": METRIC_NAMES,
        "metrics": unified_full_metrics,
        "note": (
            "This official Java output is run on the synchronized BioASQ evaluator input: "
            "YesNo rows come from the strict prior-guided self-correction 0.9146 artifact, while "
            "Factoid/List rows come from the uploaded legacy gold-aware standardized input."
        ),
        "generation_gold_used_for_yesno": False,
        "gold_used_for_factoid_list_standardization": True,
    }
    write_json(args.unified_official_metrics_json, payload)


def write_unified_outputs(
    args: argparse.Namespace,
    metrics: dict[str, dict[str, float]],
    unified_full_metrics: dict[str, dict[str, float]] | None = None,
) -> None:
    yesno_payload: dict[str, Any] = {}
    if args.yesno_metrics_json.exists():
        yesno_payload = read_json(args.yesno_metrics_json)
    elif not args.allow_missing_yesno:
        raise FileNotFoundError(f"Strict YesNo metrics are missing: {args.yesno_metrics_json}")

    rows: list[dict[str, Any]] = []
    legacy_method = Path(args.input).name.replace(".results.json", "")
    for edition, metric_map in metrics.items():
        for metric in REPORT_METRIC_NAMES:
            value = metric_map.get(metric)
            if not isinstance(value, (int, float)):
                continue
            task = _task_from_metric(metric)
            rows.append(
                {
                    "dataset": "BioASQ",
                    "task": task,
                    "metric": metric,
                    "method": legacy_method,
                    "variant": "legacy_goldaware_standardization",
                    "edition": str(edition),
                    "value": value,
                    "status": "legacy_goldaware_not_strict",
                    "protocol": "legacy_goldaware_factoid_list",
                    "path": display_path(args.metrics_json),
                }
            )

    primary_edition = str(args.editions[0])
    yesno_metrics = (
        {
            metric: unified_full_metrics[primary_edition][metric]
            for metric in STRICT_YESNO_METRIC_NAMES
            if unified_full_metrics
            and primary_edition in unified_full_metrics
            and metric in unified_full_metrics[primary_edition]
        }
        if unified_full_metrics
        else yesno_payload.get("metrics", {}) if yesno_payload else {}
    )
    for metric in STRICT_YESNO_METRIC_NAMES:
        value = yesno_metrics.get(metric)
        if not isinstance(value, (int, float)):
            continue
        rows.append(
            {
                "dataset": "BioASQ",
                "task": "YesNo",
                "metric": metric,
                "method": "deepseek_v4_pro_D",
                "variant": YESNO_VARIANT,
                "value": value,
                "status": YESNO_STATUS,
                "protocol": YESNO_PROTOCOL,
                "path": display_path(args.yesno_metrics_json),
            }
        )

    unified_metrics: dict[str, Any] = {
        "YesNo": {metric: yesno_metrics[metric] for metric in STRICT_YESNO_METRIC_NAMES if metric in yesno_metrics},
        "Factoid": {
            edition: {
                metric: value
                for metric, value in metric_map.items()
                if metric.startswith("Factoid_") and metric in REPORT_METRIC_NAMES
            }
            for edition, metric_map in metrics.items()
        },
        "List": {
            edition: {
                metric: value
                for metric, value in metric_map.items()
                if metric.startswith("List_") and metric in REPORT_METRIC_NAMES
            }
            for edition, metric_map in metrics.items()
        },
    }
    payload = {
        "mode": "bioasq_unified_task_protocols",
        "method": "deepseek_v4_pro_D",
        "warning": (
            "This is one BioASQ result bundle with task-level protocols. YesNo is strict/no-gold "
            "during generation but not held-out. Factoid/List are legacy-comparable and gold-aware "
            "during standardization, so they are not strict no-gold claims."
        ),
        "task_protocols": {
            "YesNo": {
                "protocol": YESNO_PROTOCOL,
                "variant": YESNO_VARIANT,
                "gold_used_during_generation": False,
                "gold_used_only_after_prediction": True,
                "uses_prior_prediction": True,
                "note": (
                    "The artifact filenames retain the historical prior_guarded name, "
                    "but the reported method label is prior_guided_self_correction "
                    "because generation uses the previous raw answer as a prior."
                ),
                "heldout_available": bool(yesno_payload.get("heldout_available", False)),
                "heldout_assessment": yesno_payload.get("heldout_assessment"),
            },
            "Factoid": {
                "protocol": "legacy_goldaware_factoid_list",
                "gold_used_during_standardization": True,
                "strict_no_gold": False,
            },
            "List": {
                "protocol": "legacy_goldaware_factoid_list",
                "gold_used_during_standardization": True,
                "strict_no_gold": False,
            },
        },
        "inputs": {
            "bioasq_prediction_input": display_path(args.input),
            "gold": display_path(args.gold),
            "legacy_standardized_eval": display_path(args.eval_output),
            "yesno_prior_guarded_eval": display_path(args.yesno_eval_json),
            "unified_eval_output": display_path(args.unified_eval_output),
            "yesno_metrics_json": display_path(args.yesno_metrics_json),
        },
        "outputs": {
            "legacy_metrics_json": display_path(args.metrics_json),
            "legacy_metrics_csv": display_path(args.metrics_csv),
            "legacy_manifest": display_path(args.manifest),
            "unified_metrics_json": display_path(args.unified_metrics_json),
            "unified_metrics_csv": display_path(args.unified_metrics_csv),
            "unified_manifest": display_path(args.unified_manifest),
            "unified_official_metrics_json": display_path(args.unified_official_metrics_json),
        },
        "metrics": unified_metrics,
        "rows": rows,
    }
    write_json(args.unified_metrics_json, payload)
    _write_rows_csv(args.unified_metrics_csv, rows)
    manifest = {
        **payload,
        "bioasq_root": display_path(args.bioasq_root),
        "extract_script": display_path(args.extract_script),
        "evaluator_jar": display_path(args.evaluator_jar),
        "standardizer_model": args.standardizer_model,
        "gold_aware_standardizer": args.gold_aware_standardizer,
        "gold_aware_list_unmatched": args.gold_aware_list_unmatched,
        "editions": args.editions,
    }
    write_json(args.unified_manifest, manifest)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one BioASQ evaluation bundle with task-level protocol labels."
    )
    parser.add_argument("--bioasq-root", type=Path, default=_default_bioasq_root())
    parser.add_argument("--input", type=Path, default=None, help="New-skill BioASQ results JSON.")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--output-prefix", type=Path, default=None)
    parser.add_argument("--standardized-output", type=Path, default=None)
    parser.add_argument("--eval-output", type=Path, default=None)
    parser.add_argument("--metrics-json", type=Path, default=None)
    parser.add_argument("--metrics-csv", type=Path, default=None)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--yesno-metrics-json", type=Path, default=None)
    parser.add_argument("--yesno-eval-json", type=Path, default=None)
    parser.add_argument("--unified-eval-output", type=Path, default=None)
    parser.add_argument("--unified-metrics-json", type=Path, default=None)
    parser.add_argument("--unified-metrics-csv", type=Path, default=None)
    parser.add_argument("--unified-manifest", type=Path, default=None)
    parser.add_argument("--unified-official-metrics-json", type=Path, default=None)
    parser.add_argument("--gold", type=Path, default=None)
    parser.add_argument("--extract-script", type=Path, default=None)
    parser.add_argument("--evaluator-jar", type=Path, default=None)
    parser.add_argument("--editions", nargs="+", default=["8", "5"], choices=["8", "5"])
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--request-timeout", type=int, default=240)
    parser.add_argument("--max-tokens", type=int, default=768)
    parser.add_argument("--standardizer-model", default="glm-5.1")
    parser.add_argument("--base-url", default="https://dcsapi.dcs.cloud/api/aigress/unified/v1")
    parser.add_argument("--api-key-env", default="BRICK_API_KEY")
    parser.add_argument("--enable-thinking", action="store_true")
    parser.add_argument(
        "--skip-standardizer",
        action="store_true",
        help="Reuse --eval-output instead of running extract_and_standardize.py.",
    )
    parser.add_argument(
        "--skip-standardizer-if-exists",
        action="store_true",
        help="Reuse standardized/eval outputs when both already exist; otherwise run the standardizer.",
    )
    parser.add_argument(
        "--gold-aware-standardizer",
        choices=["off", "v2", "oracle"],
        default="v2",
        help="Use v2 to match the previous BioASQ D-group evaluation style.",
    )
    parser.add_argument(
        "--gold-aware-list-unmatched",
        choices=["keep", "drop_if_matched", "drop"],
        default="drop_if_matched",
    )
    parser.add_argument(
        "--allow-missing-yesno",
        action="store_true",
        help="Write a legacy-only bundle if the strict YesNo metrics JSON is absent.",
    )
    parser.add_argument("--dry-run", action="store_true")
    return _resolve_args(parser.parse_args())


def _cli_main() -> int:
    try:
        args = parse_args()
        _check_paths(args)
        run_legacy_goldaware_standardizer(args)
        metrics = run_evaluator(args)
        if not args.dry_run:
            write_unified_eval_input(args)
            unified_full_metrics = run_evaluator(args, eval_path=args.unified_eval_output)
            write_unified_official_metrics(args, unified_full_metrics)
            write_metrics(args, metrics)
            write_manifest(args, metrics)
            write_unified_outputs(args, metrics, unified_full_metrics)
            edition = args.editions[0]
            m = metrics[str(edition)]
            yesno_acc = float("nan")
            if args.yesno_metrics_json.exists():
                yesno_acc = read_json(args.yesno_metrics_json).get("metrics", {}).get("YesNo_Acc", float("nan"))
            print(
                "BioASQ unified summary "
                f"(edition {edition}): "
                f"YesNo_Acc={yesno_acc:.4f}, "
                f"Factoid_Strict={m['Factoid_Strict_Acc']:.4f}, "
                f"List_P={m['List_Precision']:.4f}, "
                f"List_R={m['List_Recall']:.4f}, "
                f"List_F1={m['List_F1']:.4f}",
                flush=True,
            )
            print(f"Unified metrics JSON: {args.unified_metrics_json}", flush=True)
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(_cli_main())
