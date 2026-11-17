#!/usr/bin/env python3
"""Evaluate BioHopR, BioASQ, and Biomix D-group outputs from one skill command."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_ROOT_FROM_SCRIPT = SCRIPT_DIR.parents[2]

KEY_METRIC_ORDER = [
    ("Biomix", "true_false", "accuracy", None),
    ("Biomix", "true_false", "macro_f1", None),
    ("Biomix", "mcq_5_choice", "accuracy", None),
    ("Biomix", "openqa_lexical", "entity_match_rate", None),
    ("Biomix", "openqa_lexical", "mean_concept_recall", None),
    ("BioHopR", "1hop_single", "hit_rate", None),
    ("BioHopR", "1hop_single", "macro_f1", None),
    ("BioHopR", "2hop_single", "hit_rate", None),
    ("BioHopR", "2hop_single", "macro_f1", None),
    ("BioASQ", "YesNo", "YesNo_Acc", None),
    ("BioASQ", "YesNo", "YesNo_macroF1", None),
    ("BioASQ", "Factoid", "Factoid_Strict_Acc", "8"),
    ("BioASQ", "List", "List_Precision", "8"),
    ("BioASQ", "List", "List_Recall", "8"),
    ("BioASQ", "List", "List_F1", "8"),
]

METRIC_DISPLAY_NAMES = {
    "accuracy": "Accuracy",
    "macro_f1": "Macro-F1",
    "entity_match_rate": "Entity hit rate",
    "mean_concept_recall": "Mean concept recall",
    "hit_rate": "Hit rate",
    "YesNo_Acc": "Accuracy",
    "YesNo_macroF1": "Macro-F1",
    "Factoid_Strict_Acc": "Strict accuracy",
    "List_Precision": "Precision",
    "List_Recall": "Recall",
    "List_F1": "F1",
}

def _default_benchmark_root() -> Path:
    candidates = [
        Path.cwd() / "biomics-benchmark",
        Path.cwd(),
        Path(r"D:\code\Biomics\biomics-benchmark"),
    ]
    for candidate in candidates:
        if (candidate / "evaluate_pipeline" / "BioHopR" / "evaluate_results.py").exists():
            return candidate.resolve()
    return candidates[-1].resolve()


def _default_bioasq_root() -> Path:
    candidates = [
        BENCHMARK_ROOT_FROM_SCRIPT / "evaluate_pipeline" / "BioASQ",
        Path.cwd() / "BioASQ Benchmark",
        Path(r"D:\code\Biomics\BioASQ Benchmark"),
    ]
    for candidate in candidates:
        if (candidate / "scripts" / "extract_and_standardize.py").exists():
            return candidate.resolve()
    return candidates[-1].resolve()


def _default_bioasq_input(root: Path) -> Path:
    for candidate in (
        BENCHMARK_ROOT_FROM_SCRIPT / "agent_answer" / "BioASQ" / "deepseek_v4_pro_D" / "deepseek_v4_pro_D.results.json",
        root / "results_old" / "deepseek_v4_pro_D.results.json",
        root / "results_old" / "deepseek_v4_pro_D_biomics_githubinstall_latestkg_20260706.results.json",
    ):
        if candidate.exists():
            return candidate.resolve()
    return (
        BENCHMARK_ROOT_FROM_SCRIPT
        / "agent_answer"
        / "BioASQ"
        / "deepseek_v4_pro_D"
        / "deepseek_v4_pro_D.results.json"
    ).resolve()


def _default_bioasq_gold(benchmark_root: Path) -> Path:
    candidates = [
        benchmark_root / "evaluate_pipeline" / "BioASQ" / "gold" / "answer.json",
        BENCHMARK_ROOT_FROM_SCRIPT / "evaluate_pipeline" / "BioASQ" / "gold" / "answer.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def _bioasq_existing_legacy_outputs(input_path: Path) -> tuple[Path | None, Path | None]:
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


def read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as f:
        return json.load(f)


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def display_path(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve())).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")


def _run(cmd: list[str], cwd: Path, *, dry_run: bool = False) -> subprocess.CompletedProcess[str] | None:
    print("RUN " + " ".join(str(part) for part in cmd), flush=True)
    if dry_run:
        return None
    completed = subprocess.run(
        cmd,
        cwd=str(cwd),
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.stdout:
        print(completed.stdout.strip(), flush=True)
    if completed.stderr:
        print(completed.stderr.strip(), file=sys.stderr, flush=True)
    return completed


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    preferred = ["dataset", "task", "metric", "method", "edition", "value", "status", "path"]
    fields = [field for field in preferred if field in fields] + [field for field in fields if field not in preferred]
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_biomix(args: argparse.Namespace) -> list[dict[str, Any]]:
    out_dir = args.output_dir / "biomix"
    metrics_json = out_dir / f"{args.biomix_method}.metrics.json"
    cmd = [
        sys.executable,
        str(SCRIPT_DIR / "biomix_eval.py"),
        "--benchmark-root",
        str(args.benchmark_root),
        "--method",
        args.biomix_method,
        "--output-dir",
        str(out_dir),
        "--metrics-json",
        str(metrics_json),
    ]
    if args.dry_run:
        cmd.append("--dry-run")
    completed = _run(cmd, args.benchmark_root, dry_run=False)
    if args.dry_run:
        return []
    assert completed is not None
    if completed.returncode != 0:
        raise RuntimeError("Biomix evaluation failed")
    payload = read_json(metrics_json)
    return list(payload.get("rows", []))


def run_biohopr(args: argparse.Namespace) -> list[dict[str, Any]]:
    out_dir = args.output_dir / "biohopr"
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics_json = out_dir / f"{args.biohopr_method}.metrics.json"
    evaluator = args.benchmark_root / "evaluate_pipeline" / "BioHopR" / "evaluate_results.py"
    cmd = [
        sys.executable,
        str(evaluator),
        "--model",
        args.biohopr_method,
        "--json-output",
        str(metrics_json),
    ]
    completed = _run(cmd, args.benchmark_root, dry_run=args.dry_run)
    if args.dry_run:
        return []
    assert completed is not None
    if completed.returncode != 0:
        raise RuntimeError("BioHopR evaluation failed")
    payload = read_json(metrics_json)
    rows: list[dict[str, Any]] = []
    for method, tasks in payload.items():
        for task, metric_map in tasks.items():
            for metric, value in metric_map.items():
                if metric.endswith("_file") or not isinstance(value, (int, float)):
                    continue
                rows.append(
                    {
                        "dataset": "BioHopR",
                        "task": task,
                        "metric": metric,
                        "method": method,
                        "value": value,
                    }
                )
    return rows


def run_bioasq(args: argparse.Namespace) -> list[dict[str, Any]]:
    out_dir = args.output_dir / "bioasq"
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics_json = out_dir / "bioasq_legacy_goldaware.metrics.json"
    metrics_csv = out_dir / "bioasq_legacy_goldaware.metrics.csv"
    manifest = out_dir / "bioasq_legacy_goldaware.manifest.json"
    unified_eval_output = out_dir / "bioasq_unified.eval.json"
    unified_metrics_json = out_dir / "bioasq_unified.metrics.json"
    unified_metrics_csv = out_dir / "bioasq_unified.metrics.csv"
    unified_manifest = out_dir / "bioasq_unified.manifest.json"
    unified_official_metrics_json = out_dir / "bioasq_unified.official_metrics.json"
    cmd = [
        sys.executable,
        str(SCRIPT_DIR / "bioasq_eval.py"),
        "--bioasq-root",
        str(args.bioasq_root),
        "--input",
        str(args.bioasq_input),
        "--gold",
        str(args.bioasq_gold),
        "--metrics-json",
        str(metrics_json),
        "--metrics-csv",
        str(metrics_csv),
        "--manifest",
        str(manifest),
        "--yesno-metrics-json",
        str(args.bioasq_yesno_metrics),
        "--yesno-eval-json",
        str(args.bioasq_yesno_eval),
        "--unified-eval-output",
        str(unified_eval_output),
        "--unified-metrics-json",
        str(unified_metrics_json),
        "--unified-metrics-csv",
        str(unified_metrics_csv),
        "--unified-manifest",
        str(unified_manifest),
        "--unified-official-metrics-json",
        str(unified_official_metrics_json),
        "--skip-standardizer-if-exists",
    ]
    standardized, eval_output = _bioasq_existing_legacy_outputs(args.bioasq_input)
    if standardized and eval_output:
        cmd.extend(
            [
                "--standardized-output",
                str(standardized),
                "--eval-output",
                str(eval_output),
                "--skip-standardizer",
            ]
        )
    completed = _run(cmd, args.bioasq_root, dry_run=args.dry_run)
    if args.dry_run:
        return []
    assert completed is not None
    if completed.returncode != 0:
        raise RuntimeError("BioASQ evaluation failed")
    payload = read_json(unified_metrics_json)
    return list(payload.get("rows", []))


def key_metric_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    wanted = {key: index for index, key in enumerate(KEY_METRIC_ORDER)}
    key_rows = []
    for row in rows:
        key = (row.get("dataset"), row.get("task"), row.get("metric"), row.get("edition"))
        key_no_edition = (row.get("dataset"), row.get("task"), row.get("metric"), None)
        order = wanted.get(key, wanted.get(key_no_edition))
        if order is not None:
            ordered = dict(row)
            ordered["_order"] = order
            key_rows.append(ordered)
    key_rows.sort(key=lambda row: row.get("_order", 10_000))
    for row in key_rows:
        row.pop("_order", None)
    return key_rows


def write_key_metrics(args: argparse.Namespace, rows: list[dict[str, Any]]) -> None:
    key_rows = key_metric_rows(rows)
    write_json(args.output_dir / "benchmark_suite_key_metrics.json", key_rows)
    _write_csv(args.output_dir / "benchmark_suite_key_metrics.csv", key_rows)


def _rows_for_task(rows: list[dict[str, Any]], dataset: str, task: str) -> list[dict[str, Any]]:
    return [row for row in rows if row.get("dataset") == dataset and row.get("task") == task]


def _headline_entry(row: dict[str, Any]) -> dict[str, Any]:
    metric = str(row.get("metric"))
    entry = dict(row)
    entry["display_metric"] = METRIC_DISPLAY_NAMES.get(metric, metric)
    return entry


def build_io_manifest(args: argparse.Namespace, rows: list[dict[str, Any]]) -> dict[str, Any]:
    key_rows = key_metric_rows(rows)
    root = args.benchmark_root
    out = args.output_dir
    biomix_method = args.biomix_method
    biohopr_method = args.biohopr_method
    bioasq_out = out / "bioasq"

    return {
        "method": "deepseek_v4_pro_D",
        "bundle": "D_main_hopaware",
        "purpose": (
            "Maps the headline metrics table to the exact question inputs, prediction outputs, "
            "gold/evaluator inputs, and metric outputs used to produce it."
        ),
        "headline_metrics": [_headline_entry(row) for row in key_rows],
        "headline_metrics_files": {
            "json": display_path(out / "benchmark_suite_key_metrics.json", root),
            "csv": display_path(out / "benchmark_suite_key_metrics.csv", root),
            "long_json": display_path(out / "benchmark_suite_metrics_long.json", root),
            "long_csv": display_path(out / "benchmark_suite_metrics_long.csv", root),
        },
        "datasets": {
            "Biomix": {
                "method": biomix_method,
                "protocol": "Existing D predictions are scored after generation; gold labels are evaluator-only.",
                "tasks": {
                    "true_false": {
                        "question_input": display_path(root / "dataset" / "biomix" / "true_false_biomix_pure_question.json", root),
                        "gold_input": display_path(root / "evaluate_pipeline" / "biomix" / "true_false_biomix.csv", root),
                        "prediction_output": display_path(root / "agent_answer" / "biomix" / biomix_method / "true_false_biomix_pure_question.json", root),
                        "evaluator": display_path(root / "evaluate_pipeline" / "biomix" / "evaluate_results.py", root),
                        "suite_wrapper": display_path(SCRIPT_DIR / "biomix_eval.py", root),
                        "metric_outputs": [
                            display_path(out / "biomix" / f"{biomix_method}.metrics.json", root),
                            display_path(out / "biomix" / f"{biomix_method}.raw_metrics.json", root),
                        ],
                        "headline_rows": _rows_for_task(key_rows, "Biomix", "true_false"),
                    },
                    "mcq_5_choice": {
                        "question_input": display_path(root / "dataset" / "biomix" / "mcq_biomix_pure_question.json", root),
                        "gold_input": display_path(root / "evaluate_pipeline" / "biomix" / "mcq_biomix.csv", root),
                        "prediction_output": display_path(root / "agent_answer" / "biomix" / biomix_method / "mcq_biomix_pure_question.json", root),
                        "evaluator": display_path(root / "evaluate_pipeline" / "biomix" / "evaluate_results.py", root),
                        "suite_wrapper": display_path(SCRIPT_DIR / "biomix_eval.py", root),
                        "metric_outputs": [
                            display_path(out / "biomix" / f"{biomix_method}.metrics.json", root),
                            display_path(out / "biomix" / f"{biomix_method}.raw_metrics.json", root),
                        ],
                        "headline_rows": _rows_for_task(key_rows, "Biomix", "mcq_5_choice"),
                    },
                    "openqa_lexical": {
                        "question_input": display_path(root / "dataset" / "biomix" / "pure_question.json", root),
                        "gold_input": display_path(root / "evaluate_pipeline" / "biomix" / "gold" / "pure_question_entity_gold.json", root),
                        "prediction_output": display_path(root / "agent_answer" / "biomix" / biomix_method / "pure_question.json", root),
                        "evaluator": display_path(root / "evaluate_pipeline" / "biomix" / "scripts" / "evaluate_pure_question_entity_gold.py", root),
                        "suite_wrapper": display_path(SCRIPT_DIR / "biomix_eval.py", root),
                        "metric_outputs": [
                            display_path(out / "biomix" / f"{biomix_method}.metrics.json", root),
                            display_path(out / "biomix" / f"{biomix_method}.openqa_lexical.json", root),
                        ],
                        "headline_rows": _rows_for_task(key_rows, "Biomix", "openqa_lexical"),
                    },
                },
            },
            "BioHopR": {
                "method": biohopr_method,
                "protocol": (
                    "Strict question-only D evaluation. Inference uses question-derived extraction, "
                    "not gold_hops, official bridge entities, official query/target types, or gold answers."
                ),
                "tasks": {
                    "1hop_single": {
                        "question_input": display_path(root / "dataset" / "BioHopR" / "questions_1hop_single.json", root),
                        "gold_input": display_path(root / "evaluate_pipeline" / "BioHopR" / "gold" / "questions_1hop_single_gold.json", root),
                        "prediction_output": display_path(root / "agent_answer" / "BioHopR" / biohopr_method / "1hop_answers.json", root),
                        "evaluator": display_path(root / "evaluate_pipeline" / "BioHopR" / "evaluate_results.py", root),
                        "metric_outputs": [
                            display_path(out / "biohopr" / f"{biohopr_method}.metrics.json", root),
                        ],
                        "headline_rows": _rows_for_task(key_rows, "BioHopR", "1hop_single"),
                    },
                    "2hop_single": {
                        "question_input": display_path(root / "dataset" / "BioHopR" / "questions_2hop_single.json", root),
                        "gold_input": display_path(root / "evaluate_pipeline" / "BioHopR" / "gold" / "questions_2hop_single_gold.json", root),
                        "prediction_output": display_path(root / "agent_answer" / "BioHopR" / biohopr_method / "2hop_answers.json", root),
                        "evaluator": display_path(root / "evaluate_pipeline" / "BioHopR" / "evaluate_results.py", root),
                        "metric_outputs": [
                            display_path(out / "biohopr" / f"{biohopr_method}.metrics.json", root),
                        ],
                        "headline_rows": _rows_for_task(key_rows, "BioHopR", "2hop_single"),
                    },
                },
            },
            "BioASQ": {
                "method": "deepseek_v4_pro_D",
                "protocol": (
                    "One BioASQ bundle with task-level protocols. YesNo uses the strict/no-gold "
                    "prior-guided self-correction artifact. Factoid/List use the uploaded legacy gold-aware "
                    "standardized evaluator input and must not be called strict no-gold."
                ),
                "shared_inputs": {
                    "question_input": display_path(root / "dataset" / "BioASQ" / "questions.json", root),
                    "gold_input": display_path(args.bioasq_gold, root),
                    "unified_official_eval_input": display_path(bioasq_out / "bioasq_unified.eval.json", root),
                    "unified_metrics": display_path(bioasq_out / "bioasq_unified.metrics.json", root),
                    "unified_official_java_metrics": display_path(bioasq_out / "bioasq_unified.official_metrics.json", root),
                },
                "tasks": {
                    "YesNo": {
                        "question_input": display_path(root / "dataset" / "BioASQ" / "questions.json", root),
                        "prior_input": display_path(args.bioasq_input, root),
                        "prediction_output": display_path(bioasq_out / "bioasq_yesno_prior_guarded.predictions.json", root),
                        "yesno_eval_input": display_path(bioasq_out / "bioasq_yesno_prior_guarded.eval.json", root),
                        "unified_official_eval_input": display_path(bioasq_out / "bioasq_unified.eval.json", root),
                        "evaluator": display_path(args.bioasq_root / "Evaluation-Measures" / "flat" / "BioASQEvaluation" / "dist" / "BioASQEvaluation.jar", root),
                        "suite_wrapper": display_path(SCRIPT_DIR / "bioasq_eval.py", root),
                        "metric_outputs": [
                            display_path(bioasq_out / "bioasq_yesno_prior_guarded.official_metrics.json", root),
                            display_path(bioasq_out / "bioasq_unified.metrics.json", root),
                            display_path(bioasq_out / "bioasq_unified.official_metrics.json", root),
                        ],
                        "leakage_assessment": (
                            "No benchmark gold is used during YesNo generation; gold is used by the evaluator "
                            "after prediction. This is prior-guided self-correction because the previous raw "
                            "answer is used as a prior, and it is not a held-out claim."
                        ),
                        "headline_rows": _rows_for_task(key_rows, "BioASQ", "YesNo"),
                    },
                    "Factoid": {
                        "raw_prediction_input": display_path(args.bioasq_input, root),
                        "legacy_standardized_eval_input": display_path(root / "agent_answer" / "BioASQ" / "deepseek_v4_pro_D" / "deepseek_v4_pro_D.legacy_goldaware_standardized.eval.json", root),
                        "unified_official_eval_input": display_path(bioasq_out / "bioasq_unified.eval.json", root),
                        "evaluator": display_path(args.bioasq_root / "Evaluation-Measures" / "flat" / "BioASQEvaluation" / "dist" / "BioASQEvaluation.jar", root),
                        "suite_wrapper": display_path(SCRIPT_DIR / "bioasq_eval.py", root),
                        "metric_outputs": [
                            display_path(bioasq_out / "bioasq_legacy_goldaware.metrics.json", root),
                            display_path(bioasq_out / "bioasq_unified.metrics.json", root),
                            display_path(bioasq_out / "bioasq_unified.official_metrics.json", root),
                        ],
                        "leakage_assessment": "Legacy gold-aware standardization uses gold exact answers; keep for historical comparison only.",
                        "headline_rows": _rows_for_task(key_rows, "BioASQ", "Factoid"),
                    },
                    "List": {
                        "raw_prediction_input": display_path(args.bioasq_input, root),
                        "legacy_standardized_eval_input": display_path(root / "agent_answer" / "BioASQ" / "deepseek_v4_pro_D" / "deepseek_v4_pro_D.legacy_goldaware_standardized.eval.json", root),
                        "unified_official_eval_input": display_path(bioasq_out / "bioasq_unified.eval.json", root),
                        "evaluator": display_path(args.bioasq_root / "Evaluation-Measures" / "flat" / "BioASQEvaluation" / "dist" / "BioASQEvaluation.jar", root),
                        "suite_wrapper": display_path(SCRIPT_DIR / "bioasq_eval.py", root),
                        "metric_outputs": [
                            display_path(bioasq_out / "bioasq_legacy_goldaware.metrics.json", root),
                            display_path(bioasq_out / "bioasq_unified.metrics.json", root),
                            display_path(bioasq_out / "bioasq_unified.official_metrics.json", root),
                        ],
                        "leakage_assessment": "Legacy gold-aware standardization uses gold exact answers; keep for historical comparison only.",
                        "headline_rows": _rows_for_task(key_rows, "BioASQ", "List"),
                    },
                },
            },
        },
        "supplemental_external_source": {
            "Biomix_GWAS_MCQ": {
                "method": "deepseek_v4_pro_D_gwas_sourceaware_20260707",
                "status": "supplemental_external_source_not_D_main",
                "question_input": display_path(root / "dataset" / "biomix" / "mcq_biomix_pure_question.json", root),
                "prediction_output": display_path(root / "agent_answer" / "biomix" / "deepseek_v4_pro_D_gwas_sourceaware_20260707" / "mcq_biomix_pure_question.json", root),
                "metric_output": display_path(root / "evaluate_pipeline" / "skill_benchmark_eval" / "biomix" / "deepseek_v4_pro_D_gwas_sourceaware_20260707.metrics.json", root),
                "note": "No benchmark-gold leakage, but this uses GWAS Catalog as an external source-aware verifier and is not the D-main KG-only result.",
            }
        },
    }


def _format_value(value: Any) -> str:
    if isinstance(value, (int, float)):
        return f"{value:.4f}"
    return str(value)


def write_io_markdown(args: argparse.Namespace, manifest: dict[str, Any]) -> None:
    lines = [
        "# D Main Test Inputs and Outputs",
        "",
        "This file is generated by `skills/biomics_brick_fix/scripts/benchmark_suite_eval.py`.",
        "It maps each headline metric to the concrete input, prediction, evaluator, and metric files used to produce it.",
        "",
        "## Headline Metrics",
        "",
        "| Dataset | Task | Metric | Value |",
        "| --- | --- | --- | ---: |",
    ]
    for row in manifest["headline_metrics"]:
        lines.append(
            f"| {row['dataset']} | {row['task']} | {row.get('display_metric', row['metric'])} | {_format_value(row['value'])} |"
        )
    lines.extend(
        [
            "",
            "## Dataset IO Map",
            "",
        ]
    )
    for dataset, dataset_payload in manifest["datasets"].items():
        lines.extend([f"### {dataset}", "", dataset_payload["protocol"], ""])
        for task, task_payload in dataset_payload["tasks"].items():
            lines.extend([f"#### {task}", ""])
            for key in (
                "question_input",
                "prior_input",
                "raw_prediction_input",
                "gold_input",
                "prediction_output",
                "yesno_eval_input",
                "legacy_standardized_eval_input",
                "unified_official_eval_input",
                "evaluator",
                "suite_wrapper",
            ):
                value = task_payload.get(key)
                if value:
                    lines.append(f"- `{key}`: `{value}`")
            metric_outputs = task_payload.get("metric_outputs", [])
            if metric_outputs:
                lines.append("- `metric_outputs`: " + ", ".join(f"`{path}`" for path in metric_outputs))
            leakage = task_payload.get("leakage_assessment")
            if leakage:
                lines.append(f"- `leakage_assessment`: {leakage}")
            headline = task_payload.get("headline_rows", [])
            if headline:
                values = ", ".join(
                    f"{METRIC_DISPLAY_NAMES.get(str(row['metric']), row['metric'])}={_format_value(row['value'])}"
                    for row in headline
                )
                lines.append(f"- `headline_values`: {values}")
            lines.append("")
    lines.extend(
        [
            "## Supplemental External Source",
            "",
            "The GWAS source-aware Biomix MCQ run is uploaded for comparison but is not part of the D-main KG-only table.",
            "",
        ]
    )
    for name, payload in manifest["supplemental_external_source"].items():
        lines.extend(
            [
                f"### {name}",
                "",
                f"- `method`: `{payload['method']}`",
                f"- `status`: `{payload['status']}`",
                f"- `question_input`: `{payload['question_input']}`",
                f"- `prediction_output`: `{payload['prediction_output']}`",
                f"- `metric_output`: `{payload['metric_output']}`",
                f"- `note`: {payload['note']}",
                "",
            ]
        )
    (args.output_dir / "TEST_INPUT_OUTPUTS.md").write_text("\n".join(lines), encoding="utf-8")


def write_io_manifest(args: argparse.Namespace, rows: list[dict[str, Any]]) -> None:
    manifest = build_io_manifest(args, rows)
    write_json(args.output_dir / "benchmark_suite_io_manifest.json", manifest)
    write_io_markdown(args, manifest)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate D-group BioHopR, BioASQ, and Biomix outputs.")
    parser.add_argument("--benchmark-root", type=Path, default=_default_benchmark_root())
    parser.add_argument("--bioasq-root", type=Path, default=_default_bioasq_root())
    parser.add_argument("--datasets", nargs="+", choices=["all", "biomix", "biohopr", "bioasq"], default=["all"])
    parser.add_argument("--biomix-method", default="deepseek_v4_pro_D")
    parser.add_argument("--biohopr-method", default="deepseek_v4_pro_D")
    parser.add_argument("--bioasq-input", type=Path, default=None)
    parser.add_argument("--bioasq-gold", type=Path, default=None)
    parser.add_argument("--bioasq-yesno-metrics", type=Path, default=None)
    parser.add_argument("--bioasq-yesno-eval", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    args.benchmark_root = args.benchmark_root.resolve()
    args.bioasq_root = args.bioasq_root.resolve()
    args.bioasq_input = (args.bioasq_input or _default_bioasq_input(args.bioasq_root)).resolve()
    args.bioasq_gold = (args.bioasq_gold or _default_bioasq_gold(args.benchmark_root)).resolve()
    args.output_dir = (
        args.output_dir
        or (args.benchmark_root / "evaluate_pipeline" / "skill_benchmark_eval" / "D_main_hopaware")
    ).resolve()
    args.bioasq_yesno_metrics = (
        args.bioasq_yesno_metrics
        or (
            args.output_dir
            / "bioasq"
            / "bioasq_yesno_prior_guarded.official_metrics.json"
        )
    ).resolve()
    args.bioasq_yesno_eval = (
        args.bioasq_yesno_eval
        or (
            args.output_dir
            / "bioasq"
            / "bioasq_yesno_prior_guarded.eval.json"
        )
    ).resolve()
    return args


def _cli_main() -> int:
    try:
        args = parse_args()
        datasets = {"biomix", "biohopr", "bioasq"} if "all" in args.datasets else set(args.datasets)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        rows: list[dict[str, Any]] = []
        if "biomix" in datasets:
            rows.extend(run_biomix(args))
        if "biohopr" in datasets:
            rows.extend(run_biohopr(args))
        if "bioasq" in datasets:
            rows.extend(run_bioasq(args))
        if not args.dry_run:
            write_json(args.output_dir / "benchmark_suite_metrics_long.json", rows)
            _write_csv(args.output_dir / "benchmark_suite_metrics_long.csv", rows)
            write_key_metrics(args, rows)
            write_io_manifest(args, rows)
            print(f"Benchmark suite metrics: {args.output_dir / 'benchmark_suite_key_metrics.csv'}", flush=True)
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(_cli_main())
