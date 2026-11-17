#!/usr/bin/env python3
"""Run Biomix evaluation from the BiOmics-BRICK skill layer.

This wrapper evaluates existing Biomix answer files with the benchmark's
official/local evaluator. It does not generate predictions by default; use the
benchmark runner for generation, then use this command to score the output.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


def _default_benchmark_root() -> Path:
    candidates = [
        Path.cwd() / "biomics-benchmark",
        Path.cwd(),
        Path(r"D:\code\Biomics\biomics-benchmark"),
    ]
    for candidate in candidates:
        if (candidate / "evaluate_pipeline" / "biomix" / "evaluate_results.py").exists():
            return candidate.resolve()
    return candidates[-1].resolve()


def _default_method(root: Path) -> str:
    answer_root = root / "agent_answer" / "biomix"
    for name in (
        "deepseek_v4_pro_D",
        "deepseek_v4_pro_D_main_20260706",
        "hermes_biomics",
    ):
        if (answer_root / name).exists():
            return name
    return "deepseek_v4_pro_D"


def read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as f:
        return json.load(f)


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _run(cmd: list[str], cwd: Path, *, dry_run: bool = False) -> subprocess.CompletedProcess[str] | None:
    print("RUN " + " ".join(str(part) for part in cmd), flush=True)
    if dry_run:
        return None
    return subprocess.run(
        cmd,
        cwd=str(cwd),
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def _metric_rows(method: str, raw_metrics: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for metric in raw_metrics.get(method, []):
        task = str(metric.get("task", ""))
        for key, value in metric.items():
            if key == "task" or not isinstance(value, (int, float)):
                continue
            rows.append(
                {
                    "dataset": "Biomix",
                    "task": task,
                    "metric": key,
                    "method": method,
                    "value": value,
                }
            )
    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["dataset", "task", "metric", "method", "value"]
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_biomix_metrics(args: argparse.Namespace) -> dict[str, Any]:
    eval_script = args.benchmark_root / "evaluate_pipeline" / "biomix" / "evaluate_results.py"
    short_gold = args.short_gold or (
        args.benchmark_root / "evaluate_pipeline" / "biomix" / "gold" / "pure_question_entity_gold.json"
    )
    if not eval_script.exists():
        raise FileNotFoundError(f"Biomix evaluator not found: {eval_script}")
    if not (args.benchmark_root / "agent_answer" / "biomix" / args.method).exists():
        raise FileNotFoundError(
            f"Biomix answer directory not found: {args.benchmark_root / 'agent_answer' / 'biomix' / args.method}"
        )
    if not short_gold.exists():
        raise FileNotFoundError(f"Biomix short-answer gold file not found: {short_gold}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_metrics_path = args.raw_metrics or (args.output_dir / f"{args.method}.raw_metrics.json")
    cmd = [
        sys.executable,
        str(eval_script),
        "--model",
        args.method,
        "--short-gold",
        str(short_gold),
        "--heuristic-short-extraction",
        "--json-output",
        str(raw_metrics_path),
    ]
    completed = _run(cmd, args.benchmark_root, dry_run=args.dry_run)
    if args.dry_run:
        return {}
    assert completed is not None
    log_path = args.output_dir / "logs" / f"biomix_{args.method}.stdout.txt"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text((completed.stdout or "") + "\n" + (completed.stderr or ""), encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"Biomix evaluator failed with exit code {completed.returncode}; see {log_path}")
    if not raw_metrics_path.exists():
        raise FileNotFoundError(f"Biomix evaluator did not write metrics: {raw_metrics_path}")
    return read_json(raw_metrics_path)


def run_openqa_lexical(args: argparse.Namespace) -> dict[str, Any] | None:
    if args.skip_openqa_lexical or args.dry_run:
        return None
    answer_path = args.benchmark_root / "agent_answer" / "biomix" / args.method / "pure_question.json"
    if not answer_path.exists():
        return None
    script = args.benchmark_root / "evaluate_pipeline" / "biomix" / "scripts" / "evaluate_pure_question_entity_gold.py"
    gold = args.short_gold or (
        args.benchmark_root / "evaluate_pipeline" / "biomix" / "gold" / "pure_question_entity_gold.json"
    )
    cmd = [sys.executable, str(script), "--gold", str(gold), "--answers", str(answer_path)]
    completed = _run(cmd, args.benchmark_root, dry_run=False)
    assert completed is not None
    out_path = args.output_dir / f"{args.method}.openqa_lexical.json"
    out_path.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"Biomix OpenQA lexical evaluator failed with exit code {completed.returncode}")
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError:
        return {"parse_error": True, "raw_stdout_path": str(out_path)}


def write_outputs(args: argparse.Namespace, raw_metrics: dict[str, Any], lexical: dict[str, Any] | None) -> None:
    rows = _metric_rows(args.method, raw_metrics)
    if lexical:
        for key in ("entity_match_rate", "mean_concept_recall"):
            value = lexical.get(key)
            if isinstance(value, (int, float)):
                rows.append(
                    {
                        "dataset": "Biomix",
                        "task": "openqa_lexical",
                        "metric": key,
                        "method": args.method,
                        "value": value,
                    }
                )
    payload = {
        "mode": "existing_predictions",
        "dataset": "Biomix",
        "method": args.method,
        "warning": "Gold files are used only by the evaluator after prediction.",
        "benchmark_root": str(args.benchmark_root),
        "raw_metrics": raw_metrics,
        "openqa_lexical": lexical,
        "rows": rows,
    }
    write_json(args.metrics_json, payload)
    _write_csv(args.metrics_csv, rows)
    write_json(
        args.manifest,
        {
            "dataset": "Biomix",
            "method": args.method,
            "benchmark_root": str(args.benchmark_root),
            "outputs": {
                "metrics_json": str(args.metrics_json),
                "metrics_csv": str(args.metrics_csv),
            },
        },
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate existing Biomix D-group/new-skill answers.")
    parser.add_argument("--benchmark-root", type=Path, default=_default_benchmark_root())
    parser.add_argument("--method", default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--raw-metrics", type=Path, default=None)
    parser.add_argument("--metrics-json", type=Path, default=None)
    parser.add_argument("--metrics-csv", type=Path, default=None)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--short-gold", type=Path, default=None)
    parser.add_argument("--skip-openqa-lexical", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    args.benchmark_root = args.benchmark_root.resolve()
    args.method = args.method or _default_method(args.benchmark_root)
    args.output_dir = (args.output_dir or (args.benchmark_root / "evaluate_pipeline" / "skill_benchmark_eval" / "biomix")).resolve()
    args.metrics_json = args.metrics_json or (args.output_dir / f"{args.method}.metrics.json")
    args.metrics_csv = args.metrics_csv or (args.output_dir / f"{args.method}.metrics.csv")
    args.manifest = args.manifest or (args.output_dir / f"{args.method}.manifest.json")
    return args


def _cli_main() -> int:
    try:
        args = parse_args()
        raw_metrics = run_biomix_metrics(args)
        lexical = run_openqa_lexical(args)
        if not args.dry_run:
            write_outputs(args, raw_metrics, lexical)
            rows = _metric_rows(args.method, raw_metrics)
            key = {(row["task"], row["metric"]): row["value"] for row in rows}
            print(
                "Biomix summary: "
                f"true_false_acc={key.get(('true_false', 'accuracy'), float('nan')):.4f}, "
                f"mcq_acc={key.get(('mcq_5_choice', 'accuracy'), float('nan')):.4f}",
                flush=True,
            )
            print(f"Metrics JSON: {args.metrics_json}", flush=True)
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(_cli_main())
