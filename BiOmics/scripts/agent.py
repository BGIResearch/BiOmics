"""brick agent - deterministic orchestration layer for BiOmics-BRICK."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from agent_checker import (
    CheckIssue,
    ValidatorResult,
    has_errors,
    has_failed,
    load_json,
    run_artifact_validators,
    validate_agent_plan,
)
from agent_runner import run_plan


def _write_json(data: Any, output: str | None = None) -> None:
    text = json.dumps(data, indent=2, ensure_ascii=False)
    if output:
        Path(output).write_text(text + "\n", encoding="utf-8")
        print(f"Wrote {output}")
    else:
        print(text)


def _split_csv(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def _print_issues(issues: list[CheckIssue]) -> None:
    if not issues:
        print("Agent plan is valid.")
        return
    for issue in issues:
        prefix = issue.level.upper()
        target = f"[{issue.step_id}] " if issue.step_id else ""
        validator = f" ({issue.validator})" if issue.validator else ""
        print(f"{prefix}: {target}{issue.message}{validator}")


def _print_results(results: list[ValidatorResult]) -> None:
    if not results:
        print("No artifact validators ran.")
        return
    for result in results:
        status = result.status.upper()
        print(f"{status}: {result.validator}: {result.message}")


def cmd_check_plan(args: argparse.Namespace) -> int:
    plan_path = Path(args.plan)
    plan = load_json(plan_path)
    issues = validate_agent_plan(plan, base_dir=plan_path.parent, workdir=args.workdir)
    payload = {"valid": not has_errors(issues), "issues": [issue.as_dict() for issue in issues]}
    if args.format == "json":
        _write_json(payload, args.output)
    else:
        _print_issues(issues)
    return 1 if has_errors(issues) else 0


def cmd_check(args: argparse.Namespace) -> int:
    validators = _split_csv(args.validators)
    validator_args: dict[str, Any] = {}
    if args.required_columns:
        validator_args["csv_columns_present"] = {"required_columns": _split_csv(args.required_columns)}
    if args.max_hops:
        validator_args["max_hops_respected"] = {"max_hops": args.max_hops}
    if args.validator_args:
        validator_args.update(json.loads(args.validator_args))
    results = run_artifact_validators(args.input, validators, validator_args=validator_args)
    payload = {"valid": not has_failed(results), "results": [result.as_dict() for result in results]}
    if args.format == "json":
        _write_json(payload, args.output)
    else:
        _print_results(results)
    return 1 if has_failed(results) else 0


def cmd_run(args: argparse.Namespace) -> int:
    summary = run_plan(
        args.plan,
        workdir=args.workdir,
        execute=args.execute,
        max_retries=args.max_retries,
        manifest=args.manifest,
    )
    if args.format == "json":
        _write_json(summary)
    else:
        status = "OK" if summary.get("ok") else "FAILED"
        print(f"Agent run {status}. Manifest: {summary.get('manifest')}")
        for step in summary.get("steps", []):
            print(f"- {step['step_id']}: {step['status']}")
    return 0 if summary.get("ok") else 1


def _read_manifest(path: Path) -> list[dict[str, Any]]:
    rows = []
    if not path.exists():
        return rows
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append(json.loads(line))
    return rows


def _evidence_strength_summary(evidence_path: Path | None) -> list[str]:
    if not evidence_path or not evidence_path.exists():
        return ["- Evidence table: not provided"]
    df = pd.read_csv(evidence_path)
    lines = [f"- Evidence rows: {len(df)}"]
    if "evidence_strength" in df.columns:
        counts = df["evidence_strength"].fillna("missing").astype(str).value_counts().to_dict()
        lines.append("- Evidence strength: " + ", ".join(f"{key}={value}" for key, value in counts.items()))
    if "manual_review_required" in df.columns:
        review = df["manual_review_required"].fillna(False).astype(str).str.lower().isin({"true", "1", "yes"}).sum()
        lines.append(f"- Manual review required rows: {int(review)}")
    if "edge_review_note" in df.columns:
        notes = df["edge_review_note"].dropna().astype(str)
        notes = [note for note in notes.unique().tolist() if note.strip()][:5]
        if notes:
            lines.append("- Review notes:")
            lines.extend(f"  - {note}" for note in notes)
    return lines


def cmd_report(args: argparse.Namespace) -> int:
    manifest_path = Path(args.manifest)
    events = _read_manifest(manifest_path)
    event_counts = Counter(event.get("event") for event in events)
    failed_validators = [
        event
        for event in events
        if event.get("event") == "validator_finished" and event.get("status") == "failed"
    ]
    step_status = [
        event
        for event in events
        if event.get("event") in {"step_success", "step_failed_final"}
    ]
    evidence_path = Path(args.evidence) if args.evidence else None
    lines = [
        "# BRICK Agent Run Report",
        "",
        "## Task Summary",
        f"- Manifest: `{manifest_path}`",
        f"- Events: {len(events)}",
        "",
        "## Execution Summary",
    ]
    if step_status:
        for event in step_status:
            status = "success" if event.get("event") == "step_success" else "failed"
            lines.append(f"- {event.get('step_id')}: {status}")
    else:
        lines.append("- No completed step events found")
    lines.extend(
        [
            "",
            "## Validator Summary",
            f"- Validator failures: {len(failed_validators)}",
        ]
    )
    for event in failed_validators[:10]:
        lines.append(f"- {event.get('step_id')} / {event.get('validator')}: {event.get('message')}")
    lines.extend(["", "## Evidence Strength Summary", *_evidence_strength_summary(evidence_path)])
    lines.extend(
        [
            "",
            "## Event Counts",
            *[f"- {name}: {count}" for name, count in sorted(event_counts.items())],
        ]
    )
    output = "\n".join(lines) + "\n"
    if args.output:
        Path(args.output).write_text(output, encoding="utf-8")
        print(f"Wrote {args.output}")
    else:
        print(output)
    return 1 if failed_validators else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="brick agent",
        description="Check, dry-run, execute, and report deterministic brick agent plans.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("check-plan", help="Validate an agent plan before execution.")
    p.add_argument("-p", "--plan", required=True)
    p.add_argument("--workdir", help="Override workdir used for output boundary checks.")
    p.add_argument("--format", choices=["table", "json"], default="table")
    p.add_argument("-o", "--output")
    p.set_defaults(func=cmd_check_plan)

    p = sub.add_parser("check", help="Run validators against an output artifact.")
    p.add_argument("-i", "--input", required=True)
    p.add_argument("--validators", required=True, help="Comma-separated validator names.")
    p.add_argument("--required-columns", help="Comma-separated columns for csv_columns_present.")
    p.add_argument("--max-hops", type=int, default=3)
    p.add_argument("--validator-args", help="JSON object keyed by validator name.")
    p.add_argument("--format", choices=["table", "json"], default="table")
    p.add_argument("-o", "--output")
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("run", help="Dry-run or execute an agent plan.")
    p.add_argument("-p", "--plan", required=True)
    p.add_argument("--workdir", help="Run workdir and manifest directory.")
    p.add_argument("--execute", action="store_true", help="Actually execute safe brick commands. Default is dry-run.")
    p.add_argument("--max-retries", type=int, help="Override per-step retry count.")
    p.add_argument("--manifest", help="Manifest JSONL output path.")
    p.add_argument("--format", choices=["table", "json"], default="table")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("report", help="Create a lightweight markdown report from a run manifest.")
    p.add_argument("--manifest", required=True)
    p.add_argument("--evidence", help="Optional scored/normalized evidence CSV.")
    p.add_argument("-o", "--output")
    p.set_defaults(func=cmd_report)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args) or 0)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


def _cli_main() -> int:
    return main()


if __name__ == "__main__":
    sys.exit(main())
