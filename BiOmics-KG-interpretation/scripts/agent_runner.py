"""Deterministic runner for brick agent plans."""

from __future__ import annotations

import json
import subprocess
import time
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from agent_checker import (
    ARTIFACT_VALIDATORS,
    command_to_display,
    has_errors,
    has_failed,
    is_safe_brick_command,
    load_json,
    normalize_command,
    run_artifact_validators,
    validate_agent_plan,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def write_event(handle, event: str, **data: Any) -> None:
    payload = {"time": utc_now(), "event": event}
    payload.update(data)
    handle.write(json.dumps(payload, ensure_ascii=False) + "\n")
    handle.flush()


def step_order(steps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_id = {str(step["id"]): step for step in steps}
    graph: dict[str, list[str]] = defaultdict(list)
    incoming = {step_id: 0 for step_id in by_id}
    for step in steps:
        step_id = str(step["id"])
        for dep in step.get("depends_on", []) or []:
            dep = str(dep)
            graph[dep].append(step_id)
            incoming[step_id] += 1

    queue = deque([step_id for step_id in by_id if incoming[step_id] == 0])
    ordered = []
    while queue:
        step_id = queue.popleft()
        ordered.append(by_id[step_id])
        for child in graph[step_id]:
            incoming[child] -= 1
            if incoming[child] == 0:
                queue.append(child)
    if len(ordered) != len(steps):
        raise ValueError("plan contains a depends_on cycle")
    return ordered


def resolve_path(path_value: str, base_dir: Path) -> Path:
    path = Path(path_value)
    return path if path.is_absolute() else base_dir / path


def _artifact_validators(validators: list[str]) -> list[str]:
    return [str(name) for name in validators if str(name) in ARTIFACT_VALIDATORS]


def _max_retries(plan: dict[str, Any], step: dict[str, Any], override: int | None) -> int:
    if override is not None:
        return max(0, int(override))
    retry = step.get("retry_policy") or {}
    if "max_retries" in retry:
        return max(0, int(retry["max_retries"]))
    execution = plan.get("execution") or {}
    return max(0, int(execution.get("max_retries", 0) or 0))


def _timeout(plan: dict[str, Any], step: dict[str, Any]) -> int | None:
    if step.get("timeout_seconds"):
        return int(step["timeout_seconds"])
    execution = plan.get("execution") or {}
    value = execution.get("timeout_seconds")
    return int(value) if value else None


def _run_command(command: Any, *, cwd: Path, timeout: int | None) -> subprocess.CompletedProcess:
    tokens = normalize_command(command)
    safe, message = is_safe_brick_command(tokens)
    if not safe:
        raise ValueError(message)
    return subprocess.run(tokens, cwd=str(cwd), text=True, capture_output=True, timeout=timeout)


def _validate_outputs(
    step: dict[str, Any],
    *,
    base_dir: Path,
    handle,
) -> bool:
    validators = _artifact_validators(step.get("validators", []) or [])
    if not validators:
        return True
    outputs = [str(path) for path in (step.get("outputs", []) or [])]
    if not outputs:
        return True
    validator_args = step.get("validator_args") or {}
    all_ok = True
    for output in outputs:
        output_path = resolve_path(output, base_dir)
        for name in validators:
            write_event(handle, "validator_started", step_id=step.get("id"), validator=name, input=str(output_path))
        results = run_artifact_validators(output_path, validators, validator_args=validator_args)
        for result in results:
            write_event(handle, "validator_finished", step_id=step.get("id"), input=str(output_path), **result.as_dict())
        if has_failed(results):
            all_ok = False
    return all_ok


def run_plan(
    plan_path: str | Path,
    *,
    workdir: str | Path | None = None,
    execute: bool = False,
    max_retries: int | None = None,
    manifest: str | Path | None = None,
) -> dict[str, Any]:
    plan_path = Path(plan_path)
    base_dir = plan_path.parent.resolve()
    plan = load_json(plan_path)
    run_workdir = Path(workdir or plan.get("workdir") or base_dir)
    if not run_workdir.is_absolute():
        run_workdir = base_dir / run_workdir
    run_workdir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(manifest) if manifest else run_workdir / "run_manifest.jsonl"
    if not manifest_path.is_absolute():
        manifest_path = base_dir / manifest_path
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    issues = validate_agent_plan(plan, base_dir=base_dir, workdir=run_workdir)
    if has_errors(issues):
        manifest_path.write_text("", encoding="utf-8")
        with manifest_path.open("a", encoding="utf-8") as handle:
            for issue in issues:
                write_event(handle, "plan_issue", **issue.as_dict())
            write_event(handle, "plan_failed", plan=str(plan_path))
        return {"ok": False, "manifest": str(manifest_path), "issues": [issue.as_dict() for issue in issues]}

    summary = {"ok": True, "manifest": str(manifest_path), "steps": []}
    ordered = step_order(plan["steps"])
    with manifest_path.open("w", encoding="utf-8") as handle:
        write_event(handle, "run_started", plan=str(plan_path), execute=execute, workdir=str(run_workdir))
        for step in ordered:
            step_id = str(step["id"])
            retries = _max_retries(plan, step, max_retries)
            attempts = retries + 1
            step_ok = False
            for attempt in range(1, attempts + 1):
                display = command_to_display(step.get("command"))
                write_event(handle, "step_started", step_id=step_id, attempt=attempt, command=display)
                if not execute:
                    write_event(handle, "step_dry_run", step_id=step_id, attempt=attempt, outputs=step.get("outputs", []))
                    step_ok = True
                    break

                start = time.time()
                try:
                    result = _run_command(step.get("command"), cwd=base_dir, timeout=_timeout(plan, step))
                    elapsed = round(time.time() - start, 3)
                    write_event(
                        handle,
                        "step_finished",
                        step_id=step_id,
                        attempt=attempt,
                        returncode=result.returncode,
                        elapsed_seconds=elapsed,
                        stdout_tail=result.stdout[-4000:],
                        stderr_tail=result.stderr[-4000:],
                    )
                    command_ok = result.returncode == 0
                except Exception as exc:
                    command_ok = False
                    write_event(handle, "step_failed", step_id=step_id, attempt=attempt, error=f"{type(exc).__name__}: {exc}")

                validator_ok = False
                if command_ok:
                    validator_ok = _validate_outputs(step, base_dir=base_dir, handle=handle)

                if command_ok and validator_ok:
                    step_ok = True
                    break
                if attempt < attempts:
                    write_event(handle, "step_retry", step_id=step_id, next_attempt=attempt + 1)

            if step_ok:
                write_event(handle, "step_success", step_id=step_id, outputs=step.get("outputs", []))
                summary["steps"].append({"step_id": step_id, "status": "success"})
                continue
            write_event(handle, "step_failed_final", step_id=step_id, allow_failure=bool(step.get("allow_failure")))
            summary["steps"].append({"step_id": step_id, "status": "failed"})
            if not step.get("allow_failure"):
                summary["ok"] = False
                write_event(handle, "run_failed", failed_step=step_id)
                return summary
        write_event(handle, "run_finished", ok=summary["ok"])
    return summary
