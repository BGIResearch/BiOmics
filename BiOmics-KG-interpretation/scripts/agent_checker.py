"""Validation helpers for brick agent plans and artifacts."""

from __future__ import annotations

import csv
import json
import os
import re
import shlex
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

try:
    from kgpolicy import resolve_relation_aliases
except Exception:  # pragma: no cover
    resolve_relation_aliases = None


PLAN_VALIDATORS = {
    "kg_policy_not_broad",
    "max_hops_respected",
    "relation_alias_resolved",
    "no_output_collision",
    "depends_on_satisfied",
    "cuda_required",
    "ranking_frozen_before_validation",
    "stable_id_primary_identity",
}

ARTIFACT_VALIDATORS = {
    "file_exists",
    "non_empty_file",
    "json_schema_valid",
    "csv_columns_present",
    "kg_edge_metadata_columns_present",
    "kg_edge_metadata_required_for_co_mentioned",
    "kg_edge_metadata_required_for_comention",
    "co_mentioned_weak_cap",
    "evidence_strength_present",
    "manual_review_notes_present",
    "no_weak_only_strong_claim",
    "max_hops_respected",
}

UNSAFE_TOKENS = {";", "&&", "||", "|", ">", "<", "`", "$(", "\n", "\r"}
EDGE_METADATA_SUFFIXES = (
    "relation",
    "info_source",
    "original_relation",
    "condition",
    "relation_confidence",
)
WEAK_STRENGTH_PREFIXES = ("weak", "review", "co_mentioned")


@dataclass
class CheckIssue:
    level: str
    message: str
    step_id: str | None = None
    validator: str | None = None

    def as_dict(self) -> dict[str, Any]:
        data = {"level": self.level, "message": self.message}
        if self.step_id:
            data["step_id"] = self.step_id
        if self.validator:
            data["validator"] = self.validator
        return data


@dataclass
class ValidatorResult:
    validator: str
    status: str
    message: str
    level: str = "error"

    def as_dict(self) -> dict[str, Any]:
        return {
            "validator": self.validator,
            "status": self.status,
            "level": self.level,
            "message": self.message,
        }


def load_json(path: str | os.PathLike[str]) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def normalize_command(command: Any) -> list[str]:
    """Return a command token list for list-style and simple string commands."""
    if isinstance(command, (list, tuple)):
        return [str(item) for item in command]
    if isinstance(command, str):
        return shlex.split(command, posix=os.name != "nt")
    return []


def command_to_display(command: Any) -> str:
    if isinstance(command, str):
        return command
    return " ".join(str(item) for item in command or [])


def is_safe_brick_command(command: Any) -> tuple[bool, str]:
    display = command_to_display(command)
    for token in UNSAFE_TOKENS:
        if token in display:
            return False, f"unsafe shell token is not allowed in agent P0: {token}"
    tokens = normalize_command(command)
    if len(tokens) < 2:
        return False, "command must include at least 'brick <subcommand>'"
    if tokens[0] not in {"brick", "brick-gpu"}:
        return False, "agent P0 only allows commands starting with 'brick' or 'brick-gpu'"
    return True, "safe brick command"


def _flag_value(tokens: list[str], *flags: str) -> str | None:
    for index, token in enumerate(tokens):
        if token in flags and index + 1 < len(tokens):
            return tokens[index + 1]
        for flag in flags:
            prefix = flag + "="
            if token.startswith(prefix):
                return token[len(prefix) :]
    return None


def _has_flag(tokens: list[str], *flags: str) -> bool:
    return any(token in flags or any(token.startswith(flag + "=") for flag in flags) for token in tokens)


def _is_brick_qr(tokens: list[str]) -> bool:
    return len(tokens) >= 2 and tokens[0] == "brick" and tokens[1] == "qr"


def _command_hops(tokens: list[str]) -> int:
    value = _flag_value(tokens, "--hops", "--multi-hop")
    if not value:
        return 1
    try:
        return int(value)
    except ValueError:
        return 0


def _inside_workdir(path_value: str, workdir: Path, base_dir: Path) -> bool:
    if not path_value or "{" in path_value:
        return True
    raw_path = Path(path_value)
    if ".." in raw_path.parts:
        return False
    resolved = raw_path if raw_path.is_absolute() else (base_dir / raw_path)
    try:
        resolved = resolved.resolve()
        workdir = workdir.resolve()
        return resolved == workdir or workdir in resolved.parents
    except OSError:
        return False


def _topological_issues(steps: list[dict[str, Any]]) -> list[CheckIssue]:
    issues: list[CheckIssue] = []
    ids = [str(step.get("id", "")) for step in steps]
    id_set = set(ids)
    if len(id_set) != len(ids):
        issues.append(CheckIssue("error", "duplicate step id found", validator="depends_on_satisfied"))

    graph: dict[str, set[str]] = {step_id: set() for step_id in ids}
    incoming: dict[str, int] = {step_id: 0 for step_id in ids}
    for step_index, step in enumerate(steps):
        step_id = str(step.get("id", ""))
        for dep in step.get("depends_on", []) or []:
            dep = str(dep)
            if dep not in id_set:
                issues.append(
                    CheckIssue(
                        "error",
                        f"depends_on references missing step: {dep}",
                        step_id=step_id,
                        validator="depends_on_satisfied",
                    )
                )
                continue
            graph[dep].add(step_id)
            incoming[step_id] += 1

    queue = deque([step_id for step_id, count in incoming.items() if count == 0])
    visited = 0
    while queue:
        current = queue.popleft()
        visited += 1
        for child in graph[current]:
            incoming[child] -= 1
            if incoming[child] == 0:
                queue.append(child)
    if visited != len(ids):
        issues.append(CheckIssue("error", "depends_on contains a cycle", validator="depends_on_satisfied"))
    return issues


def _validate_kg_command(step: dict[str, Any]) -> list[CheckIssue]:
    issues: list[CheckIssue] = []
    step_id = str(step.get("id", ""))
    tokens = normalize_command(step.get("command"))
    if not _is_brick_qr(tokens):
        return issues

    broad = not (
        _has_flag(tokens, "--tt", "--target-type")
        or _has_flag(tokens, "-r", "--relation")
        or _has_flag(tokens, "--info-source")
        or _has_flag(tokens, "--target")
        or _has_flag(tokens, "--allow-broad-query")
    )
    if broad:
        issues.append(
            CheckIssue(
                "error",
                "KG query is too broad: add target type, relation, target, info_source, or explicit --allow-broad-query",
                step_id=step_id,
                validator="kg_policy_not_broad",
            )
        )

    hops = _command_hops(tokens)
    if hops < 1 or hops > 3:
        issues.append(
            CheckIssue(
                "error",
                f"KG query hops must be within 1-3, got {hops}",
                step_id=step_id,
                validator="max_hops_respected",
            )
        )

    relation = _flag_value(tokens, "-r", "--relation")
    if relation and resolve_relation_aliases is not None:
        _resolved, unknown, aliases = resolve_relation_aliases(relation)
        if aliases:
            issues.append(
                CheckIssue(
                    "error",
                    "relation aliases must be resolved before execution: " + ",".join(sorted(aliases)),
                    step_id=step_id,
                    validator="relation_alias_resolved",
                )
            )
        if unknown:
            issues.append(
                CheckIssue(
                    "error",
                    "unknown relation names: " + ",".join(unknown),
                    step_id=step_id,
                    validator="relation_alias_resolved",
                )
            )
    return issues


def validate_agent_plan(
    plan: dict[str, Any],
    *,
    base_dir: str | os.PathLike[str] | None = None,
    workdir: str | os.PathLike[str] | None = None,
) -> list[CheckIssue]:
    """Validate the deterministic P0 agent plan contract."""
    base = Path(base_dir or ".").resolve()
    plan_workdir = Path(workdir or plan.get("workdir") or base)
    if not plan_workdir.is_absolute():
        plan_workdir = base / plan_workdir

    issues: list[CheckIssue] = []
    if str(plan.get("schema_version", "")) not in {"1.0", "1", ""}:
        issues.append(CheckIssue("error", "unsupported schema_version; expected 1.0"))

    execution = plan.get("execution") or {}
    if execution.get("allow_shell"):
        issues.append(CheckIssue("error", "allow_shell=true is not supported in agent P0"))

    steps = plan.get("steps")
    if not isinstance(steps, list) or not steps:
        issues.append(CheckIssue("error", "agent plan must contain a non-empty steps list"))
        return issues

    issues.extend(_topological_issues(steps))

    output_owner: dict[str, str] = {}
    normalized_commands = [normalize_command(step.get("command")) for step in steps]
    freeze_indices = [
        index
        for index, tokens in enumerate(normalized_commands)
        if tokens[:3] == ["brick", "drug", "freeze"]
        or tokens[:3] == ["brick-gpu", "drug", "grace-rank"]
    ]
    for step_index, step in enumerate(steps):
        step_id = str(step.get("id") or "")
        if not step_id:
            issues.append(CheckIssue("error", "step is missing id", validator="depends_on_satisfied"))
            step_id = "<missing>"
        if not step.get("type"):
            issues.append(CheckIssue("error", "step is missing type", step_id=step_id))
        if "command" not in step:
            issues.append(CheckIssue("error", "step is missing command", step_id=step_id))
            continue

        safe, message = is_safe_brick_command(step.get("command"))
        if not safe:
            issues.append(CheckIssue("error", message, step_id=step_id))

        validators = step.get("validators", [])
        if not isinstance(validators, list):
            issues.append(CheckIssue("error", "validators must be a list", step_id=step_id))
        else:
            unknown = sorted(set(str(v) for v in validators) - PLAN_VALIDATORS - ARTIFACT_VALIDATORS)
            for name in unknown:
                issues.append(CheckIssue("error", f"unknown validator: {name}", step_id=step_id, validator=name))

        tokens = normalize_command(step.get("command"))
        if "cuda_required" in validators:
            if tokens[:3] != ["brick-gpu", "drug", "grace-rank"] or _flag_value(tokens, "--device") != "cuda":
                issues.append(CheckIssue("error", "CUDA GRACE steps must use 'brick-gpu drug grace-rank --device cuda'", step_id=step_id, validator="cuda_required"))
        if "stable_id_primary_identity" in validators and tokens[:3] != ["brick", "drug", "validate"]:
            issues.append(CheckIssue("error", "stable-ID evaluation must use 'brick drug validate'", step_id=step_id, validator="stable_id_primary_identity"))
        if "ranking_frozen_before_validation" in validators:
            if not any(index < step_index for index in freeze_indices):
                issues.append(CheckIssue("error", "ranking validation occurs before a ranking freeze", step_id=step_id, validator="ranking_frozen_before_validation"))

        issues.extend(_validate_kg_command(step))

        for output in step.get("outputs", []) or []:
            output = str(output)
            key = str((base / output).resolve()) if not Path(output).is_absolute() else str(Path(output).resolve())
            if key in output_owner and output_owner[key] != step_id:
                issues.append(
                    CheckIssue(
                        "error",
                        f"output collision with step {output_owner[key]}: {output}",
                        step_id=step_id,
                        validator="no_output_collision",
                    )
                )
            output_owner[key] = step_id
            if not _inside_workdir(output, plan_workdir, base):
                issues.append(
                    CheckIssue(
                        "error",
                        f"output must stay inside workdir: {output}",
                        step_id=step_id,
                        validator="no_output_collision",
                    )
                )
    return issues


def _read_csv_columns(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        return next(reader, [])


def _empty(value: Any) -> bool:
    if value is None:
        return True
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "none", "null"}


def _family_columns(columns: Iterable[str], suffix: str) -> list[str]:
    direct = {
        "relation": {"relation"},
        "info_source": {"info_source", "source"},
        "original_relation": {"original_relation", "evidence_text"},
        "condition": {"condition"},
        "relation_confidence": {"relation_confidence", "confidence"},
    }.get(suffix, {suffix})
    result = []
    for column in columns:
        if column in direct or column.endswith("." + suffix):
            result.append(column)
    return result


def _relation_columns(columns: Iterable[str]) -> list[str]:
    return _family_columns(columns, "relation")


def _co_mentioned_mask(df: pd.DataFrame) -> pd.Series:
    mask = pd.Series(False, index=df.index)
    relation_cols = _relation_columns(df.columns)
    review_cols = [
        column
        for column in df.columns
        if column in {"edge_support_class", "evidence_strength"} or column.endswith(".relation")
    ]
    for column in relation_cols + review_cols:
        values = df[column].fillna("").astype(str).str.lower()
        mask = mask | values.str.contains("co_mentioned|co-mentioned|cooccurrence|co-occurrence", regex=True)
    return mask


def _same_edge_field(relation_col: str, suffix: str, columns: Iterable[str]) -> str | None:
    if "." in relation_col:
        candidate = relation_col.rsplit(".", 1)[0] + "." + suffix
        if candidate in columns:
            return candidate
    direct = {
        "info_source": ["info_source", "source"],
        "original_relation": ["original_relation", "evidence_text"],
        "condition": ["condition"],
        "relation_confidence": ["relation_confidence", "confidence"],
    }.get(suffix, [suffix])
    for candidate in direct:
        if candidate in columns:
            return candidate
    family = _family_columns(columns, suffix)
    return family[0] if family else None


def _validator_file_exists(path: Path, _args: dict[str, Any]) -> ValidatorResult:
    if path.exists():
        return ValidatorResult("file_exists", "passed", f"file exists: {path}")
    return ValidatorResult("file_exists", "failed", f"file does not exist: {path}")


def _validator_non_empty_file(path: Path, _args: dict[str, Any]) -> ValidatorResult:
    if not path.exists():
        return ValidatorResult("non_empty_file", "failed", f"file does not exist: {path}")
    if path.stat().st_size > 0:
        return ValidatorResult("non_empty_file", "passed", f"file is non-empty: {path}")
    return ValidatorResult("non_empty_file", "failed", f"file is empty: {path}")


def _validator_csv_columns_present(path: Path, args: dict[str, Any]) -> ValidatorResult:
    required = args.get("required_columns") or args.get("columns") or []
    if isinstance(required, str):
        required = [item.strip() for item in required.split(",") if item.strip()]
    if not required:
        return ValidatorResult("csv_columns_present", "passed", "no required columns supplied")
    columns = set(_read_csv_columns(path))
    missing = [column for column in required if column not in columns]
    if missing:
        return ValidatorResult("csv_columns_present", "failed", "missing CSV columns: " + ",".join(missing))
    return ValidatorResult("csv_columns_present", "passed", "required CSV columns are present")


def _validator_json_schema_valid(path: Path, _args: dict[str, Any]) -> ValidatorResult:
    try:
        json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return ValidatorResult("json_schema_valid", "failed", f"invalid JSON: {exc}")
    return ValidatorResult("json_schema_valid", "passed", "JSON parses successfully")


def _validator_kg_edge_metadata_columns_present(path: Path, _args: dict[str, Any]) -> ValidatorResult:
    columns = _read_csv_columns(path)
    missing = [suffix for suffix in EDGE_METADATA_SUFFIXES if not _family_columns(columns, suffix)]
    if missing:
        return ValidatorResult(
            "kg_edge_metadata_columns_present",
            "failed",
            "missing KG edge metadata column family: " + ",".join(missing),
        )
    return ValidatorResult("kg_edge_metadata_columns_present", "passed", "KG edge metadata columns are present")


def _validator_kg_edge_metadata_required_for_co_mentioned(path: Path, _args: dict[str, Any]) -> ValidatorResult:
    df = pd.read_csv(path)
    if df.empty:
        return ValidatorResult("kg_edge_metadata_required_for_co_mentioned", "passed", "empty CSV has no co_mentioned rows")
    co_mask = _co_mentioned_mask(df)
    if not co_mask.any():
        return ValidatorResult("kg_edge_metadata_required_for_co_mentioned", "passed", "no co_mentioned rows found")

    bad_rows = set()
    relation_cols = _relation_columns(df.columns)
    for relation_col in relation_cols:
        rel_mask = df[relation_col].fillna("").astype(str).str.lower().str.contains("co_mentioned|co-mentioned", regex=True)
        if not rel_mask.any():
            continue
        info_col = _same_edge_field(relation_col, "info_source", df.columns)
        text_col = _same_edge_field(relation_col, "original_relation", df.columns)
        if not info_col or not text_col:
            bad_rows.update(df.index[rel_mask].tolist())
            continue
        missing = rel_mask & (df[info_col].apply(_empty) | df[text_col].apply(_empty))
        bad_rows.update(df.index[missing].tolist())

    if bad_rows:
        sample = ",".join(str(int(idx) + 2) for idx in sorted(bad_rows)[:10])
        return ValidatorResult(
            "kg_edge_metadata_required_for_co_mentioned",
            "failed",
            "co_mentioned rows require non-empty original_relation/evidence_text and info_source; CSV lines: " + sample,
        )
    return ValidatorResult(
        "kg_edge_metadata_required_for_co_mentioned",
        "passed",
        "co_mentioned rows have reviewable edge metadata",
    )


def _validator_co_mentioned_weak_cap(path: Path, _args: dict[str, Any]) -> ValidatorResult:
    df = pd.read_csv(path)
    if df.empty:
        return ValidatorResult("co_mentioned_weak_cap", "passed", "empty CSV has no co_mentioned rows")
    co_mask = _co_mentioned_mask(df)
    if not co_mask.any():
        return ValidatorResult("co_mentioned_weak_cap", "passed", "no co_mentioned rows found")
    if "evidence_strength" not in df.columns:
        return ValidatorResult("co_mentioned_weak_cap", "failed", "co_mentioned rows need evidence_strength before weak-cap validation")
    strength = df.loc[co_mask, "evidence_strength"].fillna("").astype(str).str.lower()
    weak = strength.str.startswith(WEAK_STRENGTH_PREFIXES)
    if "manual_review_required" in df.columns:
        manual = df.loc[co_mask, "manual_review_required"].fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})
        weak = weak | manual
    if not weak.all():
        return ValidatorResult("co_mentioned_weak_cap", "failed", "some co_mentioned rows are not capped as weak/review-required")
    return ValidatorResult("co_mentioned_weak_cap", "passed", "co_mentioned rows are capped as weak/review-required")


def _validator_evidence_strength_present(path: Path, _args: dict[str, Any]) -> ValidatorResult:
    df = pd.read_csv(path)
    if "evidence_strength" not in df.columns:
        return ValidatorResult("evidence_strength_present", "failed", "missing evidence_strength column")
    if df["evidence_strength"].apply(_empty).any():
        return ValidatorResult("evidence_strength_present", "failed", "evidence_strength contains empty values")
    return ValidatorResult("evidence_strength_present", "passed", "evidence_strength is present")


def _validator_manual_review_notes_present(path: Path, _args: dict[str, Any]) -> ValidatorResult:
    df = pd.read_csv(path)
    if df.empty:
        return ValidatorResult("manual_review_notes_present", "passed", "empty CSV has no review rows")
    weak_mask = pd.Series(False, index=df.index)
    if "evidence_strength" in df.columns:
        weak_mask = weak_mask | df["evidence_strength"].fillna("").astype(str).str.lower().str.startswith(WEAK_STRENGTH_PREFIXES)
    weak_mask = weak_mask | _co_mentioned_mask(df)
    if not weak_mask.any():
        return ValidatorResult("manual_review_notes_present", "passed", "no weak/co_mentioned rows require notes")
    note_ok = pd.Series(False, index=df.index)
    if "edge_review_note" in df.columns:
        note_ok = note_ok | ~df["edge_review_note"].apply(_empty)
    if "manual_review_required" in df.columns:
        note_ok = note_ok | df["manual_review_required"].fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})
    if not note_ok[weak_mask].all():
        return ValidatorResult("manual_review_notes_present", "failed", "weak/co_mentioned rows require edge_review_note or manual_review_required=true")
    return ValidatorResult("manual_review_notes_present", "passed", "weak/co_mentioned rows include review notes")


def _validator_no_weak_only_strong_claim(path: Path, _args: dict[str, Any]) -> ValidatorResult:
    text = path.read_text(encoding="utf-8", errors="ignore").lower()
    strong_claim = re.search(r"\b(strong|confirmed|causal|proves|validated)\b", text)
    weak_only = "weak_co_mentioned_review_required" in text and "direct" not in text
    if strong_claim and weak_only:
        return ValidatorResult("no_weak_only_strong_claim", "failed", "report appears to make a strong claim from weak-only evidence")
    return ValidatorResult("no_weak_only_strong_claim", "passed", "no weak-only strong claim detected")


def _validator_max_hops_respected(path: Path, args: dict[str, Any]) -> ValidatorResult:
    max_hops = int(args.get("max_hops", 3))
    if path.suffix.lower() != ".csv":
        return ValidatorResult("max_hops_respected", "passed", "non-CSV artifact skipped")
    df = pd.read_csv(path, nrows=500)
    if "path_hops" in df.columns:
        observed = pd.to_numeric(df["path_hops"], errors="coerce").max()
        observed = int(observed) if pd.notna(observed) else 1
    else:
        path_indices = []
        for column in df.columns:
            match = re.match(r"path\.(\d+)\.", str(column))
            if match:
                path_indices.append(int(match.group(1)))
        observed = max([idx for idx in path_indices if idx % 2 == 1], default=1)
        observed = (observed + 1) // 2
    if observed > max_hops:
        return ValidatorResult("max_hops_respected", "failed", f"observed {observed} hops exceeds max_hops={max_hops}")
    return ValidatorResult("max_hops_respected", "passed", f"observed hops <= {max_hops}")


VALIDATOR_FUNCTIONS = {
    "file_exists": _validator_file_exists,
    "non_empty_file": _validator_non_empty_file,
    "json_schema_valid": _validator_json_schema_valid,
    "csv_columns_present": _validator_csv_columns_present,
    "kg_edge_metadata_columns_present": _validator_kg_edge_metadata_columns_present,
    "kg_edge_metadata_required_for_co_mentioned": _validator_kg_edge_metadata_required_for_co_mentioned,
    "kg_edge_metadata_required_for_comention": _validator_kg_edge_metadata_required_for_co_mentioned,
    "co_mentioned_weak_cap": _validator_co_mentioned_weak_cap,
    "evidence_strength_present": _validator_evidence_strength_present,
    "manual_review_notes_present": _validator_manual_review_notes_present,
    "no_weak_only_strong_claim": _validator_no_weak_only_strong_claim,
    "max_hops_respected": _validator_max_hops_respected,
}


def run_artifact_validators(
    input_path: str | os.PathLike[str],
    validators: Iterable[str],
    *,
    validator_args: dict[str, Any] | None = None,
) -> list[ValidatorResult]:
    path = Path(input_path)
    args_by_name = defaultdict(dict)
    if validator_args:
        for key, value in validator_args.items():
            args_by_name[str(key)] = value or {}

    results: list[ValidatorResult] = []
    for name in validators:
        name = str(name)
        if name in PLAN_VALIDATORS and name not in VALIDATOR_FUNCTIONS:
            continue
        func = VALIDATOR_FUNCTIONS.get(name)
        if not func:
            results.append(ValidatorResult(name, "failed", f"unknown validator: {name}"))
            continue
        try:
            result = func(path, args_by_name.get(name, {}))
        except Exception as exc:
            result = ValidatorResult(name, "failed", f"validator raised {type(exc).__name__}: {exc}")
        results.append(result)
    return results


def has_errors(issues: Iterable[CheckIssue]) -> bool:
    return any(issue.level == "error" for issue in issues)


def has_failed(results: Iterable[ValidatorResult]) -> bool:
    return any(result.status == "failed" and result.level == "error" for result in results)
