"""Evaluation adapters keyed by BiOmics-BRICK task contracts.

This module intentionally stays small: benchmark-specific runners can compute
official metrics, while BRICK exposes common normalization and exact-match
helpers for contract-level diagnostics.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Optional

try:
    from task_contract import TaskContract
except Exception:  # pragma: no cover
    TaskContract = None


def _as_contract(contract: Any):
    if TaskContract is not None and isinstance(contract, TaskContract):
        return contract
    if TaskContract is not None and isinstance(contract, dict):
        return TaskContract.from_dict(contract)
    return contract


def normalize_scalar(value: Any) -> str:
    if isinstance(value, list):
        if len(value) == 1:
            return normalize_scalar(value[0])
        return "|".join(normalize_scalar(item) for item in value)
    if isinstance(value, dict):
        for key in ["answer", "text", "value", "option"]:
            if key in value:
                return normalize_scalar(value.get(key))
    return " ".join(str(value or "").strip().casefold().split())


def normalize_answer_set(value: Any) -> set[str]:
    if value is None:
        return set()
    if isinstance(value, dict):
        value = value.get("answers", value.get("answer", []))
    if not isinstance(value, list):
        value = [value]
    out = set()
    for item in value:
        if isinstance(item, list):
            out.update(normalize_answer_set(item))
        else:
            text = normalize_scalar(item)
            if text:
                out.add(text)
    return out


def exact_match(prediction: Any, gold: Any, contract: Any) -> dict[str, Any]:
    contract = _as_contract(contract)
    task_type = getattr(contract, "task_type", "")
    if task_type in {"EntityList", "PathAnswerSet"}:
        pred_set = normalize_answer_set(prediction)
        gold_set = normalize_answer_set(gold)
        correct = pred_set == gold_set
        return {
            "metric": "set_exact_match",
            "correct": correct,
            "prediction": sorted(pred_set),
            "gold": sorted(gold_set),
        }
    pred = normalize_scalar(prediction)
    target = normalize_scalar(gold)
    return {
        "metric": "exact_match",
        "correct": pred == target,
        "prediction": pred,
        "gold": target,
    }


def selected_strong_overlap(
    rows: list[dict[str, Any]],
    native_answers: dict[int, Any],
    gold_answers: dict[int, Any],
) -> dict[str, int]:
    counts = {
        "native_correct_evidence_correct": 0,
        "native_wrong_evidence_correct": 0,
        "native_correct_evidence_wrong": 0,
        "native_wrong_evidence_wrong": 0,
        "selected_strong_questions": 0,
    }
    for row in rows:
        kg = row.get("kg", {})
        if not kg.get("selected_option_has_strong_evidence"):
            continue
        qid = int(row.get("id"))
        native_ok = normalize_scalar(native_answers.get(qid)) == normalize_scalar(gold_answers.get(qid))
        evidence_ok = normalize_scalar(row.get("answer")) == normalize_scalar(gold_answers.get(qid))
        counts["selected_strong_questions"] += 1
        if native_ok and evidence_ok:
            counts["native_correct_evidence_correct"] += 1
        elif not native_ok and evidence_ok:
            counts["native_wrong_evidence_correct"] += 1
        elif native_ok and not evidence_ok:
            counts["native_correct_evidence_wrong"] += 1
        else:
            counts["native_wrong_evidence_wrong"] += 1
    return counts


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Evaluate a prediction with a task contract")
    parser.add_argument("--contract", required=True)
    parser.add_argument("--prediction", required=True)
    parser.add_argument("--gold", required=True)
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    contract = json.loads(Path(args.contract).read_text(encoding="utf-8"))
    try:
        prediction = json.loads(args.prediction)
    except json.JSONDecodeError:
        prediction = args.prediction
    try:
        gold = json.loads(args.gold)
    except json.JSONDecodeError:
        gold = args.gold
    print(json.dumps(exact_match(prediction, gold, contract), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
