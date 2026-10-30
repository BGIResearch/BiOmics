"""Task-specific answer renderers for BiOmics-BRICK contracts."""

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


def _answer_text(answer: Any) -> str:
    if isinstance(answer, dict):
        for key in ["answer", "text", "option", "value"]:
            if key in answer:
                return str(answer.get(key) or "").strip()
    if isinstance(answer, list):
        return str(answer[0]).strip() if answer else ""
    return str(answer or "").strip()


def render_closed_choice(answer: Any, contract: Any) -> str:
    contract = _as_contract(contract)
    text = _answer_text(answer)
    options = list(getattr(contract, "options", []) or [])
    for option in options:
        if text.casefold() == str(option).casefold():
            return str(option)
    return text


def render_binary_claim(answer: Any, contract: Any = None, dataset_format: str = "") -> str:
    text = _answer_text(answer).casefold()
    if text in {"yes", "true", "1", "y"}:
        return "yes" if dataset_format == "bioasq" else "true"
    if text in {"no", "false", "0", "n"}:
        return "no" if dataset_format == "bioasq" else "false"
    return _answer_text(answer)


def render_bioasq_factoid(answer: Any) -> list[list[str]]:
    if isinstance(answer, list):
        values = answer
    elif isinstance(answer, dict) and isinstance(answer.get("answers"), list):
        values = answer["answers"]
    else:
        values = [_answer_text(answer)]
    return [[str(value).strip()] for value in values if str(value).strip()]


def render_bioasq_list(answer: Any) -> list[list[str]]:
    if isinstance(answer, dict) and isinstance(answer.get("answers"), list):
        values = answer["answers"]
    elif isinstance(answer, list):
        values = answer
    else:
        values = [part.strip() for part in _answer_text(answer).split(",")]
    return [[str(value).strip()] for value in values if str(value).strip()]


def render_answer(answer: Any, contract: Any, dataset_format: str = "") -> Any:
    contract = _as_contract(contract)
    task_type = getattr(contract, "task_type", "")
    output_schema = dict(getattr(contract, "output_schema", {}) or {})
    fmt = dataset_format or str(output_schema.get("format") or "")
    if task_type == "ClosedChoice":
        return render_closed_choice(answer, contract)
    if task_type == "BinaryClaim":
        return render_binary_claim(answer, contract, dataset_format=fmt)
    if fmt == "bioasq_factoid":
        return render_bioasq_factoid(answer)
    if fmt == "bioasq_list":
        return render_bioasq_list(answer)
    if task_type == "PathAnswerSet" and isinstance(answer, dict):
        return {"answers": answer.get("answers", []), "paths": answer.get("paths", [])}
    return _answer_text(answer)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Render an answer through a task contract")
    parser.add_argument("--contract", required=True, help="TaskContract JSON")
    parser.add_argument("--answer", required=True, help="Raw answer JSON or text")
    parser.add_argument("--dataset-format", default="")
    parser.add_argument("-o", "--output")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    contract = json.loads(Path(args.contract).read_text(encoding="utf-8"))
    try:
        answer = json.loads(args.answer)
    except json.JSONDecodeError:
        answer = args.answer
    rendered = render_answer(answer, contract, dataset_format=args.dataset_format)
    output = json.dumps(rendered, ensure_ascii=False, indent=2)
    if args.output:
        Path(args.output).write_text(output + "\n", encoding="utf-8")
    else:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
