"""Task contracts for dataset-agnostic BiOmics-BRICK reasoning.

The contract layer maps benchmark rows or user questions into a small set of
task types. Dataset runners should provide IO/cache/evaluation, while BRICK
uses the contract to choose retrieval, evidence boundary, reasoning, and
rendering behavior.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal, Optional


TaskType = Literal[
    "ClosedChoice",
    "BinaryClaim",
    "ShortEntity",
    "EntityList",
    "EntityProfile",
    "MechanisticSummary",
    "PathAnswerSet",
]


@dataclass
class TaskContract:
    task_type: TaskType
    question: str
    expected_answer_type: Optional[str] = None
    source_mentions: list[str] = field(default_factory=list)
    target_mentions: list[str] = field(default_factory=list)
    options: list[str] = field(default_factory=list)
    output_schema: dict[str, Any] = field(default_factory=dict)
    evidence_policy: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "TaskContract":
        return cls(
            task_type=data["task_type"],
            question=data.get("question", ""),
            expected_answer_type=data.get("expected_answer_type"),
            source_mentions=list(data.get("source_mentions") or []),
            target_mentions=list(data.get("target_mentions") or []),
            options=list(data.get("options") or []),
            output_schema=dict(data.get("output_schema") or {}),
            evidence_policy=dict(data.get("evidence_policy") or {}),
            metadata=dict(data.get("metadata") or {}),
        )


def _split_options(text: str) -> list[str]:
    return [part.strip() for part in re.split(r"\s*,\s*", text or "") if part.strip()]


def _split_targets(text: str) -> list[str]:
    text = (text or "").strip()
    if not text:
        return []
    if " and " in text:
        left, right = text.rsplit(" and ", 1)
        return [left.strip(), right.strip()]
    return [text]


def _entity_type_from_label(label: str, default: str = "Gene") -> str:
    label = (label or "").casefold()
    if "drug" in label or "chemical" in label:
        return "Chemical"
    if "disease" in label:
        return "Disease"
    if "mutation" in label or "variant" in label or "snp" in label:
        return "Mutation"
    if "protein" in label:
        return "Protein"
    return default


def infer_biomix_mcq_contract(row: dict[str, Any]) -> TaskContract:
    question = str(row.get("question") or row.get("text") or "")
    match = re.search(
        r"which\s+([A-Za-z_/ -]+?)\s+is\s+associated\s+with\s+(.+?)\.\s*Given list is:\s*(.+)$",
        question,
        flags=re.I | re.S,
    )
    if match:
        expected_type = _entity_type_from_label(match.group(1), default="Gene")
        targets = _split_targets(match.group(2))
        options = _split_options(match.group(3))
    else:
        expected_type = "Gene"
        targets = []
        options = _split_options(question)
    return TaskContract(
        task_type="ClosedChoice",
        question=question,
        expected_answer_type=expected_type,
        source_mentions=options,
        target_mentions=targets,
        options=options,
        output_schema={"format": "option_text"},
        evidence_policy={
            "option_anchor": "hard",
            "weak_evidence_policy": "snippet_supported_only",
            "scope_aggregation": "option",
            "candidate_stage_policy": "target_retained_only",
            "evidence_ranking": "on",
            "score_mode": "dual_disease_support",
            "hub_penalty": True,
            "fallback": "hermes_if_low_margin",
        },
        metadata={"dataset": "biomix", "task": "mcq", "id": row.get("id")},
    )


def infer_biomix_true_false_contract(row: dict[str, Any]) -> TaskContract:
    question = str(row.get("question") or row.get("text") or "")
    return TaskContract(
        task_type="BinaryClaim",
        question=question,
        output_schema={"format": "true_false"},
        evidence_policy={
            "solver": "kg_proposition_verifier",
            "fallback": "hermes_biomics",
            "weak_evidence_policy": "context_only",
            "use_co_mentioned_as_evidence": False,
            "kg_no_hit_is_negative": False,
            "reasoning_policy": "binary_claim_kg_verifier_then_fallback",
        },
        metadata={"dataset": "biomix", "task": "true_false", "id": row.get("id")},
    )


def infer_biomix_openqa_contract(row: dict[str, Any]) -> TaskContract:
    question = str(row.get("question") or row.get("text") or "")
    match = re.search(r"\bwhat\s+is\s+(.+?)\??$", question.strip(), flags=re.I)
    entity = match.group(1).strip(" ?.") if match else ""
    source_type = "Gene" if re.fullmatch(r"[A-Z0-9-]{2,}", entity) else None
    return TaskContract(
        task_type="EntityProfile",
        question=question,
        expected_answer_type=source_type,
        source_mentions=[entity] if entity else [],
        output_schema={"format": "profile_text"},
        evidence_policy={
            "weak_evidence_policy": "context_only",
            "profile_sections": ["aliases", "definition", "functions", "diseases", "pathways", "cells"],
        },
        metadata={"dataset": "biomix", "task": "open_qa", "id": row.get("id")},
    )


def infer_bioasq_contract(row: dict[str, Any]) -> TaskContract:
    question = str(row.get("body") or row.get("question") or "")
    qtype = str(row.get("type") or row.get("question_type") or "").casefold()
    mapping: dict[str, TaskType] = {
        "yesno": "BinaryClaim",
        "factoid": "ShortEntity",
        "list": "EntityList",
        "summary": "MechanisticSummary",
    }
    task_type = mapping.get(qtype, "ShortEntity")
    schema = {"format": qtype or "text"}
    if task_type == "BinaryClaim":
        schema = {"format": "yes_no"}
    elif task_type == "ShortEntity":
        schema = {"format": "bioasq_factoid"}
    elif task_type == "EntityList":
        schema = {"format": "bioasq_list"}
    evidence_policy: dict[str, Any] = {"weak_evidence_policy": "context_only"}
    if task_type in {"ShortEntity", "EntityList"}:
        evidence_policy = {
            "kg_relation_policy": "curated_preferred",
            "weak_evidence_policy": "snippet_supported",
            "allow_weak_final_types": "factoid,list",
            "co_mentioned_max_share": 0.15,
            "co_mentioned_max_rows_per_question": 3,
            "candidate_stage_policy": "target_retained_only",
        }
    elif task_type == "BinaryClaim":
        evidence_policy = {
            "weak_evidence_policy": "context_only",
            "kg_no_hit_is_negative": False,
        }
    return TaskContract(
        task_type=task_type,
        question=question,
        output_schema=schema,
        evidence_policy=evidence_policy,
        metadata={"dataset": "bioasq", "type": qtype, "id": row.get("id")},
    )


def _biohopr_node_type(label: str | None) -> str | None:
    key = (label or "").strip().casefold()
    if not key:
        return None
    if "gene" in key or "protein" in key:
        return "Gene"
    if "drug" in key or "chemical" in key:
        return "Chemical"
    if "effect" in key or "phenotype" in key:
        return "Phenotype"
    if "disease" in key:
        return "Disease"
    if "pathway" in key:
        return "Pathway"
    if "function" in key:
        return "Function"
    if "process" in key:
        return "Process"
    if "cell" in key:
        return "Cell"
    if "tissue" in key:
        return "Tissue"
    return None


def _strip_biohopr_prompt_tail(text: str) -> str:
    return re.split(r"\n\s*Just give me|\n\s*Answer\s*:", text or "", maxsplit=1, flags=re.I)[0].strip()


def _parse_biohopr_question(question: str) -> dict[str, Any]:
    text = _strip_biohopr_prompt_tail(question).strip()

    def build(
        *,
        hop_count: int,
        target_type_label: str,
        source_type_label: str,
        source: str,
        intermediate_type_label: str | None = None,
    ) -> dict[str, Any]:
        out: dict[str, Any] = {
            "hop_count": hop_count,
            "target_type_label": target_type_label.strip(),
            "source_type_label": source_type_label.strip(),
            "source": source.strip(),
        }
        if intermediate_type_label:
            out["intermediate_type_label"] = intermediate_type_label.strip()
        return out

    # BioHopR stores only natural-language questions in this benchmark copy.
    # Keep these template parses deterministic so graph traversal receives a
    # typed source entity instead of issuing an empty/ill-typed KG query.
    two_hop_templates: list[tuple[str, Any]] = [
        (
            r"^Name\s+a[n]?\s+(.+?)\s+that\s+is\s+related\s+to\s+a[n]?\s+(.+?)\s+that\s+is\s+associated\s+with\s+(.+?)\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label=m.group(1),
                intermediate_type_label=m.group(2),
                source_type_label=m.group(3),
                source=m.group(4),
            ),
        ),
        (
            r"^Name\s+a\s+drug\s+that\s+can\s+treat\s+a\s+disease\s+that\s+is\s+associated\s+with\s+(.+?)\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label="drug",
                intermediate_type_label="disease",
                source_type_label=m.group(1),
                source=m.group(2),
            ),
        ),
        (
            r"^Name\s+a\s+drug\s+that\s+can\s+treat\s+a\s+disease\s+that\s+has\s+phenotype\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label="drug",
                intermediate_type_label="disease",
                source_type_label="effect/phenotype",
                source=m.group(1),
            ),
        ),
        (
            r"^Name\s+a[n]?\s+gene/protein\s+that\s+is\s+related\s+to\s+disease\s+that\s+is\s+treated\s+by\s+drug\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label="gene/protein",
                intermediate_type_label="disease",
                source_type_label="drug",
                source=m.group(1),
            ),
        ),
        (
            r"^Name\s+a[n]?\s+effect/phenotype\s+which\s+is\s+a\s+side\s+effect\s+of\s+drug\s+that\s+can\s+treat\s+disease\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label="effect/phenotype",
                intermediate_type_label="drug",
                source_type_label="disease",
                source=m.group(1),
            ),
        ),
        (
            r"^Name\s+a[n]?\s+disease\s+that\s+is\s+treated\s+by\s+a\s+drug\s+that\s+has\s+side\s+effect\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label="disease",
                intermediate_type_label="drug",
                source_type_label="effect/phenotype",
                source=m.group(1),
            ),
        ),
        (
            r"^Name\s+a[n]?\s+disease\s+that\s+is\s+treated\s+by\s+a\s+drug\s+that\s+is\s+associated\s+with\s+gene/protein\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label="disease",
                intermediate_type_label="drug",
                source_type_label="gene/protein",
                source=m.group(1),
            ),
        ),
        (
            r"^Name\s+a\s+drug\s+that\s+is\s+associated\s+with\s+gene/protein\s+that\s+is\s+related\s+to\s+disease\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label="drug",
                intermediate_type_label="gene/protein",
                source_type_label="disease",
                source=m.group(1),
            ),
        ),
        (
            r"^Name\s+a\s+drug\s+that\s+has\s+a\s+side\s+effect\s+which\s+is\s+also\s+a\s+phenotype\s+of\s+a\s+disease\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label="drug",
                intermediate_type_label="effect/phenotype",
                source_type_label="disease",
                source=m.group(1),
            ),
        ),
        (
            r"^Name\s+a[n]?\s+effect/phenotype\s+of\s+a\s+disease\s+that\s+is\s+treated\s+by\s+drug\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label="effect/phenotype",
                intermediate_type_label="disease",
                source_type_label="drug",
                source=m.group(1),
            ),
        ),
        (
            r"^Name\s+a[n]?\s+gene/protein\s+that\s+is\s+associated\s+with\s+drug\s+that\s+can\s+treat\s+disease\s+(.+?)\.$",
            lambda m: build(
                hop_count=2,
                target_type_label="gene/protein",
                intermediate_type_label="drug",
                source_type_label="disease",
                source=m.group(1),
            ),
        ),
    ]
    for pattern, factory in two_hop_templates:
        match = re.search(pattern, text, flags=re.I | re.S)
        if match:
            return factory(match)

    one_hop_templates: list[tuple[str, Any]] = [
        (
            r"^Name\s+a\s+drug\s+that\s+can\s+treat\s+disease\s+(.+?)\.$",
            lambda m: build(hop_count=1, target_type_label="drug", source_type_label="disease", source=m.group(1)),
        ),
        (
            r"^Name\s+a[n]?\s+effect/phenotype\s+which\s+is\s+a\s+side\s+effect\s+of\s+drug\s+(.+?)\.$",
            lambda m: build(hop_count=1, target_type_label="effect/phenotype", source_type_label="drug", source=m.group(1)),
        ),
        (
            r"^Name\s+a[n]?\s+disease\s+that\s+is\s+treated\s+by\s+drug\s+(.+?)\.$",
            lambda m: build(hop_count=1, target_type_label="disease", source_type_label="drug", source=m.group(1)),
        ),
        (
            r"^Name\s+a\s+drug\s+that\s+has\s+side\s+effect\s+(.+?)\.$",
            lambda m: build(hop_count=1, target_type_label="drug", source_type_label="effect/phenotype", source=m.group(1)),
        ),
        (
            r"^Name\s+a\s+drug\s+that\s+is\s+associated\s+with\s+gene/protein\s+(.+?)\.$",
            lambda m: build(hop_count=1, target_type_label="drug", source_type_label="gene/protein", source=m.group(1)),
        ),
        (
            r"^Name\s+a[n]?\s+effect/phenotype\s+of\s+disease\s+(.+?)\.$",
            lambda m: build(hop_count=1, target_type_label="effect/phenotype", source_type_label="disease", source=m.group(1)),
        ),
        (
            r"^Name\s+a[n]?\s+gene/protein\s+that\s+is\s+associated\s+with\s+drug\s+(.+?)\.$",
            lambda m: build(hop_count=1, target_type_label="gene/protein", source_type_label="drug", source=m.group(1)),
        ),
        (
            r"^Name\s+a[n]?\s+(.+?)\s+that\s+is\s+related\s+to\s+(.+?)\s+(.+?)\.$",
            lambda m: build(hop_count=1, target_type_label=m.group(1), source_type_label=m.group(2), source=m.group(3)),
        ),
    ]
    for pattern, factory in one_hop_templates:
        match = re.search(pattern, text, flags=re.I | re.S)
        if match:
            return factory(match)
    return {"hop_count": row_hops(question), "source": ""}


def row_hops(question: str) -> int:
    text = (question or "").casefold()
    return 2 if "that is associated with" in text and "that is related to a" in text else 1


def infer_biohopr_contract(row: dict[str, Any]) -> TaskContract:
    question = str(row.get("question") or row.get("text") or "")
    parsed = _parse_biohopr_question(question)
    hop_count = int(row.get("hops") or row.get("hop") or parsed.get("hop_count") or 1)
    source_type = _biohopr_node_type(row.get("source_type") or row.get("query_type") or parsed.get("source_type_label"))
    target_type = _biohopr_node_type(row.get("target_type") or parsed.get("target_type_label"))
    intermediate_type = _biohopr_node_type(row.get("intermediate_type") or row.get("bridge_type") or parsed.get("intermediate_type_label"))
    source = row.get("source") or row.get("query_entity") or parsed.get("source") or ""
    record_id = row.get("id") if row.get("id") is not None else row.get("sample_id")
    return TaskContract(
        task_type="PathAnswerSet",
        question=question,
        expected_answer_type=target_type,
        source_mentions=[str(source)] if source else [],
        target_mentions=[str(row.get("target"))] if row.get("target") else [],
        output_schema={"format": "answers_with_paths"},
        evidence_policy={
            "max_hops": hop_count,
            "source_type": source_type,
            "target_type": target_type,
            "intermediate_type": intermediate_type,
            "deterministic_kg_path": True,
            "kg_relation_policy": "curated_preferred",
            "weak_evidence_policy": "context_only",
            "co_mentioned_main_path_edge": False,
            "co_mentioned_fallback": True,
            "fallback_requires_type_match": True,
            "path_template_rerank": True,
            "hub_penalty": True,
            "generic_disease_penalty": True,
        },
        metadata={
            "dataset": "biohopr",
            "id": record_id,
            "sample_id": row.get("sample_id"),
            "structured_query_source": row.get("structured_query_source") or row.get("query_source") or "question",
            "query_extraction_method": "biohopr_question_template_parser",
            "query_entity": row.get("query_entity"),
            "query_type": row.get("query_type"),
            "bridge_entity": row.get("bridge_entity"),
            "bridge_type": row.get("bridge_type"),
            **parsed,
        },
    )


def infer_task_contract(row: dict[str, Any], dataset: str = "", task: str = "") -> TaskContract:
    dataset_key = (dataset or str(row.get("dataset") or "")).casefold()
    task_key = (task or str(row.get("task") or row.get("type") or "")).casefold()
    if dataset_key == "biomix" and ("mcq" in task_key or "given list is:" in str(row.get("question", "")).casefold()):
        return infer_biomix_mcq_contract(row)
    if dataset_key == "biomix" and (task_key in {"true_false", "tf", "binary_claim"} or "true_false" in task_key):
        return infer_biomix_true_false_contract(row)
    if dataset_key == "biomix" and (task_key in {"open_qa", "openqa", "pure_question", "short"} or "what is" in str(row.get("question", "")).casefold()):
        return infer_biomix_openqa_contract(row)
    if dataset_key == "bioasq" or task_key in {"yesno", "factoid", "list", "summary"}:
        return infer_bioasq_contract(row)
    if dataset_key == "biohopr" or "hop" in task_key:
        return infer_biohopr_contract(row)
    question = str(row.get("question") or row.get("text") or "")
    return TaskContract(
        task_type="ShortEntity",
        question=question,
        output_schema={"format": "text"},
        evidence_policy={"weak_evidence_policy": "context_only"},
        metadata={"dataset": dataset, "task": task, "id": row.get("id")},
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Infer a BiOmics-BRICK task contract")
    parser.add_argument("-i", "--input", help="JSON row or JSONL file")
    parser.add_argument("--dataset", default="")
    parser.add_argument("--task", default="")
    parser.add_argument("-o", "--output")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.input:
        path = Path(args.input)
        text = path.read_text(encoding="utf-8")
        rows = [json.loads(line) for line in text.splitlines() if line.strip()] if path.suffix == ".jsonl" else json.loads(text)
    else:
        rows = json.load(__import__("sys").stdin)
    if isinstance(rows, dict):
        result: Any = infer_task_contract(rows, dataset=args.dataset, task=args.task).to_dict()
    else:
        result = [infer_task_contract(row, dataset=args.dataset, task=args.task).to_dict() for row in rows]
    output = json.dumps(result, ensure_ascii=False, indent=2)
    if args.output:
        Path(args.output).write_text(output + "\n", encoding="utf-8")
    else:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
