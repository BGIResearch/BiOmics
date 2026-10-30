"""Evidence packet construction for BiOmics-BRICK KG rows.

This module is the boundary between KG retrieval and downstream answer
generation. QueryGraph should retrieve rows; reasoning/runners should consume
evidence packets where weak trace rows cannot silently become final answers.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any, Optional

try:  # pandas is a project dependency, but keep imports test-friendly.
    import pandas as pd
except Exception:  # pragma: no cover
    pd = None

try:
    from kgpolicy import canonicalize_info_source, load_kg_schema
except Exception:  # pragma: no cover
    def canonicalize_info_source(value: Optional[Any]) -> str:
        return ""

    def load_kg_schema(schema_file: Optional[str] = None) -> dict[str, Any]:
        path = Path(schema_file) if schema_file else Path(__file__).resolve().parents[1] / "references" / "kg_schema.json"
        return json.loads(path.read_text(encoding="utf-8"))

try:
    from task_contract import TaskContract
except Exception:  # pragma: no cover
    TaskContract = None


DEFAULT_STRONG_CURATED = {
    "causally_related_to",
    "disease_relation",
    "pharmaco_effect_on",
    "regulated",
    "produce",
    "participate_in",
    "interact_with",
    "functionally_related_to",
    "is_a",
}
DEFAULT_WEAK_TRACE = {"co_mentioned", "mentioned", "no_relation"}
DEFAULT_KG_RELATION_POLICY = "curated_preferred"
DEFAULT_WEAK_EVIDENCE_POLICY = "snippet_supported"
TARGET_RETAINED_CANDIDATE_STAGES = {"initial", "fallback_no_relation", "fallback_relax_target_type"}
DROP_RELATIONS = {"", "no_relation"}
KG_RELATION_POLICIES = {
    "curated_preferred",
    "curated_only",
    "no_co_mentioned",
    "p2_current",
    "p2_no_co_mentioned",
    "p2_curated_only",
    "p2_curated_preferred_snippet",
}
WEAK_EVIDENCE_POLICIES = {
    "snippet_supported",
    "context_only",
    "allow_weak_candidates",
    "include",
    "drop",
}
WEAK_FINAL_TASK_ALIASES = {
    "factoid": "shortentity",
    "short_entity": "shortentity",
    "shortentity": "shortentity",
    "list": "entitylist",
    "entity_list": "entitylist",
    "entitylist": "entitylist",
}
DEFAULT_ALLOW_WEAK_FINAL_TYPES = {"shortentity", "entitylist"}
SNIPPET_RELATION_CUES = (
    "associate",
    "associated",
    "association",
    "correlat",
    "relat",
    "caus",
    "induc",
    "inhibit",
    "activat",
    "regulat",
    "express",
    "marker",
    "biomarker",
    "target",
    "treat",
    "therap",
    "bind",
    "interact",
    "mutat",
    "variant",
    "risk",
    "response",
    "predict",
    "mediate",
    "involve",
    "enrich",
    "upregulat",
    "downregulat",
    "pathway",
    "phenotype",
)
SNIPPET_NEGATION_PATTERNS = (
    r"\bno\s+(?:significant\s+)?association\b",
    r"\bnot\s+(?:significantly\s+)?associated\b",
    r"\bnot\s+related\b",
    r"\bnot\s+correlated\b",
    r"\black(?:ed|s)?\s+of\s+association\b",
    r"\bfailed\s+to\s+(?:show|detect|demonstrate)\b",
    r"\bwithout\s+(?:evidence|association)\b",
)
HIGH_PRIORITY_RELATIONS = {
    "disease_relation",
    "causally_related_to",
    "pharmaco_effect_on",
    "regulated",
    "produce",
    "interact_with",
}
MEDIUM_PRIORITY_RELATIONS = {"participate_in", "functionally_related_to"}
LOW_PRIORITY_RELATIONS = {"is_a", "has_part", "equals_to", "chemical_relation", "chemicals_relation"}
HIGH_RELIABILITY_SOURCES = {
    "CLINVAR",
    "GWAS_CATALOG",
    "GWAS_CATAOG",
    "OMIM",
    "PHARMGKB",
    "DRUGBANK",
    "DRUGCENTRAL",
    "ORPHANET",
    "MEDGEN",
    "MONDO",
    "HPO",
}
MEDIUM_RELIABILITY_SOURCES = {
    "AGRKB",
    "NCIT",
    "GO",
    "GOA",
    "REACTOME",
    "KEGG",
    "EFO",
    "GWAS",
}
STAGE_SPECIFICITY = {
    "initial": 1.0,
    "fallback_no_relation": 0.85,
    "fallback_relax_target_type": 0.75,
    "fallback_drop_target": 0.25,
    "fallback_drop_target_no_relation": 0.10,
}


def _clean(value: Any) -> str:
    if value is None:
        return ""
    if pd is not None:
        try:
            if pd.isna(value):
                return ""
        except (TypeError, ValueError):
            pass
    text = str(value).strip()
    return "" if text.casefold() in {"nan", "none", "null", "[]", "{}"} else text


def _norm(value: Any) -> str:
    return _clean(value).casefold().replace("-", "_").replace(" ", "_")


def _normalize_kg_relation_policy(policy: Optional[str]) -> str:
    value = _norm(policy or DEFAULT_KG_RELATION_POLICY)
    aliases = {
        "p2_curated_preferred_snippet": "curated_preferred",
        "p2_curated_only": "curated_only",
        "p2_no_co_mentioned": "no_co_mentioned",
        "p2_current": "curated_preferred",
    }
    value = aliases.get(value, value)
    if value not in {"curated_preferred", "curated_only", "no_co_mentioned"}:
        raise ValueError(
            "kg_relation_policy must be curated_preferred, curated_only, "
            "no_co_mentioned, p2_current, p2_no_co_mentioned, "
            "p2_curated_only, or p2_curated_preferred_snippet"
        )
    return value


def _normalize_weak_policy(policy: Optional[str]) -> str:
    value = _norm(policy or DEFAULT_WEAK_EVIDENCE_POLICY)
    if value == "include":
        value = "allow_weak_candidates"
    if value not in WEAK_EVIDENCE_POLICIES:
        raise ValueError(
            "weak_evidence_policy must be snippet_supported, context_only, "
            "allow_weak_candidates, include, or drop"
        )
    return value


def _row_get(row: Mapping[str, Any], candidates: Iterable[str]) -> str:
    for column in candidates:
        if column in row:
            value = _clean(row.get(column))
            if value:
                return value
    return ""


def _row_dict(row: Any) -> dict[str, Any]:
    if hasattr(row, "to_dict"):
        raw = row.to_dict()
    elif isinstance(row, Mapping):
        raw = dict(row)
    else:
        raw = dict(row)
    return {str(key): value for key, value in raw.items()}


def _rows_iter(rows: Any) -> list[dict[str, Any]]:
    if rows is None:
        return []
    if pd is not None and isinstance(rows, pd.DataFrame):
        return [_row_dict(row) for _, row in rows.iterrows()]
    if isinstance(rows, Mapping):
        return [_row_dict(rows)]
    return [_row_dict(row) for row in rows]


def _split_terms(values: Any) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str):
        raw = [values]
    elif isinstance(values, Iterable):
        raw = [str(value) for value in values if _clean(value)]
    else:
        raw = [str(values)]
    terms: list[str] = []
    for value in raw:
        for piece in re.split(r"[;,|]+", value):
            term = re.sub(r"\s+", " ", piece).strip()
            if len(term) >= 2:
                terms.append(term)
    seen = set()
    out = []
    for term in terms:
        key = term.casefold()
        if key not in seen:
            seen.add(key)
            out.append(term)
    return out


def _contains_any(text: str, terms: list[str]) -> bool:
    if not terms:
        return True
    haystack = text.casefold()
    return any(term.casefold() in haystack for term in terms if term)


def _text_contains_term(text: str, term: str) -> bool:
    term = re.sub(r"\s+", " ", _clean(term)).strip()
    if len(term) < 2:
        return False
    haystack = text.casefold()
    needle = term.casefold()
    if needle in haystack:
        return True
    compact_haystack = re.sub(r"[^a-z0-9]+", "", haystack)
    compact_needle = re.sub(r"[^a-z0-9]+", "", needle)
    return len(compact_needle) >= 3 and compact_needle in compact_haystack


def _contains_any_term(text: str, terms: Iterable[str]) -> bool:
    return any(_text_contains_term(text, term) for term in terms if _clean(term))


def _has_relation_cue(text: str) -> bool:
    lowered = text.casefold()
    return any(cue in lowered for cue in SNIPPET_RELATION_CUES)


def _has_negation_or_contradiction(text: str) -> bool:
    lowered = text.casefold()
    return any(re.search(pattern, lowered) for pattern in SNIPPET_NEGATION_PATTERNS)


def _task_key(value: Any) -> str:
    key = _norm(value)
    return WEAK_FINAL_TASK_ALIASES.get(key, key)


def _split_allowed_tasks(values: Any) -> set[str]:
    terms = _split_terms(values)
    if not terms:
        return set(DEFAULT_ALLOW_WEAK_FINAL_TYPES)
    return {_task_key(term) for term in terms}


def _weak_final_cap_for_task(task_type: str) -> int:
    task = _task_key(task_type)
    if task == "shortentity":
        return 1
    if task == "entitylist":
        return 3
    return 0


def _task_allows_weak_final(task_type: str, allow_weak_final_types: Any) -> bool:
    task = _task_key(task_type)
    return bool(task and task in _split_allowed_tasks(allow_weak_final_types))


def _scope_matches(text: str, terms: list[str], mode: str = "any") -> bool:
    if not terms:
        return True
    haystack = text.casefold()
    if mode == "all":
        return all(term.casefold() in haystack for term in terms if term)
    if mode == "any":
        return any(term.casefold() in haystack for term in terms if term)
    raise ValueError("question_scope_match must be any or all")


def _matched_scope_terms(text: str, terms: list[str]) -> list[str]:
    if not terms:
        return []
    haystack = text.casefold()
    return [term for term in terms if term and term.casefold() in haystack]


def _serialized_scope_haystack(row: Mapping[str, Any]) -> str:
    fields = extract_evidence_fields(row)
    values = [
        fields.get(key, "")
        for key in [
            "source",
            "source_type",
            "relation",
            "target",
            "target_type",
            "info_source",
            "condition",
            "evidence_text",
            "original_relation",
        ]
    ]
    values.extend(str(row.get(key, "")) for key in ["boundary_reason", "_matched_targets"])
    return " ".join(value for value in values if value)


def _snippet_supported_weak(
    fields: Mapping[str, str],
    *,
    scope_terms: list[str],
    candidate_terms: list[str],
    candidate_side: str,
) -> tuple[bool, str]:
    snippet = " ".join(
        _clean(fields.get(key))
        for key in ["original_relation", "evidence_text", "condition"]
        if _clean(fields.get(key))
    )
    if not snippet or snippet.casefold() in {"undef", "['undef']", "[undef]"}:
        return False, "weak relation has no reviewable snippet"
    if _has_negation_or_contradiction(snippet):
        return False, "snippet contains negation or contradictory association language"
    if not _has_relation_cue(snippet):
        return False, "snippet lacks an explicit relation cue"

    source = fields.get("source", "")
    target = fields.get("target", "")
    candidate_name = target if candidate_side == "target" else source
    if candidate_side == "either":
        candidate_name = " ".join(x for x in [source, target] if x)
    anchor_terms = list(scope_terms)
    if candidate_side == "target":
        anchor_terms.append(source)
    elif candidate_side == "source":
        anchor_terms.append(target)
    else:
        anchor_terms.extend([source, target])

    candidate_term_list = list(candidate_terms) or [candidate_name]
    candidate_hit = _contains_any_term(snippet, candidate_term_list)
    anchor_hit = not anchor_terms or _contains_any_term(snippet, anchor_terms)
    if not candidate_hit:
        return False, "snippet does not mention the candidate answer/entity"
    if not anchor_hit:
        return False, "snippet does not mention the question anchor or scope"
    return True, "co_mentioned snippet directly mentions candidate, anchor/scope, and relation cue"


def _schema_sets(schema_file: Optional[str] = None) -> tuple[set[str], set[str], dict[str, Any]]:
    schema = load_kg_schema(schema_file)
    strong = {_norm(x) for x in schema.get("strong_curated_edge_types", DEFAULT_STRONG_CURATED)}
    weak = set(_norm(x) for x in schema.get("weak_trace_edge_types", []))
    if not weak:
        weak = set(_norm(x) for x in schema.get("weak_evidence_edge_types", DEFAULT_WEAK_TRACE))
        weak.update(_norm(x) for x in schema.get("blocked_query_edge_types", ["no_relation"]))
    return strong, weak, schema


def _query_stage(row: Mapping[str, Any]) -> str:
    return _row_get(row, ["query.stage", "query_stage", "selected_stage", "stage"])


def _candidate_stage_allowed(row: Mapping[str, Any], candidate_stage_policy: str) -> bool:
    if candidate_stage_policy == "any":
        return True
    if candidate_stage_policy != "target_retained_only":
        raise ValueError("candidate_stage_policy must be any or target_retained_only")
    stage = _query_stage(row)
    return not stage or stage in TARGET_RETAINED_CANDIDATE_STAGES


def _relation_priority(relation: Any) -> float:
    rel = _norm(relation)
    if rel in HIGH_PRIORITY_RELATIONS:
        return 1.0
    if rel in MEDIUM_PRIORITY_RELATIONS:
        return 0.65
    if rel in LOW_PRIORITY_RELATIONS:
        return 0.35
    return 0.20


def _source_reliability(info_source: Any) -> float:
    source = _clean(info_source).upper()
    if not source:
        return 0.35
    if any(token in source for token in HIGH_RELIABILITY_SOURCES):
        return 1.0
    if any(token in source for token in MEDIUM_RELIABILITY_SOURCES):
        return 0.65
    if "PUBMED" in source or "TEXT" in source or "CO_MENTION" in source:
        return 0.25
    return 0.45


def _stage_specificity(stage: Any) -> float:
    stage_key = _clean(stage)
    return STAGE_SPECIFICITY.get(stage_key, 0.60 if stage_key else 0.70)


def _score_option_packet(packet: Mapping[str, Any], required_targets: list[str]) -> dict[str, Any]:
    rows = list(packet.get("strong_curated", []))
    covered_targets = list(packet.get("diagnostics", {}).get("covered_targets", []))
    if not rows:
        return {
            "score": 0.0,
            "coverage": 0.0,
            "relation": 0.0,
            "source": 0.0,
            "stage": 0.0,
            "strong_rows": 0,
        }
    coverage = 1.0 if not required_targets else len({x.casefold() for x in covered_targets}) / max(1, len(required_targets))
    relation_score = max(_relation_priority(row.get("relation")) for row in rows)
    source_score = max(_source_reliability(row.get("info_source")) for row in rows)
    stage_score = max(_stage_specificity(row.get("query.stage") or row.get("query_stage")) for row in rows)
    score = (0.45 * coverage) + (0.25 * relation_score) + (0.15 * source_score) + (0.15 * stage_score)
    return {
        "score": round(score, 4),
        "coverage": round(coverage, 4),
        "relation": round(relation_score, 4),
        "source": round(source_score, 4),
        "stage": round(stage_score, 4),
        "strong_rows": len(rows),
    }


def extract_evidence_fields(row: Mapping[str, Any]) -> dict[str, str]:
    """Extract canonical fields from raw brick qr rows or evidence-unit rows."""
    info_source = _row_get(row, ["info_source", "path.1.info_source"])
    canonical_info_source = canonicalize_info_source(info_source)
    return {
        "source": _row_get(row, ["source", "path.0.name", "path.0.id"]),
        "source_type": _row_get(row, ["source_type", "path.0.type"]),
        "relation": _row_get(row, ["relation", "path.1.relation", "path.1.name", "path.1"]),
        "target": _row_get(row, ["target", "path.2.name", "path.2.id"]),
        "target_type": _row_get(row, ["target_type", "path.2.type"]),
        "info_source": canonical_info_source or info_source,
        "condition": _row_get(row, ["condition", "path.1.condition"]),
        "evidence_text": _row_get(row, ["evidence_text", "path.1.original_relation", "original_relation", "snippet"]),
        "original_relation": _row_get(row, ["original_relation", "path.1.original_relation", "evidence_text"]),
        "confidence": _row_get(row, ["confidence", "path.1.relation_confidence", "relation_confidence"]),
    }


def classify_evidence_row(
    row: Mapping[str, Any],
    expected_answer_type: Optional[str] = None,
    question_scope_terms: Optional[Iterable[str]] = None,
    question_scope_match: str = "any",
    candidate_terms: Optional[Iterable[str]] = None,
    candidate_side: str = "target",
    weak_evidence_policy: str = DEFAULT_WEAK_EVIDENCE_POLICY,
    schema_file: Optional[str] = None,
    candidate_stage_policy: str = "any",
    kg_relation_policy: str = DEFAULT_KG_RELATION_POLICY,
    task_type: str = "",
    allow_weak_final_types: Any = None,
) -> dict[str, Any]:
    """Classify one KG row into evidence tiers.

    ``candidate_side`` controls which endpoint is considered the answer
    candidate. MCQ option-wise verification usually uses ``source`` because
    the option is the hard anchor; list/factoid candidate ranking usually uses
    ``target``.
    """
    strong_relations, weak_relations, _schema = _schema_sets(schema_file)
    kg_policy = _normalize_kg_relation_policy(kg_relation_policy)
    weak_policy = _normalize_weak_policy(weak_evidence_policy)
    fields = extract_evidence_fields(row)
    relation = _norm(fields["relation"])

    if candidate_side not in {"source", "target", "either"}:
        raise ValueError("candidate_side must be source, target, or either")
    candidate_name = fields["target"] if candidate_side == "target" else fields["source"]
    candidate_type = fields["target_type"] if candidate_side == "target" else fields["source_type"]
    if candidate_side == "either":
        candidate_name = " ".join(x for x in [fields["source"], fields["target"]] if x)
        candidate_type = " ".join(x for x in [fields["source_type"], fields["target_type"]] if x)

    haystack = " ".join(
        fields.get(key, "")
        for key in ["source", "source_type", "relation", "target", "target_type", "info_source", "condition", "evidence_text", "original_relation"]
    )
    scope_terms = _split_terms(question_scope_terms)
    candidate_term_list = _split_terms(candidate_terms)
    expected_types = {_norm(x) for x in _split_terms(expected_answer_type)}

    type_match = not expected_types or any(_norm(piece) in expected_types for piece in re.split(r"[,;| ]+", candidate_type) if piece)
    scope_match = _scope_matches(haystack, scope_terms, question_scope_match)
    candidate_match = not candidate_term_list or _contains_any(candidate_name, candidate_term_list)
    stage_allowed = _candidate_stage_allowed(row, candidate_stage_policy)
    stage = _query_stage(row)

    is_weak = relation in weak_relations or "co_mention" in relation or relation == "mentioned"
    is_drop_relation = relation in DROP_RELATIONS
    is_strong_relation = relation in strong_relations

    if is_drop_relation:
        evidence_group = "drop"
        answer_role = "dropped"
        candidate_eligible = False
        reason = f"non-evidence relation: {fields['relation'] or 'missing'}"
    elif is_weak:
        snippet_supported, snippet_reason = _snippet_supported_weak(
            fields,
            scope_terms=scope_terms,
            candidate_terms=candidate_term_list,
            candidate_side=candidate_side,
        )
        if kg_policy in {"curated_only", "no_co_mentioned"} or weak_policy == "drop":
            evidence_group = "drop"
            answer_role = "dropped"
            candidate_eligible = False
            reason = f"weak relation dropped by policy: {fields['relation'] or 'missing'}"
        elif weak_policy == "snippet_supported" and snippet_supported:
            evidence_group = "snippet_supported_weak"
            answer_role = "support_evidence"
            candidate_eligible = False
            reason = snippet_reason
        elif weak_policy == "allow_weak_candidates" and snippet_supported:
            evidence_group = "snippet_supported_weak"
            answer_role = "support_evidence"
            candidate_eligible = False
            reason = snippet_reason
        else:
            evidence_group = "weak_trace"
            answer_role = "weak_trace"
            candidate_eligible = False
            reason = f"weak trace relation: {fields['relation'] or 'missing'}; {snippet_reason}"
    elif is_strong_relation and type_match and scope_match and candidate_match and stage_allowed:
        evidence_group = "strong_curated"
        answer_role = "candidate_answer"
        candidate_eligible = True
        reason = "schema-curated relation with matching answer boundary"
    else:
        evidence_group = "background"
        answer_role = "background_context"
        candidate_eligible = False
        if not is_strong_relation:
            reason = f"relation is not in strong_curated set: {fields['relation'] or 'missing'}"
        elif not type_match:
            reason = f"candidate type {candidate_type or 'missing'} does not match expected answer type"
        elif not scope_match:
            reason = "row did not match question scope terms"
        else:
            reason = "row did not match candidate anchor" if not candidate_match else f"query stage is not candidate-eligible: {stage or 'missing'}"

    return {
        "evidence_group": evidence_group,
        "answer_role": answer_role,
        "candidate_eligible": candidate_eligible,
        "evidence_tier": {
            "strong_curated": "E3",
            "snippet_supported_weak": "E2",
            "weak_trace": "E1",
            "background": "E0",
            "drop": "E0",
        }.get(evidence_group, "E0"),
        "candidate_side": candidate_side,
        "candidate_name": candidate_name,
        "candidate_type": candidate_type,
        "relation": fields["relation"],
        "query_stage": stage,
        "candidate_stage_policy": candidate_stage_policy,
        "candidate_stage_eligible": stage_allowed,
        "kg_relation_policy": kg_policy,
        "weak_evidence_policy": weak_policy,
        "task_type": task_type,
        "weak_final_allowed_task": _task_allows_weak_final(task_type, allow_weak_final_types),
        "reason": reason,
    }


def _serialize_row(row: Mapping[str, Any], classification: Mapping[str, Any]) -> dict[str, Any]:
    fields = extract_evidence_fields(row)
    out = {
        "source": fields["source"],
        "source_type": fields["source_type"],
        "relation": fields["relation"],
        "target": fields["target"],
        "target_type": fields["target_type"],
        "info_source": fields["info_source"],
        "condition": fields["condition"],
        "evidence_text": fields["evidence_text"],
        "confidence": fields["confidence"],
        "evidence_group": classification["evidence_group"],
        "answer_role": classification["answer_role"],
        "evidence_tier": classification.get("evidence_tier", ""),
        "candidate_eligible": classification["candidate_eligible"],
        "candidate_name": classification["candidate_name"],
        "candidate_type": classification["candidate_type"],
        "query_stage": classification.get("query_stage", ""),
        "candidate_stage_policy": classification.get("candidate_stage_policy", "any"),
        "candidate_stage_eligible": classification.get("candidate_stage_eligible", True),
        "kg_relation_policy": classification.get("kg_relation_policy", DEFAULT_KG_RELATION_POLICY),
        "weak_evidence_policy": classification.get("weak_evidence_policy", DEFAULT_WEAK_EVIDENCE_POLICY),
        "boundary_reason": classification["reason"],
    }
    for key in [
        "evidence_id",
        "path_id",
        "edge_position",
        "edge_support_class",
        "evidence_strength",
        "evidence_score",
        "query.stage",
        "query.target_query_mode",
        "query.source_stage",
        "query.info_source_filter",
    ]:
        if key in row and _clean(row.get(key)):
            out[key] = _clean(row.get(key))
    return out


def build_evidence_packet(
    rows: Any,
    expected_answer_type: Optional[str] = None,
    question_scope_terms: Optional[Iterable[str]] = None,
    question_scope_match: str = "any",
    candidate_terms: Optional[Iterable[str]] = None,
    candidate_side: str = "target",
    weak_evidence_policy: str = DEFAULT_WEAK_EVIDENCE_POLICY,
    max_candidate_rows: Optional[int] = None,
    schema_file: Optional[str] = None,
    selected_stage: str = "",
    target_query_mode: str = "",
    relation_alias_map: Optional[Mapping[str, Any]] = None,
    candidate_stage_policy: str = "any",
    kg_relation_policy: str = DEFAULT_KG_RELATION_POLICY,
    task_type: str = "",
    allow_weak_final_types: Any = None,
    co_mentioned_max_share: float = 0.15,
    co_mentioned_max_rows_per_question: int = 3,
) -> dict[str, Any]:
    """Build a structured packet from KG rows.

    The packet is safe to pass to a prompt renderer because rows are separated
    into answer-eligible strong evidence and non-answer weak/background context.
    """
    row_list = _rows_iter(rows)
    weak_policy = _normalize_weak_policy(weak_evidence_policy)
    kg_policy = _normalize_kg_relation_policy(kg_relation_policy)
    packet = {
        "strong_curated": [],
        "snippet_supported_weak": [],
        "weak_trace": [],
        "background": [],
        "drop": [],
        "selected_stage": selected_stage,
        "target_query_mode": target_query_mode,
        "relation_alias_map": dict(relation_alias_map or {}),
        "diagnostics": {
            "row_count": len(row_list),
            "kg_raw_rows": len(row_list),
            "strong_curated_rows": 0,
            "snippet_supported_weak_rows": 0,
            "weak_trace_rows": 0,
            "background_rows": 0,
            "drop_rows": 0,
            "candidate_eligible_rows": 0,
            "kg_relation_policy": kg_policy,
            "weak_evidence_policy": weak_policy,
            "candidate_side": candidate_side,
            "expected_answer_type": expected_answer_type or "",
            "candidate_stage_policy": candidate_stage_policy,
            "task_type": task_type,
            "allow_weak_final_types": sorted(_split_allowed_tasks(allow_weak_final_types)),
            "co_mentioned_max_share": float(co_mentioned_max_share),
            "co_mentioned_max_rows_per_question": int(co_mentioned_max_rows_per_question),
        },
    }
    for row in row_list:
        classification = classify_evidence_row(
            row,
            expected_answer_type=expected_answer_type,
            question_scope_terms=question_scope_terms,
            question_scope_match=question_scope_match,
            candidate_terms=candidate_terms,
            candidate_side=candidate_side,
            weak_evidence_policy=weak_policy,
            schema_file=schema_file,
            candidate_stage_policy=candidate_stage_policy,
            kg_relation_policy=kg_policy,
            task_type=task_type,
            allow_weak_final_types=allow_weak_final_types,
        )
        group = classification["evidence_group"]
        serialized = _serialize_row(row, classification)
        if group == "strong_curated" and max_candidate_rows is not None:
            current_candidates = sum(1 for item in packet["strong_curated"] if item.get("candidate_eligible"))
            if current_candidates >= max_candidate_rows:
                serialized["candidate_eligible"] = False
                serialized["answer_role"] = "support_evidence"
                serialized["boundary_reason"] = "candidate cap reached; kept as support evidence"
        packet[group].append(serialized)

    weak_final_allowed = _task_allows_weak_final(task_type, allow_weak_final_types)
    weak_final_cap = _weak_final_cap_for_task(task_type)
    weak_final_cap = min(weak_final_cap, max(0, int(co_mentioned_max_rows_per_question)))
    strong_candidate_count = sum(1 for item in packet["strong_curated"] if item.get("candidate_eligible"))
    if strong_candidate_count > 0 and co_mentioned_max_share >= 0:
        share_cap = int((float(co_mentioned_max_share) * strong_candidate_count) / max(1e-9, 1.0 - float(co_mentioned_max_share)))
        weak_final_cap = min(weak_final_cap, max(0, share_cap))
    if not weak_final_allowed:
        weak_final_cap = 0
    snippet_rows = sorted(
        packet["snippet_supported_weak"],
        key=lambda item: (
            -_source_reliability(item.get("info_source")),
            -_relation_priority(item.get("relation")),
            str(item.get("candidate_name", "")).casefold(),
        ),
    )
    allowed_ids = {id(item) for item in snippet_rows[:weak_final_cap]}
    for row in packet["snippet_supported_weak"]:
        if id(row) in allowed_ids:
            row["candidate_eligible"] = True
            row["answer_role"] = "candidate_answer"
            row["boundary_reason"] = f"{row.get('boundary_reason', '')}; snippet-supported weak evidence allowed by capped {task_type or 'task'} policy"
        else:
            row["candidate_eligible"] = False
            row["answer_role"] = "support_evidence"
            if weak_final_allowed:
                row["boundary_reason"] = f"{row.get('boundary_reason', '')}; snippet-supported weak evidence kept as support after cap"

    for group in ["strong_curated", "snippet_supported_weak", "weak_trace", "background", "drop"]:
        packet["diagnostics"][f"{group}_rows"] = len(packet[group])
    packet["diagnostics"]["background_reason_counts"] = dict(
        Counter(row.get("boundary_reason", "") for row in packet["background"] if row.get("boundary_reason"))
    )
    packet["diagnostics"]["drop_reason_counts"] = dict(
        Counter(row.get("boundary_reason", "") for row in packet["drop"] if row.get("boundary_reason"))
    )
    packet["diagnostics"]["candidate_eligible_rows"] = sum(
        1
        for group in ["strong_curated", "snippet_supported_weak"]
        for item in packet[group]
        if item.get("candidate_eligible")
    )
    packet["diagnostics"]["has_strong_evidence"] = bool(packet["strong_curated"])
    packet["diagnostics"]["has_snippet_supported_weak"] = bool(packet["snippet_supported_weak"])
    packet["diagnostics"]["has_weak_only_evidence"] = (
        bool(packet["weak_trace"] or packet["snippet_supported_weak"]) and not packet["strong_curated"]
    )
    co_raw = sum(1 for row in row_list if _norm(extract_evidence_fields(row).get("relation")) == "co_mentioned")
    co_final = sum(
        1
        for item in packet["snippet_supported_weak"]
        if item.get("candidate_eligible") and _norm(item.get("relation")) == "co_mentioned"
    )
    used_rows = packet["strong_curated"] + [row for row in packet["snippet_supported_weak"] if row.get("candidate_eligible")]
    packet["diagnostics"].update(
        {
            "kg_used_rows": len(used_rows),
            "kg_final_answer_sources": sum(1 for row in used_rows if row.get("candidate_eligible")),
            "co_mentioned_raw_rows": co_raw,
            "co_mentioned_used_rows": co_final,
            "co_mentioned_final_sources": co_final,
            "ordinary_co_mentioned_final_sources": 0,
            "co_mentioned_raw_share": round(co_raw / len(row_list), 4) if row_list else 0.0,
            "co_mentioned_used_share": round(co_final / len(used_rows), 4) if used_rows else 0.0,
            "co_mentioned_final_source_share": round(co_final / max(1, sum(1 for row in used_rows if row.get("candidate_eligible"))), 4)
            if used_rows
            else 0.0,
            "evidence_tier_counts": dict(
                Counter(
                    row.get("evidence_tier", "E0")
                    for group in ["strong_curated", "snippet_supported_weak", "weak_trace", "background", "drop"]
                    for row in packet[group]
                )
            ),
            "answer_role_counts": dict(
                Counter(
                    row.get("answer_role", "")
                    for group in ["strong_curated", "snippet_supported_weak", "weak_trace", "background", "drop"]
                    for row in packet[group]
                    if row.get("answer_role")
                )
            ),
        }
    )
    return packet


def _contract_dict(contract: Any) -> dict[str, Any]:
    if TaskContract is not None and isinstance(contract, TaskContract):
        return contract.to_dict()
    if isinstance(contract, Mapping):
        return dict(contract)
    return {
        "task_type": getattr(contract, "task_type", ""),
        "question": getattr(contract, "question", ""),
        "expected_answer_type": getattr(contract, "expected_answer_type", None),
        "source_mentions": getattr(contract, "source_mentions", []),
        "target_mentions": getattr(contract, "target_mentions", []),
        "options": getattr(contract, "options", []),
        "evidence_policy": getattr(contract, "evidence_policy", {}),
        "metadata": getattr(contract, "metadata", {}),
    }


def build_evidence_packet_from_contract(
    rows: Any,
    contract: Any,
    *,
    schema_file: Optional[str] = None,
    selected_stage: str = "",
    target_query_mode: str = "",
    relation_alias_map: Optional[Mapping[str, Any]] = None,
    **overrides: Any,
) -> dict[str, Any]:
    """Build a packet using a dataset-agnostic TaskContract.

    Runners should prefer this wrapper when they already have a contract. It
    keeps weak-evidence and answer-boundary defaults in BRICK instead of
    spreading them across benchmark-specific code.
    """
    data = _contract_dict(contract)
    policy = dict(data.get("evidence_policy") or {})
    policy.update({key: value for key, value in overrides.items() if value is not None})
    task_type = data.get("task_type", "")
    weak_policy = policy.get("weak_evidence_policy", DEFAULT_WEAK_EVIDENCE_POLICY)
    kg_relation_policy = policy.get("kg_relation_policy", DEFAULT_KG_RELATION_POLICY)
    allow_weak_final_types = policy.get("allow_weak_final_types")
    expected_answer_type = data.get("expected_answer_type")
    target_mentions = data.get("target_mentions") or []
    source_mentions = data.get("source_mentions") or []

    if task_type == "ClosedChoice" and isinstance(rows, Mapping):
        return build_optionwise_evidence_packets(
            rows,
            expected_answer_type=expected_answer_type,
            question_scope_terms=target_mentions,
            question_scope_match=policy.get("question_scope_match") or ("all" if len(target_mentions) > 1 else "any"),
            weak_evidence_policy=weak_policy,
            schema_file=schema_file,
            scope_aggregation=policy.get("scope_aggregation", "row"),
            candidate_stage_policy=policy.get("candidate_stage_policy", "any"),
            evidence_ranking=policy.get("evidence_ranking", "off") in {True, "on", "true", "1"},
            evidence_margin=float(policy.get("evidence_margin") or 0.0),
            kg_relation_policy=kg_relation_policy,
            allow_weak_final_types=allow_weak_final_types,
            co_mentioned_max_share=float(policy.get("co_mentioned_max_share", 0.15)),
            co_mentioned_max_rows_per_question=int(policy.get("co_mentioned_max_rows_per_question", 3)),
        )

    candidate_side = policy.get("candidate_side")
    if not candidate_side:
        candidate_side = "source" if task_type in {"ClosedChoice", "BinaryClaim"} else "target"

    return build_evidence_packet(
        rows,
        expected_answer_type=expected_answer_type,
        question_scope_terms=policy.get("question_scope_terms", target_mentions),
        question_scope_match=policy.get("question_scope_match", "any"),
        candidate_terms=policy.get("candidate_terms", source_mentions if task_type == "BinaryClaim" else None),
        candidate_side=candidate_side,
        weak_evidence_policy=weak_policy,
        max_candidate_rows=policy.get("max_candidate_rows"),
        schema_file=schema_file,
        selected_stage=selected_stage,
        target_query_mode=target_query_mode,
        relation_alias_map=relation_alias_map,
        candidate_stage_policy=policy.get("candidate_stage_policy", "any"),
        kg_relation_policy=kg_relation_policy,
        task_type=task_type,
        allow_weak_final_types=allow_weak_final_types,
        co_mentioned_max_share=float(policy.get("co_mentioned_max_share", 0.15)),
        co_mentioned_max_rows_per_question=int(policy.get("co_mentioned_max_rows_per_question", 3)),
    )


def build_optionwise_evidence_packets(
    option_rows: Mapping[str, Any],
    expected_answer_type: Optional[str] = None,
    question_scope_terms: Optional[Iterable[str]] = None,
    question_scope_match: str = "any",
    weak_evidence_policy: str = "context_only",
    schema_file: Optional[str] = None,
    scope_aggregation: str = "option",
    candidate_stage_policy: str = "any",
    evidence_ranking: bool = False,
    evidence_margin: float = 0.0,
    kg_relation_policy: str = DEFAULT_KG_RELATION_POLICY,
    allow_weak_final_types: Any = None,
    co_mentioned_max_share: float = 0.15,
    co_mentioned_max_rows_per_question: int = 3,
) -> dict[str, Any]:
    """Build MCQ-style packets where each option is a hard source anchor."""
    if question_scope_match not in {"any", "all"}:
        raise ValueError("question_scope_match must be any or all")
    if scope_aggregation not in {"row", "option"}:
        raise ValueError("scope_aggregation must be row or option")
    if candidate_stage_policy not in {"any", "target_retained_only"}:
        raise ValueError("candidate_stage_policy must be any or target_retained_only")
    required_targets = _split_terms(question_scope_terms)
    row_scope_match = "any" if question_scope_match == "all" and scope_aggregation == "option" else question_scope_match
    packets = {}
    option_strong_counts = {}
    option_weak_counts = {}
    option_background_counts = {}
    option_partial_strong_counts = {}
    option_scope_pass = {}
    option_required_targets = {}
    option_covered_targets = {}
    option_missing_targets = {}
    for option, rows in option_rows.items():
        packet = build_evidence_packet(
            rows,
            expected_answer_type=expected_answer_type,
            question_scope_terms=question_scope_terms,
            question_scope_match=row_scope_match,
            candidate_terms=[option],
            candidate_side="source",
            weak_evidence_policy=weak_evidence_policy,
            schema_file=schema_file,
            candidate_stage_policy=candidate_stage_policy,
            kg_relation_policy=kg_relation_policy,
            task_type="ClosedChoice",
            allow_weak_final_types=allow_weak_final_types,
            co_mentioned_max_share=co_mentioned_max_share,
            co_mentioned_max_rows_per_question=co_mentioned_max_rows_per_question,
        )
        strong_rows = list(packet.get("strong_curated", []))
        covered_targets: list[str] = []
        for row in strong_rows:
            matched = _matched_scope_terms(_serialized_scope_haystack(row), required_targets)
            row["_matched_targets"] = matched
            for target in matched:
                if target not in covered_targets:
                    covered_targets.append(target)

        if scope_aggregation == "row":
            scope_pass = bool(strong_rows)
        elif not required_targets:
            scope_pass = bool(strong_rows)
        elif question_scope_match == "all":
            covered_keys = {target.casefold() for target in covered_targets}
            scope_pass = all(target.casefold() in covered_keys for target in required_targets)
        else:
            scope_pass = bool(covered_targets)

        partial_strong: list[dict[str, Any]] = []
        if scope_aggregation == "option" and not scope_pass and strong_rows:
            covered_keys = {target.casefold() for target in covered_targets}
            missing_targets = [target for target in required_targets if target.casefold() not in covered_keys]
            partial_reason = (
                "option-level scope incomplete"
                + (f"; missing targets: {', '.join(missing_targets)}" if missing_targets else "")
            )
            for row in strong_rows:
                partial = dict(row)
                partial["candidate_eligible"] = False
                partial["answer_role"] = "support_evidence"
                partial["evidence_group"] = "partial_strong"
                partial["boundary_reason"] = partial_reason
                partial_strong.append(partial)
            packet["strong_curated"] = []

        packet["partial_strong"] = partial_strong
        covered_keys = {target.casefold() for target in covered_targets}
        missing_targets = [target for target in required_targets if target.casefold() not in covered_keys]
        packet["diagnostics"].update(
            {
                "scope_match_mode": question_scope_match,
                "row_scope_match_mode": row_scope_match,
                "scope_aggregation": scope_aggregation,
                "required_targets": list(required_targets),
                "covered_targets": list(covered_targets),
                "missing_targets": missing_targets,
                "option_scope_pass": scope_pass,
                "partial_strong_rows": len(partial_strong),
                "strong_curated_rows": len(packet["strong_curated"]),
                "candidate_eligible_rows": sum(1 for item in packet["strong_curated"] if item.get("candidate_eligible")),
                "has_strong_evidence": bool(packet["strong_curated"]),
                "snippet_supported_weak_rows": len(packet["snippet_supported_weak"]),
                "drop_rows": len(packet["drop"]),
                "has_weak_only_evidence": bool(packet["weak_trace"] or packet["snippet_supported_weak"]) and not packet["strong_curated"],
            }
        )
        packets[str(option)] = packet
        option_strong_counts[str(option)] = packet["diagnostics"]["strong_curated_rows"]
        option_weak_counts[str(option)] = packet["diagnostics"]["weak_trace_rows"] + packet["diagnostics"]["snippet_supported_weak_rows"]
        option_background_counts[str(option)] = packet["diagnostics"]["background_rows"]
        option_partial_strong_counts[str(option)] = packet["diagnostics"]["partial_strong_rows"]
        option_scope_pass[str(option)] = scope_pass
        option_required_targets[str(option)] = list(required_targets)
        option_covered_targets[str(option)] = list(covered_targets)
        option_missing_targets[str(option)] = list(missing_targets)
    has_strong = any(count > 0 for count in option_strong_counts.values())
    option_evidence_scores = {
        option: _score_option_packet(packet, required_targets)
        for option, packet in packets.items()
    }
    ranked_options = sorted(
        [
            {
                "option": option,
                **score,
            }
            for option, score in option_evidence_scores.items()
            if option_strong_counts.get(option, 0) > 0
        ],
        key=lambda item: (-float(item.get("score", 0.0)), str(item.get("option", ""))),
    )
    best_option = str(ranked_options[0]["option"]) if ranked_options else ""
    second_option = str(ranked_options[1]["option"]) if len(ranked_options) > 1 else ""
    best_score = float(ranked_options[0]["score"]) if ranked_options else 0.0
    second_score = float(ranked_options[1]["score"]) if len(ranked_options) > 1 else 0.0
    margin = round(best_score - second_score, 4) if ranked_options else 0.0
    margin_pass = bool(has_strong and (not evidence_ranking or len(ranked_options) == 1 or margin >= evidence_margin))
    prompt_options = [option for option, count in option_strong_counts.items() if count > 0]
    if evidence_ranking:
        prompt_options = [best_option] if margin_pass and best_option else []
    fallback_native = not has_strong or (evidence_ranking and not margin_pass)
    return {
        "mode": "mcq_direct_only",
        "weak_evidence_policy": weak_evidence_policy,
        "kg_relation_policy": kg_relation_policy,
        "scope_aggregation": scope_aggregation,
        "candidate_stage_policy": candidate_stage_policy,
        "evidence_ranking": bool(evidence_ranking),
        "evidence_margin": float(evidence_margin),
        "packets": packets,
        "option_strong_counts": option_strong_counts,
        "option_weak_counts": option_weak_counts,
        "option_background_counts": option_background_counts,
        "option_partial_strong_counts": option_partial_strong_counts,
        "option_scope_pass": option_scope_pass,
        "option_required_targets": option_required_targets,
        "option_covered_targets": option_covered_targets,
        "option_missing_targets": option_missing_targets,
        "option_evidence_scores": option_evidence_scores,
        "option_evidence_rank": ranked_options,
        "best_option": best_option,
        "second_option": second_option,
        "evidence_score_margin": margin,
        "evidence_margin_pass": margin_pass,
        "prompt_options": prompt_options,
        "fallback_native": fallback_native,
        "prompt_groups": ["strong_curated"] if not fallback_native else [],
        "diagnostics": {
            "option_count": len(option_rows),
            "options_with_strong_evidence": sum(1 for count in option_strong_counts.values() if count > 0),
            "options_with_weak_trace": sum(1 for count in option_weak_counts.values() if count > 0),
            "options_with_partial_strong": sum(1 for count in option_partial_strong_counts.values() if count > 0),
            "scope_match_mode": question_scope_match,
            "scope_aggregation": scope_aggregation,
            "candidate_stage_policy": candidate_stage_policy,
            "evidence_ranking": bool(evidence_ranking),
            "evidence_margin": float(evidence_margin),
            "evidence_margin_pass": margin_pass,
        },
    }


def packet_to_rows(packet: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for group in ["strong_curated", "partial_strong", "snippet_supported_weak", "weak_trace", "background", "drop"]:
        for row in packet.get(group, []):
            out = dict(row)
            out["packet_group"] = group
            rows.append(out)
    return rows


def render_prompt_evidence(packet: Mapping[str, Any], max_rows: int = 12) -> str:
    """Render answer-eligible evidence for a model prompt.

    E3 strong rows are always prompt-eligible. E2 snippet-supported weak rows
    are rendered only when their packet-level caps marked them candidate
    eligible; ordinary co_mentioned rows are never rendered here.
    """
    lines = []
    prompt_rows = list(packet.get("strong_curated", []))
    prompt_rows.extend(row for row in packet.get("snippet_supported_weak", []) if row.get("candidate_eligible"))
    for row in prompt_rows[:max_rows]:
        tier = row.get("evidence_tier") or "E3"
        lines.append(
            f"- [{tier}] {row.get('source')} --{row.get('relation')}--> {row.get('target')}"
            f" | source={row.get('info_source') or 'NA'}"
            f" | text={row.get('evidence_text') or row.get('condition') or 'NA'}"
        )
    return "\n".join(lines)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build KG evidence packets with answer-boundary roles")
    parser.add_argument("-i", "--input", required=True, help="Evidence-unit or brick qr CSV")
    parser.add_argument("--expected-answer-type")
    parser.add_argument("--scope-term", action="append", help="Question scope term; may be repeated")
    parser.add_argument("--scope-match", choices=["any", "all"], default="any")
    parser.add_argument("--candidate-term", action="append", help="Candidate anchor term; may be repeated")
    parser.add_argument("--candidate-side", choices=["source", "target", "either"], default="target")
    parser.add_argument(
        "--kg-relation-policy",
        choices=sorted(KG_RELATION_POLICIES),
        default="curated_preferred",
    )
    parser.add_argument(
        "--weak-evidence-policy",
        choices=sorted(WEAK_EVIDENCE_POLICIES),
        default=DEFAULT_WEAK_EVIDENCE_POLICY,
    )
    parser.add_argument("--allow-weak-final-types", default="factoid,list")
    parser.add_argument("--co-mentioned-max-share", type=float, default=0.15)
    parser.add_argument("--co-mentioned-max-rows-per-question", type=int, default=3)
    parser.add_argument("--task-type", default="")
    parser.add_argument("--candidate-stage-policy", choices=["any", "target_retained_only"], default="any")
    parser.add_argument("--max-candidate-rows", type=int)
    parser.add_argument("--kg-schema")
    parser.add_argument("--format", choices=["json", "csv", "prompt"], default="json")
    parser.add_argument("-o", "--output")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if pd is None:
        raise RuntimeError("pandas is required for evidence packet CLI")
    df = pd.read_csv(args.input)
    packet = build_evidence_packet(
        df,
        expected_answer_type=args.expected_answer_type,
        question_scope_terms=args.scope_term,
        question_scope_match=args.scope_match,
        candidate_terms=args.candidate_term,
        candidate_side=args.candidate_side,
        weak_evidence_policy=args.weak_evidence_policy,
        max_candidate_rows=args.max_candidate_rows,
        schema_file=args.kg_schema,
        candidate_stage_policy=args.candidate_stage_policy,
        kg_relation_policy=args.kg_relation_policy,
        task_type=args.task_type,
        allow_weak_final_types=args.allow_weak_final_types,
        co_mentioned_max_share=args.co_mentioned_max_share,
        co_mentioned_max_rows_per_question=args.co_mentioned_max_rows_per_question,
    )
    if args.format == "prompt":
        text = render_prompt_evidence(packet)
    elif args.format == "csv":
        text = pd.DataFrame(packet_to_rows(packet)).to_csv(index=False)
    else:
        text = json.dumps(packet, ensure_ascii=False, indent=2)
    if args.output:
        Path(args.output).write_text(text + ("\n" if not text.endswith("\n") else ""), encoding="utf-8")
        print(f"Wrote {args.output}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
