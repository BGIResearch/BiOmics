"""brick reason - offline reasoning helpers over brick query/rank outputs."""

import argparse
import ast
import hashlib
import json
import math
import re
import sqlite3
import sys
from difflib import SequenceMatcher
from pathlib import Path

import networkx as nx
import pandas as pd

try:
    from kgpolicy import source_quality_tier, source_quality_weight
except Exception:  # pragma: no cover
    def source_quality_tier(info_source):
        return "unknown"
    def source_quality_weight(info_source):
        return 0.55

try:
    from evidence_packet import build_evidence_packet, classify_evidence_row, packet_to_rows, render_prompt_evidence
except Exception:  # pragma: no cover
    build_evidence_packet = None
    classify_evidence_row = None
    packet_to_rows = None
    render_prompt_evidence = None

try:
    from task_contract import TaskContract
except Exception:  # pragma: no cover
    TaskContract = None


def _contract_value(contract, key, default=None):
    if isinstance(contract, dict):
        return contract.get(key, default)
    return getattr(contract, key, default)


def reason_closed_choice(contract, evidence_packet, native_answer=None):
    """Choose a closed-choice answer from a contract and optionwise packet.

    Evidence can override native reasoning only when the packet says it should
    be injected. Otherwise the native answer is preserved.
    """
    if not evidence_packet or evidence_packet.get("fallback_native"):
        return {"answer": native_answer or "", "reasoning_policy": "closed_choice_native_fallback"}
    if evidence_packet.get("evidence_ranking") is False:
        return {
            "answer": native_answer or "",
            "reasoning_policy": "closed_choice_llm_with_strong_evidence",
            "prompt_options": evidence_packet.get("prompt_options") or [],
        }
    prompt_options = list(evidence_packet.get("prompt_options") or [])
    best_option = evidence_packet.get("best_option") or (prompt_options[0] if prompt_options else "")
    if best_option:
        return {
            "answer": best_option,
            "reasoning_policy": "closed_choice_evidence_rank",
            "evidence_margin": evidence_packet.get("evidence_score_margin"),
            "prompt_options": prompt_options,
        }
    return {"answer": native_answer or "", "reasoning_policy": "closed_choice_native_fallback"}


WEAK_CLAIM_RELATIONS = {"", "co_mentioned", "mentioned", "no_relation"}
ENTITY_TYPE_LABELS = {
    "gene": "Gene",
    "protein": "Protein",
    "disease": "Disease",
    "phenotype": "Phenotype",
    "chemical": "Chemical",
    "drug": "Chemical",
    "pathway": "Pathway",
}


def _normalize_claim_text(value):
    text = "" if value is None else str(value)
    text = text.casefold().strip().strip("\"'`.,;:()[]{}")
    return re.sub(r"\s+", " ", text)


def _contract_policy(contract):
    policy = _contract_value(contract, "evidence_policy", {}) or {}
    return dict(policy) if isinstance(policy, dict) else {}


def _contract_schema_format(contract):
    schema = _contract_value(contract, "output_schema", {}) or {}
    if isinstance(schema, dict):
        return str(schema.get("format") or "")
    return ""


def _verdict_to_contract_answer(contract, verdict):
    fmt = _contract_schema_format(contract).casefold()
    if verdict not in {"true", "false"}:
        return "unknown"
    if fmt == "yes_no":
        return "yes" if verdict == "true" else "no"
    return verdict


def _split_labeled_entity(text, default_type=None):
    cleaned = re.sub(r"\s+", " ", text or "").strip(" .?;:,\"'`")
    match = re.match(r"^(gene|protein|disease|phenotype|chemical|drug|pathway)\s+(.+)$", cleaned, flags=re.I)
    if match:
        return match.group(2).strip(" .?;:,\"'`"), ENTITY_TYPE_LABELS.get(match.group(1).casefold())
    return cleaned, default_type


def _kg_query_for_claim(source, source_type=None, target=None, target_type=None, relation=None, limit=40):
    try:
        from querygraph import query_relation

        return query_relation(
            source_entity_set=source,
            target_entity_set=target,
            relation=relation,
            multi_hop=1,
            directed=False,
            source_entity_type=source_type,
            target_entity_type=target_type,
            limit_number=limit,
            return_type="dataframe",
        ), ""
    except Exception as exc:  # pragma: no cover - live KG failures are environment dependent
        return pd.DataFrame(), repr(exc)


def _claim_relation_counts(df):
    if df is None or df.empty:
        return {"rows": 0, "direct_rows": 0, "co_mentioned_rows": 0, "relations": {}}
    rel_col = "path.1.relation" if "path.1.relation" in df.columns else "relation"
    if rel_col not in df.columns:
        return {"rows": int(len(df)), "direct_rows": 0, "co_mentioned_rows": 0, "relations": {}}
    rel = df[rel_col].fillna("").astype(str)
    weak = rel.str.casefold().isin(WEAK_CLAIM_RELATIONS)
    return {
        "rows": int(len(df)),
        "direct_rows": int((~weak).sum()),
        "co_mentioned_rows": int((rel == "co_mentioned").sum()),
        "relations": {str(k): int(v) for k, v in rel.value_counts().to_dict().items()},
    }


def _claim_df_text(df, max_rows=80):
    if df is None or df.empty:
        return ""
    cols = [
        col
        for col in df.columns
        if any(token in col.casefold() for token in ("name", "id", "xref", "synonym", "relation", "original_relation", "condition"))
    ]
    if not cols:
        cols = list(df.columns[:12])
    values = []
    for _, row in df.head(max_rows).iterrows():
        for col in cols:
            value = row.get(col)
            try:
                missing = bool(pd.isna(value))
            except Exception:
                missing = False
            if value is not None and not missing:
                values.append(str(value))
    return "\n".join(values)


def _claim_text_mentions(text, term):
    needle = _normalize_claim_text(term)
    return bool(needle) and needle in _normalize_claim_text(text)


def _claim_verifier_result(verdict, template, stats, error="", **extra):
    out = {
        "verdict": verdict,
        "claim_template": template,
        "fallback_native": verdict == "unknown",
        "selected_stage": "kg_proposition_verifier" if verdict != "unknown" else "kg_proposition_unknown",
        "error": error,
        **stats,
    }
    out.update(extra)
    return out


def verify_binary_claim_with_kg(question, kg_limit=40):
    """Verify simple binary biomedical claims with direct KG evidence only.

    Returns a dict with verdict=true/false/unknown. co_mentioned/mentioned/no_relation
    rows are never treated as truth evidence.
    """
    text = re.sub(r"\s+", " ", question or "").strip(" .")

    doid = re.match(r"^Disease ontology identifier for\s+(.+?)\s+is\s+(DOID:\d+)$", text, flags=re.I)
    if doid:
        disease, claimed = doid.group(1).strip(), doid.group(2).upper()
        df, error = _kg_query_for_claim(disease, source_type="Disease", limit=min(kg_limit, 40))
        stats = _claim_relation_counts(df)
        haystack = _claim_df_text(df).upper()
        if claimed in haystack:
            return _claim_verifier_result("true", "doid", stats, error, claim_subject=disease, claimed_doid=claimed)
        if df is not None and not df.empty:
            return _claim_verifier_result("false", "doid", stats, error, claim_subject=disease, claimed_doid=claimed)
        return _claim_verifier_result("unknown", "doid", stats, error, claim_subject=disease, claimed_doid=claimed)

    templates = [
        {
            "name": "associated_with",
            "pattern": r"^(?P<subject>.+?)\s+is\s+(?P<neg>not\s+)?associated\s+with\s+(?P<object>.+?)$",
            "relations": None,
            "source_default": None,
            "target_default": None,
        },
        {
            "name": "is_a",
            "pattern": r"^(?P<subject>.+?)\s+is\s+(?P<neg>not\s+)?(?:a|an)?\s*(?P<object>.+?)$",
            "relations": ["is_a"],
            "source_default": "Disease",
            "target_default": "Disease",
        },
    ]
    for template in templates:
        match = re.match(template["pattern"], text, flags=re.I)
        if not match:
            continue
        source, source_type = _split_labeled_entity(match.group("subject"), template["source_default"])
        target, target_type = _split_labeled_entity(match.group("object"), template["target_default"])
        if not source or not target or len(source.split()) > 14 or len(target.split()) > 14:
            return _claim_verifier_result("unknown", template["name"], {"rows": 0, "direct_rows": 0, "co_mentioned_rows": 0, "relations": {}}, claim_subject=source, claim_target=target)
        df, error = _kg_query_for_claim(
            source,
            source_type=source_type,
            target=target,
            target_type=target_type,
            relation=template["relations"],
            limit=min(kg_limit, 50),
        )
        stats = _claim_relation_counts(df)
        if df is None or df.empty or stats.get("direct_rows", 0) <= 0:
            return _claim_verifier_result("unknown", template["name"], stats, error, claim_subject=source, claim_target=target)
        if not _claim_text_mentions(_claim_df_text(df), target):
            return _claim_verifier_result("unknown", template["name"], stats, error, claim_subject=source, claim_target=target)
        negated = bool(match.groupdict().get("neg"))
        verdict = "false" if negated else "true"
        return _claim_verifier_result(
            verdict,
            template["name"],
            stats,
            error,
            claim_subject=source,
            claim_subject_type=source_type or "",
            claim_target=target,
            claim_target_type=target_type or "",
            claim_negated=negated,
        )
    return _claim_verifier_result("unknown", "unparsed", {"rows": 0, "direct_rows": 0, "co_mentioned_rows": 0, "relations": {}})


def reason_binary_claim(contract, evidence_packet, native_answer=None):
    """Return yes/no/unknown-style claim reasoning from evidence roles."""
    evidence_packet = evidence_packet or {}
    explicit_verdict = evidence_packet.get("binary_claim_verdict")
    if isinstance(explicit_verdict, dict) and explicit_verdict.get("verdict") in {"true", "false"}:
        verdict = explicit_verdict.get("verdict")
        return {
            "answer": _verdict_to_contract_answer(contract, verdict),
            "reasoning_policy": "binary_claim_kg_proposition_verifier",
            "verdict": verdict,
            "verifier": explicit_verdict,
        }

    policy = _contract_policy(contract)
    if (
        policy.get("solver") == "kg_proposition_verifier"
        and not evidence_packet.get("kg_verifier_checked")
    ):
        verification = verify_binary_claim_with_kg(_contract_value(contract, "question", ""))
        if verification.get("verdict") in {"true", "false"}:
            verdict = verification.get("verdict")
            return {
                "answer": _verdict_to_contract_answer(contract, verdict),
                "reasoning_policy": "binary_claim_kg_proposition_verifier",
                "verdict": verdict,
                "verifier": verification,
            }
        if native_answer:
            return {
                "answer": native_answer,
                "reasoning_policy": "binary_claim_kg_unknown_fallback",
                "verdict": "unknown",
                "verifier": verification,
            }

    strong = list((evidence_packet or {}).get("strong_curated", []))
    weak = list((evidence_packet or {}).get("weak_trace", []))
    if strong:
        return {"answer": native_answer or "yes", "reasoning_policy": "binary_claim_supported"}
    if weak:
        return {"answer": native_answer or "unknown", "reasoning_policy": "binary_claim_weak_trace_only"}
    return {"answer": native_answer or "unknown", "reasoning_policy": "binary_claim_insufficient"}


def reason_short_entity(contract, evidence_packet, native_answer=None):
    """Use strong KG evidence for canonicalization, not forced invention."""
    strong = list((evidence_packet or {}).get("strong_curated", []))
    if strong:
        candidate = strong[0].get("candidate_name") or strong[0].get("target") or native_answer or ""
        return {"answer": candidate, "reasoning_policy": "short_entity_strong_canonical"}
    snippet = [
        row
        for row in list((evidence_packet or {}).get("snippet_supported_weak", []))
        if row.get("candidate_eligible")
    ]
    if snippet:
        candidate = snippet[0].get("candidate_name") or snippet[0].get("target") or native_answer or ""
        return {"answer": candidate, "reasoning_policy": "short_entity_snippet_supported_weak_canonical"}
    return {"answer": native_answer or "", "reasoning_policy": "short_entity_native_fallback"}


def reason_entity_list(contract, evidence_packet, native_answer=None, max_kg_added_answers=3):
    """Direct candidates first; KG may add only capped strong candidates."""
    answers = []
    if isinstance(native_answer, list):
        answers.extend(str(item).strip() for item in native_answer if str(item).strip())
    elif native_answer:
        answers.extend(part.strip() for part in str(native_answer).split(",") if part.strip())
    kg_rows = list((evidence_packet or {}).get("strong_curated", []))
    kg_rows.extend(
        row
        for row in list((evidence_packet or {}).get("snippet_supported_weak", []))
        if row.get("candidate_eligible")
    )
    for row in kg_rows[:max_kg_added_answers]:
        candidate = row.get("candidate_name") or row.get("target")
        if candidate and str(candidate).strip() not in answers:
            answers.append(str(candidate).strip())
    return {"answer": answers, "reasoning_policy": "entity_list_direct_then_capped_kg"}


def reason_entity_profile(contract, evidence_packet, native_answer=None):
    return {"answer": native_answer or "", "reasoning_policy": "entity_profile_packet_context"}


def reason_path_answer_set(contract, evidence_packet, native_answer=None):
    answer = native_answer or {"answers": [], "paths": []}
    if isinstance(answer, dict):
        return {
            "answer": answer,
            "answers": answer.get("answers", []),
            "paths": answer.get("paths", []),
            "reasoning_policy": "path_answer_set_deterministic",
        }
    return {
        "answer": answer,
        "reasoning_policy": "path_answer_set_deterministic",
    }


def reason_with_contract(contract, evidence_packet=None, native_answer=None):
    task_type = _contract_value(contract, "task_type", "")
    if task_type == "ClosedChoice":
        return reason_closed_choice(contract, evidence_packet, native_answer=native_answer)
    if task_type == "BinaryClaim":
        return reason_binary_claim(contract, evidence_packet, native_answer=native_answer)
    if task_type == "ShortEntity":
        return reason_short_entity(contract, evidence_packet, native_answer=native_answer)
    if task_type == "EntityList":
        return reason_entity_list(contract, evidence_packet, native_answer=native_answer)
    if task_type == "EntityProfile":
        return reason_entity_profile(contract, evidence_packet, native_answer=native_answer)
    if task_type == "PathAnswerSet":
        return reason_path_answer_set(contract, evidence_packet, native_answer=native_answer)
    return {"answer": native_answer or "", "reasoning_policy": "native_fallback"}


def _relation_column(df):
    for col in ("path.1.relation", "path.1.name", "relation", "path.1"):
        if col in df.columns:
            return col
    return None


def _first_existing(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col
    return None


def _filter_df(df, args):
    if getattr(args, "source", None) and "path.0.name" in df:
        df = df[df["path.0.name"].astype(str).str.contains(args.source, case=False, na=False)]
    if getattr(args, "target", None) and "path.2.name" in df:
        df = df[df["path.2.name"].astype(str).str.contains(args.target, case=False, na=False)]
    relation_col = _relation_column(df)
    if getattr(args, "relation", None) and relation_col:
        df = df[df[relation_col].astype(str).str.contains(args.relation, case=False, na=False)]
    return df


def _write_or_print(text, output=None):
    if output:
        Path(output).write_text(text, encoding="utf-8")
        print(f"Wrote {output}")
    else:
        print(text)


def _series_or_default(df, column, default=None):
    if column and column in df.columns:
        return df[column]
    return pd.Series([default] * len(df), index=df.index)


def _has_text(value):
    if value is None:
        return False
    try:
        if pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return bool(text) and text.casefold() not in {"nan", "none", "null", "[]", "{}"}


def _clean_text(value):
    return str(value).strip() if _has_text(value) else ""


def _entity_terms(value):
    text = _clean_text(value)
    if not text:
        return []
    text = re.sub(r"[\[\]\(\){}\"']", " ", text)
    pieces = re.split(r"[;,|/]+|\s+\bor\b\s+", text)
    terms = []
    for piece in pieces:
        term = re.sub(r"\s+", " ", piece).strip()
        if len(term) >= 2 and term.casefold() not in {"nan", "none", "null"}:
            terms.append(term)
    return terms


def _entity_in_text(entity, text):
    haystack = _clean_text(text).casefold()
    if not haystack:
        return False
    for term in _entity_terms(entity):
        needle = term.casefold()
        if len(needle) >= 2 and needle in haystack:
            return True
    return False


WEAK_RELATION_PATTERNS = [
    r"co[-_ ]?mention",
    r"co[-_ ]?occurr",
    r"co[-_ ]?expression",
    r"text[-_ ]?min",
    r"literature[-_ ]?min",
    r"abstract",
]


def _weak_evidence_reason(row):
    """Return a reason string when an edge should be treated as weak co-mention-like evidence."""
    haystack = " ".join(
        _clean_text(row.get(col))
        for col in ["relation", "evidence_text", "condition", "info_source", "edge_support_class"]
    ).casefold()
    for pattern in WEAK_RELATION_PATTERNS:
        if re.search(pattern, haystack):
            return f"matched weak-evidence pattern: {pattern}"
    return ""


def _co_mentioned_edge_review(row):
    """Classify weak co-mentioned/text-mined edges using edge metadata, not relation type alone."""
    weak_reason = _weak_evidence_reason(row)
    if not weak_reason:
        return "structured_relation", 1.0, "Structured KG relation; inspect edge provenance as needed."

    snippet = _clean_text(row.get("evidence_text"))
    condition = _clean_text(row.get("condition"))
    info_source = _clean_text(row.get("info_source"))
    confidence = _to_float(row.get("confidence"), None) if "_to_float" in globals() else None
    edge_text = " ".join(x for x in [snippet, condition] if x)
    source_hit = _entity_in_text(row.get("source"), edge_text)
    target_hit = _entity_in_text(row.get("target"), edge_text)

    detail_score = 0.0
    if snippet:
        detail_score += 0.40
    if source_hit and target_hit:
        detail_score += 0.25
    elif source_hit or target_hit:
        detail_score += 0.10
    if condition:
        detail_score += 0.10
    if info_source:
        detail_score += 0.10
    if confidence is not None:
        detail_score += 0.15 * _bounded(confidence, 0.5) if "_bounded" in globals() else 0.10
    detail_score = max(0.0, min(1.0, detail_score))

    if snippet and source_hit and target_hit:
        return (
            "co_mentioned_edge_text_with_endpoints",
            detail_score,
            "Review original_relation: edge text mentions both endpoints, but co_mentioned still needs semantic validation.",
        )
    if snippet:
        return (
            "co_mentioned_edge_text",
            detail_score,
            "Review original_relation: edge has text, but endpoint/semantic support must be checked.",
        )
    if condition or info_source or confidence is not None:
        return (
            "co_mentioned_metadata_only",
            detail_score,
            "Only edge metadata is available; do not treat as direct biological support.",
        )
    return (
        "co_mentioned_unresolved",
        detail_score,
        "No useful edge text or metadata found; treat as weak background co-occurrence.",
    )


def _edge_is_co_mentioned(row):
    """Detect co-mention/co-occurrence evidence from relation and metadata."""
    haystack = " ".join(
        _clean_text(row.get(col))
        for col in ["relation", "original_relation", "evidence_text", "condition", "info_source", "edge_support_class"]
    ).casefold()
    patterns = [
        "co_mentioned",
        "co-mentioned",
        "co mentioned",
        "cooccur",
        "co-occur",
        "co occurrence",
        "co-occurrence",
        "mentioned together",
    ]
    return any(pattern in haystack for pattern in patterns)


CO_MENTIONED_RELATION_RE = re.compile(r"co[-_ ]?mention", flags=re.I)

RELATION_CUE_PATTERNS = [
    r"\bmarker\b",
    r"\bregulat",
    r"\bexpress",
    r"\bassociat",
    r"\bcaus",
    r"\binhibit",
    r"\btreat",
    r"\bbind",
    r"\binteract",
    r"\brole in\b",
    r"\bessential for\b",
    r"\binvolved in\b",
    r"\brelated to\b",
    r"\bmediat",
    r"\btarget",
    r"\bpathway",
    r"\bsignal",
    r"\bdevelop",
    r"\bdifferentiat",
    r"\bconvert",
    r"\btransit",
    r"\bpharmaco",
    r"\bindicat",
    r"\boverexpress",
    r"\bknockdown",
    r"\bknockout",
]


def _is_co_mentioned_relation(value):
    return bool(CO_MENTIONED_RELATION_RE.search(_clean_text(value)))


def _snippet_shows_relation_cue(snippet):
    text = _clean_text(snippet)
    if not text:
        return False
    for pattern in RELATION_CUE_PATTERNS:
        if re.search(pattern, text, flags=re.I):
            return True
    return False


def _edge_is_weak_co_mentioned(edge):
    if _is_co_mentioned_relation(edge.get("relation")):
        return True
    pseudo = {
        "relation": edge.get("relation"),
        "evidence_text": edge.get("original_relation"),
        "condition": edge.get("condition"),
        "info_source": edge.get("info_source"),
        "edge_support_class": "",
    }
    return bool(_weak_evidence_reason(pseudo))


def _co_mentioned_edge_keep_decision(
    edge,
    *,
    min_detail_score=0.55,
    require_both_endpoints=True,
    require_relation_cue=True,
    mode="strict",
):
    """Return keep, reason, support_class, detail_score for one path edge."""
    if not _edge_is_weak_co_mentioned(edge):
        return True, "Structured KG relation; no co-mentioned filtering applied.", "structured_relation", 1.0

    pseudo = {
        "relation": edge.get("relation"),
        "evidence_text": edge.get("original_relation"),
        "condition": edge.get("condition"),
        "info_source": edge.get("info_source"),
        "confidence": edge.get("confidence"),
        "source": edge.get("source"),
        "target": edge.get("target"),
    }
    support_class, detail_score, note = _co_mentioned_edge_review(pseudo)
    snippet = edge.get("original_relation")

    if mode == "permissive":
        if support_class in {"co_mentioned_unresolved", "co_mentioned_metadata_only"}:
            return False, note, support_class, detail_score
        return True, note or "Co-mentioned edge retained in permissive mode.", support_class, detail_score

    if not _has_text(snippet):
        return False, "Missing path.*.original_relation snippet; exclude co-mentioned row.", support_class, detail_score

    if require_both_endpoints and not (
        _entity_in_text(edge.get("source"), snippet) and _entity_in_text(edge.get("target"), snippet)
    ):
        return (
            False,
            "original_relation does not mention both endpoints; exclude co-mentioned row.",
            support_class,
            detail_score,
        )

    if support_class in {"co_mentioned_unresolved", "co_mentioned_metadata_only"}:
        return False, note, support_class, detail_score

    if detail_score < float(min_detail_score):
        return (
            False,
            f"Co-mentioned detail score {detail_score:.2f} below threshold {float(min_detail_score):.2f}.",
            support_class,
            detail_score,
        )

    if require_relation_cue and not _snippet_shows_relation_cue(snippet):
        return (
            False,
            "original_relation shows co-occurrence only; no surface relationship cue in the sentence.",
            support_class,
            detail_score,
        )

    return True, note or "Co-mentioned edge retained after original_relation review.", support_class, detail_score


CURATED_SOURCE_HINTS = {
    "cellmarker", "panglaodb", "celltaxonomy", "drugbank", "drugcentral", "iuphar", "pharmgkb",
    "clinvar", "omim", "medgen", "orphanet", "hpo", "mondo", "doid", "reactome", "kegg", "goa", "rhea",
    "biogrid", "intact", "corum", "uniprot", "chebi", "hmdb", "foodon", "metacyc",
    "gwas_catalog", "gwas catalog", "gwas_cataog",
}
NETWORK_OR_PREDICTED_SOURCE_HINTS = {
    "string", "stringdb", "stitch", "biosnap", "primekg", "sider", "repodb", "scenic",
}
MAPPING_RELATIONS = {"equals_to", "same_as", "xref", "mapped_to"}
NON_ASSERTIVE_RELATIONS = {"no_relation"}


def _source_quality(row):
    """Coarse source-quality tier from edge-level info_source metadata."""
    text = _clean_text(row.get("info_source")).casefold()
    if not text:
        return "unspecified_source"
    if any(hint in text for hint in CURATED_SOURCE_HINTS):
        return "curated_or_manual_source"
    if any(hint in text for hint in NETWORK_OR_PREDICTED_SOURCE_HINTS):
        return "network_textmined_or_predicted_source"
    return "other_source"


def _evidence_strength(row, support_class):
    relation = _clean_text(row.get("relation")).casefold()
    quality = _source_quality(row)
    confidence_text = _clean_text(row.get("confidence"))
    try:
        confidence = float(confidence_text) if confidence_text else None
        if confidence is not None and confidence > 1.0 and confidence <= 100.0:
            confidence /= 100.0
    except Exception:
        confidence = None

    if relation in NON_ASSERTIVE_RELATIONS:
        return "not_biological_support"
    if relation in MAPPING_RELATIONS:
        return "mapping_only"
    if _edge_is_co_mentioned(row) or str(support_class).startswith("co_mentioned"):
        return "weak_co_mentioned_review_required"
    if confidence is not None and confidence < 0.35:
        return "weak_low_confidence"
    if quality == "curated_or_manual_source" and relation not in {"is_a"}:
        return "strong_curated"
    if quality == "network_textmined_or_predicted_source":
        return "medium_network_or_predicted"
    return "medium_structured"


def _ensure_edge_support_columns(evidence):
    evidence = evidence.copy()
    classes = []
    detail_scores = []
    notes = []
    source_qualities = []
    strengths = []
    review_required = []
    evidence_groups = []
    answer_roles = []
    boundary_notes = []
    for _, row in evidence.iterrows():
        support_class, detail_score, note = _co_mentioned_edge_review(row)
        # Do not let non-relation metadata hide co-mention/co-occurrence evidence.
        if support_class == "structured_relation" and _edge_is_co_mentioned(row):
            support_class = "co_mentioned_metadata_detected"
            detail_score = min(float(detail_score or 0.0), 0.35)
            note = "Co-mention/co-occurrence language appears in edge metadata; treat as weak candidate evidence and review original_relation/info_source/condition/confidence."
        quality = _source_quality(row)
        strength = _evidence_strength(row, support_class)
        classes.append(support_class)
        detail_scores.append(round(float(detail_score or 0.0), 6))
        notes.append(note)
        source_qualities.append(quality)
        strengths.append(strength)
        review_required.append(strength.startswith("weak") or strength in {"mapping_only", "not_biological_support"})
        if classify_evidence_row is not None:
            boundary = classify_evidence_row(row.to_dict())
            evidence_groups.append(boundary.get("evidence_group", "background"))
            answer_roles.append(boundary.get("answer_role", "background_context"))
            boundary_notes.append(boundary.get("reason", ""))
        else:
            evidence_groups.append("background")
            answer_roles.append("background_context")
            boundary_notes.append("")
    evidence["edge_support_class"] = classes
    evidence["edge_detail_score"] = detail_scores
    evidence["edge_review_note"] = notes
    evidence["source_quality"] = source_qualities
    evidence["evidence_strength"] = strengths
    evidence["manual_review_required"] = review_required
    evidence["evidence_group"] = evidence_groups
    evidence["answer_role"] = answer_roles
    evidence["evidence_boundary_note"] = boundary_notes
    return evidence

def _evidence_id(row):
    key = "|".join(
        str(row.get(col, ""))
        for col in [
            "path_id",
            "edge_position",
            "source",
            "relation",
            "target",
            "info_source",
            "condition",
            "original_relation",
            "evidence_text",
            "provenance_file",
        ]
    )
    return hashlib.sha1(key.encode("utf-8")).hexdigest()[:16]

def _path_indices(df):
    node_indices = sorted(
        {int(match.group(1)) for col in df.columns for match in [re.match(r"^path\.(\d+)\.name$", col)] if match}
    )
    edge_indices = sorted(
        {int(match.group(1)) for col in df.columns for match in [re.match(r"^path\.(\d+)\.relation$", col)] if match}
    )
    return node_indices, edge_indices


def _cell(row, column, default=None):
    if column in row.index:
        return row.get(column, default)
    return default


def _normalize_path_rows(df, provenance_file=None):
    """Normalize every edge in a path.N KG result into an evidence unit."""
    _, edge_indices = _path_indices(df)
    if not edge_indices:
        return None

    rows = []
    path_length = max(1, len(edge_indices))
    for row_idx, row in df.iterrows():
        path_source = _cell(row, "path.0.name", _cell(row, "path.0.id", ""))
        path_target_index = 2 * path_length
        path_target = _cell(row, f"path.{path_target_index}.name", _cell(row, f"path.{path_target_index}.id", ""))
        for hop_number, edge_idx in enumerate(edge_indices, start=1):
            left_idx = edge_idx - 1
            right_idx = edge_idx + 1
            if f"path.{left_idx}.name" not in df.columns or f"path.{right_idx}.name" not in df.columns:
                continue
            out = {
                "source": _cell(row, f"path.{left_idx}.name", _cell(row, f"path.{left_idx}.id", "")),
                "source_type": _cell(row, f"path.{left_idx}.type", ""),
                "source_id": _cell(row, f"path.{left_idx}.id", ""),
                "relation": _cell(row, f"path.{edge_idx}.relation", _cell(row, f"path.{edge_idx}.name", "")),
                "target": _cell(row, f"path.{right_idx}.name", _cell(row, f"path.{right_idx}.id", "")),
                "target_type": _cell(row, f"path.{right_idx}.type", ""),
                "target_id": _cell(row, f"path.{right_idx}.id", ""),
                "confidence": _cell(row, f"path.{edge_idx}.relation_confidence", _cell(row, "relation_confidence", "")),
                "info_source": _cell(row, f"path.{edge_idx}.info_source", _cell(row, "info_source", "")),
                "info_source_count": _cell(row, f"path.{edge_idx}.info_source_length", _cell(row, "info_source_count", "")),
                "condition": _cell(row, f"path.{edge_idx}.condition", _cell(row, "condition", "")),
                "evidence_text": _cell(row, f"path.{edge_idx}.original_relation", _cell(row, "original_relation", _cell(row, "evidence_text", ""))),
                "provenance_file": provenance_file or "",
                "evidence_kind": "kg_path_edge" if path_length > 1 else "kg_edge",
                "direction": "source_to_target",
                "path_row_id": row_idx,
                "hop_index": hop_number,
                "path_length": path_length,
                "path_source": path_source,
                "path_target": path_target,
            }
            rows.append(out)

    if not rows:
        return None
    evidence = pd.DataFrame(rows)
    evidence["evidence_id"] = evidence.apply(_evidence_id, axis=1)
    columns = [
        "evidence_id",
        "source",
        "source_type",
        "source_id",
        "relation",
        "target",
        "target_type",
        "target_id",
        "confidence",
        "info_source",
        "info_source_count",
        "condition",
        "evidence_text",
        "provenance_file",
        "evidence_kind",
        "direction",
        "path_row_id",
        "hop_index",
        "path_length",
        "path_source",
        "path_target",
    ]
    return _ensure_edge_support_columns(evidence[columns])


def _row_first(row, candidates, default=""):
    """Return the first non-empty value from a row, falling back to the first present value."""
    present = []
    for col in candidates:
        if col in row.index:
            value = row.get(col)
            present.append(value)
            if _has_text(value):
                return value
    return present[0] if present else default


def _path_edge_indices(df):
    indices = []
    for col in df.columns:
        match = re.match(r"^path\.(\d+)\.(?:relation|name)$", str(col))
        if not match:
            continue
        idx = int(match.group(1))
        if idx % 2 == 1 and (f"path.{idx - 1}.name" in df.columns or f"path.{idx - 1}.id" in df.columns) and (f"path.{idx + 1}.name" in df.columns or f"path.{idx + 1}.id" in df.columns):
            indices.append(idx)
    return sorted(set(indices))


def _path_fingerprint(row, edge_indices):
    bits = []
    for edge_idx in edge_indices:
        for col in [
            f"path.{edge_idx - 1}.id",
            f"path.{edge_idx - 1}.name",
            f"path.{edge_idx}.relation",
            f"path.{edge_idx}.info_source",
            f"path.{edge_idx + 1}.id",
            f"path.{edge_idx + 1}.name",
        ]:
            bits.append(str(row.get(col, "")))
    return hashlib.sha1("|".join(bits).encode("utf-8")).hexdigest()[:16]


def _active_path_edges(row, edge_indices):
    active = []
    for edge_idx in edge_indices:
        source_idx = edge_idx - 1
        target_idx = edge_idx + 1
        relation = _row_first(row, [f"path.{edge_idx}.relation", f"path.{edge_idx}.name", "relation"])
        source_value = _row_first(row, [f"path.{source_idx}.name", f"path.{source_idx}.id", "source"])
        target_value = _row_first(row, [f"path.{target_idx}.name", f"path.{target_idx}.id", "target"])
        if _has_text(source_value) and _has_text(relation) and _has_text(target_value):
            active.append(edge_idx)
    return active


def _row_edge_snapshot(row, edge_idx):
    source_idx = edge_idx - 1
    target_idx = edge_idx + 1
    return {
        "edge_index": edge_idx,
        "source": _row_first(row, [f"path.{source_idx}.name", f"path.{source_idx}.id", "source"]),
        "target": _row_first(row, [f"path.{target_idx}.name", f"path.{target_idx}.id", "target"]),
        "relation": _row_first(row, [f"path.{edge_idx}.relation", f"path.{edge_idx}.name", "relation"]),
        "original_relation": _row_first(row, [f"path.{edge_idx}.original_relation", "original_relation"]),
        "condition": _row_first(row, [f"path.{edge_idx}.condition", "condition"]),
        "info_source": _row_first(row, [f"path.{edge_idx}.info_source", "info_source"]),
        "confidence": _row_first(row, [f"path.{edge_idx}.relation_confidence", "confidence", "relation_confidence"]),
    }


def _query_row_co_mentioned_decision(row, edge_indices, **options):
    """Decide whether a brick qr row should be kept after co-mentioned review."""
    active_edges = _active_path_edges(row, edge_indices) if edge_indices else []
    if not active_edges:
        active_edges = [1]

    co_edges = [edge_idx for edge_idx in active_edges if _edge_is_weak_co_mentioned(_row_edge_snapshot(row, edge_idx))]
    if not co_edges:
        return True, "No co-mentioned edges on path.", "structured_relation", 1.0

    last_support = "co_mentioned_edge_text_with_endpoints"
    last_score = 0.0
    kept_notes = []
    for edge_idx in co_edges:
        keep, reason, support_class, detail_score = _co_mentioned_edge_keep_decision(
            _row_edge_snapshot(row, edge_idx),
            **options,
        )
        last_support = support_class
        last_score = detail_score
        if not keep:
            return False, f"path.{edge_idx}: {reason}", support_class, detail_score
        kept_notes.append(f"path.{edge_idx}: {reason}")

    return True, " | ".join(kept_notes[:3]), last_support, last_score


def filter_co_mentioned_df(
    df,
    *,
    min_detail_score=0.55,
    require_both_endpoints=True,
    require_relation_cue=True,
    mode="strict",
):
    """Filter brick qr rows, dropping co-mentioned edges whose snippets fail review."""
    edge_indices = _path_edge_indices(df)
    keep_flags = []
    review_notes = []
    support_classes = []
    detail_scores = []

    for _, row in df.iterrows():
        keep, reason, support_class, detail_score = _query_row_co_mentioned_decision(
            row,
            edge_indices,
            min_detail_score=min_detail_score,
            require_both_endpoints=require_both_endpoints,
            require_relation_cue=require_relation_cue,
            mode=mode,
        )
        keep_flags.append(bool(keep))
        review_notes.append(reason)
        support_classes.append(support_class)
        detail_scores.append(detail_score)

    audit = df.copy()
    audit["co_mentioned_keep"] = keep_flags
    audit["co_mentioned_review_note"] = review_notes
    audit["co_mentioned_support_class"] = support_classes
    audit["co_mentioned_detail_score"] = detail_scores

    kept = audit[audit["co_mentioned_keep"]].copy()
    rejected = audit[~audit["co_mentioned_keep"]].copy()
    return kept, rejected, audit


def cmd_filter_co_mentioned(args):
    df = pd.read_csv(args.input)
    df = _filter_df(df, args)

    kept, rejected, audit = filter_co_mentioned_df(
        df,
        min_detail_score=args.min_detail_score,
        require_both_endpoints=not args.allow_missing_endpoint,
        require_relation_cue=not args.keep_bare_co_occurrence,
        mode=args.mode,
    )

    audit_cols = [
        "co_mentioned_keep",
        "co_mentioned_review_note",
        "co_mentioned_support_class",
        "co_mentioned_detail_score",
    ]

    if args.audit_output:
        audit.to_csv(args.audit_output, index=False)
        print(f"Wrote audit table {args.audit_output}")

    if args.rejected_output:
        rejected.to_csv(args.rejected_output, index=False)
        print(f"Wrote rejected rows {args.rejected_output} ({len(rejected)} rows)")

    if not args.output:
        raise ValueError("-o/--output is required")

    output_df = kept if args.include_audit_columns else kept.drop(columns=audit_cols, errors="ignore")
    output_df.to_csv(args.output, index=False)

    dropped = len(df) - len(output_df)
    print(
        f"Wrote {args.output} ({len(output_df)} kept / {dropped} dropped co-mentioned rows from {len(df)} input rows)"
    )


def _generic_evidence_record(row, provenance_file=None):
    relation_col = _first_existing(pd.DataFrame([row]), ["relation", "path.1.relation", "path.1.name", "path.1"])
    return {
        "source": _row_first(row, ["source", "path.0.name", "path.0.id"]),
        "source_type": _row_first(row, ["source_type", "path.0.type"]),
        "source_id": _row_first(row, ["source_id", "path.0.id"]),
        "relation": _row_first(row, [relation_col, "relation", "path.1.relation", "path.1.name", "path.1"]),
        "target": _row_first(row, ["target", "path.2.name", "path.2.id"]),
        "target_type": _row_first(row, ["target_type", "path.2.type"]),
        "target_id": _row_first(row, ["target_id", "path.2.id"]),
        "confidence": _row_first(row, ["confidence", "path.1.relation_confidence", "relation_confidence"]),
        "info_source": _row_first(row, ["info_source", "path.1.info_source"]),
        "info_source_count": _row_first(row, ["info_source_count", "path.1.info_source_length", "info_source_length"]),
        "condition": _row_first(row, ["condition", "path.1.condition"]),
        "original_relation": _row_first(row, ["original_relation", "path.1.original_relation"]),
        "evidence_text": _row_first(row, ["evidence_text", "path.1.original_relation", "original_relation", "snippet", "evidence"]),
        "provenance_file": provenance_file or "",
        "evidence_kind": "kg_edge",
        "direction": "source_to_target",
        "path_id": "",
        "path_hops": 1,
        "edge_index": 1,
        "edge_position": 1,
    }


def normalize_evidence_df(df, provenance_file=None):
    """Convert brick query-like rows into canonical evidence units.

    Multi-hop ``brick qr --hops 2/3`` rows are expanded into one evidence unit
    per edge.  This prevents a path such as A-r1-B-r2-C from being compressed
    into a misleading one-hop A->B relation and keeps every edge's
    original_relation/info_source/condition/confidence available for review.
    """
    edge_indices = _path_edge_indices(df)
    records = []

    for _, row in df.iterrows():
        if not edge_indices:
            records.append(_generic_evidence_record(row, provenance_file=provenance_file))
            continue

        active_edges = _active_path_edges(row, edge_indices)
        if not active_edges:
            records.append(_generic_evidence_record(row, provenance_file=provenance_file))
            continue

        path_id = _path_fingerprint(row, active_edges)
        path_hops = len(active_edges)
        for position, edge_idx in enumerate(active_edges, start=1):
            source_idx = edge_idx - 1
            target_idx = edge_idx + 1
            relation = _row_first(row, [f"path.{edge_idx}.relation", f"path.{edge_idx}.name", "relation"])
            source_value = _row_first(row, [f"path.{source_idx}.name", f"path.{source_idx}.id", "source"])
            target_value = _row_first(row, [f"path.{target_idx}.name", f"path.{target_idx}.id", "target"])
            if not (_has_text(source_value) and _has_text(relation) and _has_text(target_value)):
                continue
            rec = {
                "source": source_value,
                "source_type": _row_first(row, [f"path.{source_idx}.type", "source_type"]),
                "source_id": _row_first(row, [f"path.{source_idx}.id", "source_id"]),
                "relation": relation,
                "target": target_value,
                "target_type": _row_first(row, [f"path.{target_idx}.type", "target_type"]),
                "target_id": _row_first(row, [f"path.{target_idx}.id", "target_id"]),
                "confidence": _row_first(row, [f"path.{edge_idx}.relation_confidence", "confidence", "relation_confidence"]),
                "info_source": _row_first(row, [f"path.{edge_idx}.info_source", "info_source"]),
                "info_source_count": _row_first(row, [f"path.{edge_idx}.info_source_length", "info_source_count", "info_source_length"]),
                "condition": _row_first(row, [f"path.{edge_idx}.condition", "condition"]),
                "original_relation": _row_first(row, [f"path.{edge_idx}.original_relation", "original_relation"]),
                "evidence_text": _row_first(row, [f"path.{edge_idx}.original_relation", "evidence_text", "original_relation", "snippet", "evidence"]),
                "provenance_file": provenance_file or "",
                "evidence_kind": "kg_path_edge" if path_hops > 1 else "kg_edge",
                "direction": "source_to_target",
                "path_id": path_id,
                "path_hops": path_hops,
                "edge_index": edge_idx,
                "edge_position": position,
            }
            records.append(rec)

    evidence = pd.DataFrame.from_records(records)
    if evidence.empty:
        evidence = pd.DataFrame(columns=[
            "source", "source_type", "source_id", "relation", "target", "target_type", "target_id",
            "confidence", "info_source", "info_source_count", "condition", "original_relation",
            "evidence_text", "provenance_file", "evidence_kind", "direction", "path_id", "path_hops",
            "edge_index", "edge_position",
        ])
    evidence["evidence_id"] = evidence.apply(_evidence_id, axis=1)
    columns = [
        "evidence_id",
        "source",
        "source_type",
        "source_id",
        "relation",
        "target",
        "target_type",
        "target_id",
        "confidence",
        "info_source",
        "info_source_count",
        "condition",
        "original_relation",
        "evidence_text",
        "provenance_file",
        "evidence_kind",
        "direction",
        "path_id",
        "path_hops",
        "edge_index",
        "edge_position",
    ]
    columns = [col for col in columns if col in evidence.columns]
    return _ensure_edge_support_columns(evidence[columns])

def _load_evidence_units(path):
    df = pd.read_csv(path)
    if {"evidence_id", "source", "relation", "target"}.issubset(df.columns):
        return df
    if {"source", "relation", "target"}.issubset(df.columns):
        evidence = df.copy()
        if "evidence_id" not in evidence.columns:
            evidence["evidence_id"] = evidence.apply(_evidence_id, axis=1)
        return _ensure_edge_support_columns(evidence)
    return normalize_evidence_df(df, provenance_file=path)


def cmd_normalize(args):
    df = pd.read_csv(args.input)
    df = _filter_df(df, args)
    if df.empty:
        raise ValueError("No evidence rows after filtering.")
    evidence = normalize_evidence_df(df, provenance_file=args.provenance_file or args.input)

    if args.min_confidence is not None and "confidence" in evidence:
        confidence = pd.to_numeric(evidence["confidence"], errors="coerce")
        evidence = evidence[confidence.isna() | (confidence >= args.min_confidence)]

    if args.deduplicate:
        dedup_cols = [c for c in ["source", "relation", "target", "info_source", "condition", "evidence_text"] if c in evidence.columns]
        evidence = evidence.drop_duplicates(subset=dedup_cols)

    if args.format == "jsonl":
        text = "\n".join(json for json in evidence.to_json(orient="records", lines=True, force_ascii=False).splitlines())
        if args.output:
            Path(args.output).write_text(text + "\n", encoding="utf-8")
            print(f"Wrote {args.output}")
        else:
            print(text)
        return

    if args.output:
        evidence.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(evidence.to_csv(index=False))


def cmd_packet(args):
    if build_evidence_packet is None:
        raise RuntimeError("evidence_packet module is unavailable")
    evidence = _load_evidence_units(args.input)
    packet = build_evidence_packet(
        evidence,
        expected_answer_type=args.expected_answer_type,
        question_scope_terms=args.scope_term,
        question_scope_match=args.scope_match,
        candidate_terms=args.candidate_term,
        candidate_side=args.candidate_side,
        weak_evidence_policy=args.weak_evidence_policy,
        max_candidate_rows=args.max_candidate_rows,
        schema_file=args.kg_schema,
    )
    if args.format == "prompt":
        text = render_prompt_evidence(packet, max_rows=args.max_rows)
    elif args.format == "csv":
        text = pd.DataFrame(packet_to_rows(packet)).to_csv(index=False)
    else:
        text = json.dumps(packet, ensure_ascii=False, indent=2)
    _write_or_print(text, args.output)


def _target_rows(df, target):
    return df[df["target"].astype(str).str.lower() == str(target).lower()]


def _top_targets_from_evidence(evidence, top_n):
    counts = (
        evidence.groupby("target", dropna=False)
        .agg(match_count=("source", lambda x: len(set(x.dropna().astype(str)))), evidence_count=("evidence_id", "count"))
        .sort_values(["match_count", "evidence_count"], ascending=False)
    )
    return [str(x) for x in counts.head(top_n).index if pd.notna(x)]


def _top_targets_from_hypotheses(path, top_n):
    df = pd.read_csv(path)
    rank_col = _first_existing(df, ["path.2.rank_voting", "path.2.enrich_fdr", "path.2.enrich_pvalue"])
    target_col = _first_existing(df, ["path.2.name", "target"])
    if not target_col:
        return [], df
    if rank_col:
        df = df.sort_values(rank_col, ascending=True)
    return [str(x) for x in df[target_col].dropna().head(top_n)], df


def _value_counts_text(series, max_items=8):
    values = [str(x) for x in series.dropna().astype(str) if _has_text(x)]
    counts = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:max_items]
    return ", ".join(f"{name}={count}" for name, count in ordered)


def _append_global_strength_summary(lines, evidence):
    strength_counts = evidence.get("evidence_strength", pd.Series(dtype=str)).dropna().astype(str).value_counts().to_dict()
    support_counts = evidence.get("edge_support_class", pd.Series(dtype=str)).dropna().astype(str).value_counts().to_dict()
    review_count = int(evidence.get("manual_review_required", pd.Series(False, index=evidence.index)).astype(bool).sum())
    co_mentioned_count = int(evidence.apply(_edge_is_co_mentioned, axis=1).sum())
    if not strength_counts and not support_counts and not review_count and not co_mentioned_count:
        return
    lines.append("## Evidence strength summary")
    for name, count in strength_counts.items():
        lines.append(f"- {name}: {count}")
    if support_counts:
        lines.append("- Edge support classes: " + "; ".join(f"{k}={v}" for k, v in support_counts.items()))
    if review_count:
        lines.append(f"- Manual-review rows: {review_count}")
    if co_mentioned_count:
        lines.append(f"- Co-mentioned/text-mined rows: {co_mentioned_count}")
    lines.append("")


def cmd_explain(args):
    evidence = _load_evidence_units(args.evidence)
    if "evidence_score" not in evidence.columns or "evidence_strength" not in evidence.columns:
        evidence = score_evidence_df(evidence)
    if args.target:
        targets = [args.target]
        hypotheses = None
    elif args.hypotheses:
        targets, hypotheses = _top_targets_from_hypotheses(args.hypotheses, args.top_n)
    else:
        targets = _top_targets_from_evidence(evidence, args.top_n)
        hypotheses = None

    if not targets:
        raise ValueError("No targets found to explain.")

    lines = ["# brick Interpretability Report", ""]
    lines.append(f"- Evidence file: {args.evidence}")
    if args.hypotheses:
        lines.append(f"- Hypothesis file: {args.hypotheses}")
    lines.append("- Evidence contract: co-mentioned/co-occurrence rows are weak candidate evidence and require review of original_relation, info_source, condition, and relation_confidence.")
    lines.append("")
    _append_global_strength_summary(lines, evidence)

    for target in targets:
        group = _target_rows(evidence, target)
        if group.empty:
            continue
        sources = sorted(set(group.get("source", pd.Series(dtype=str)).dropna().astype(str)))
        relations = sorted(set(group.get("relation", pd.Series(dtype=str)).dropna().astype(str)))
        info_count = pd.to_numeric(group.get("info_source_count", pd.Series(dtype=float)), errors="coerce").sum()
        confidence = pd.to_numeric(group.get("confidence", pd.Series(dtype=float)), errors="coerce").mean()
        evidence_score = pd.to_numeric(group.get("evidence_score", pd.Series(dtype=float)), errors="coerce").mean()
        strength_counts = group.get("evidence_strength", pd.Series(dtype=str)).dropna().astype(str).value_counts().to_dict()
        review_count = int(group.get("manual_review_required", pd.Series(False, index=group.index)).astype(bool).sum())
        co_mentioned_count = int(group.apply(_edge_is_co_mentioned, axis=1).sum())

        lines.append(f"## {target}")
        lines.append(f"- Matched source entities: {len(sources)}")
        lines.append(f"- Evidence rows: {len(group)}")
        if pd.notna(info_count):
            lines.append(f"- Information source count: {info_count:g}")
        if pd.notna(confidence):
            lines.append(f"- Mean confidence: {confidence:.3f}")
        if pd.notna(evidence_score):
            lines.append(f"- Mean evidence score: {evidence_score:.3f}")
        if review_count:
            lines.append(f"- Manual-review rows: {review_count}")
        if co_mentioned_count:
            lines.append(f"- Co-mentioned/co-occurrence rows: {co_mentioned_count} (weak candidates, not direct support)")
        if co_mentioned_count and co_mentioned_count == len(group):
            lines.append("- Review required: this claim is supported only by co-mentioned/text-mined evidence. Do not state it as a confirmed mechanistic relation.")
        if strength_counts:
            lines.append("- Evidence strength: " + "; ".join(f"{k}={v}" for k, v in strength_counts.items()))
        if relations:
            lines.append(f"- Relations: {', '.join(relations)}")
        if sources:
            lines.append(f"- Source examples: {', '.join(sources[:args.max_sources])}")

        if args.hypotheses and hypotheses is not None:
            target_col = _first_existing(hypotheses, ["path.2.name", "target"])
            rank_cols = [
                col for col in [
                    "path.2.rank_voting",
                    "path.2.match_count",
                    "path.2.info_source_count",
                    "path.2.enrich_pvalue",
                    "path.2.enrich_fdr",
                    "net_support_score",
                    "support_score",
                ]
                if col in hypotheses.columns
            ]
            if target_col:
                hrow = hypotheses[hypotheses[target_col].astype(str).str.lower() == str(target).lower()]
                if not hrow.empty and rank_cols:
                    metrics = [f"{col}={hrow.iloc[0].get(col)}" for col in rank_cols]
                    lines.append(f"- Ranking/enrichment metrics: {'; '.join(metrics)}")

        edge_rows = group.copy()
        if "evidence_score" in edge_rows:
            edge_rows = edge_rows.sort_values("evidence_score", ascending=False)
        edge_rows = edge_rows.head(args.max_snippets)
        if len(edge_rows):
            lines.append("- Evidence rows to inspect:")
            for _, edge in edge_rows.iterrows():
                details = []
                for label, col in [
                    ("relation", "relation"),
                    ("strength", "evidence_strength"),
                    ("score", "evidence_score"),
                    ("support", "edge_support_class"),
                    ("source", "info_source"),
                    ("condition", "condition"),
                    ("confidence", "confidence"),
                    ("path", "path_id"),
                    ("edge_pos", "edge_position"),
                ]:
                    value = edge.get(col)
                    if _has_text(value):
                        details.append(f"{label}={value}")
                lines.append(f"  - [{' ; '.join(details)}] {str(edge.get('evidence_text', ''))[:300]}")
                note = edge.get("edge_review_note")
                if _has_text(note) and (edge.get("manual_review_required") or _edge_is_co_mentioned(edge)):
                    lines.append(f"    Review note: {note}")
        lines.append("")

    _write_or_print("\n".join(lines), args.output)

def _normalize_label(value):
    text = "" if pd.isna(value) else str(value)
    text = text.strip().casefold()
    text = re.sub(r"\s*\([^)]*\)", "", text)
    text = re.sub(r"[^a-z0-9]+", "", text)
    return text


def _split_aliases(value):
    if pd.isna(value):
        return []
    return [x.strip() for x in re.split(r"[|;,]", str(value)) if x.strip()]


def _load_lexicon(path):
    if not path:
        return [], {}
    df = pd.read_csv(path)
    name_col = _first_existing(df, ["canonical_name", "name", "label", "entity_name"])
    id_col = _first_existing(df, ["canonical_id", "entity_id", "id"])
    type_col = _first_existing(df, ["entity_type", "type"])
    alias_col = _first_existing(df, ["aliases", "alias", "synonyms", "synonym"])
    if not name_col:
        raise ValueError("Lexicon must contain a name/canonical_name/label/entity_name column.")

    records = []
    alias_map = {}
    for _, row in df.iterrows():
        name = str(row.get(name_col, "")).strip()
        if not name:
            continue
        record = {
            "canonical_name": name,
            "canonical_id": str(row.get(id_col, name)).strip() if id_col else name,
            "entity_type": str(row.get(type_col, "")).strip() if type_col else "",
            "aliases": [],
        }
        records.append(record)
        aliases = [name]
        if alias_col:
            aliases.extend(_split_aliases(row.get(alias_col)))
        record["aliases"] = sorted(set(x for x in aliases if x))
        for alias in aliases:
            key = _normalize_label(alias)
            if key:
                alias_map.setdefault(key, []).append(record)
    return records, alias_map


CURATED_QUERY_EXPANSIONS = {
    "amd": [("age-related macular degeneration", "Disease"), ("age related macular degeneration", "Disease")],
    "cms": [("consensus molecular subtype", ""), ("consensus molecular subtype classification", "")],
    "crc": [("colorectal cancer", "Disease"), ("colon cancer", "Disease"), ("rectal cancer", "Disease")],
    "dmd": [("Duchenne muscular dystrophy", "Disease"), ("Duchenne muscular dystrophy disease", "Disease")],
    "gbm": [("GBM", "Disease"), ("glioblastoma", "Disease"), ("glioblastoma multiforme", "Disease")],
    "mes": [("MES", "Phenotype"), ("mesenchymal", "Phenotype"), ("mesenchymal-like", "Phenotype"), ("MES-like", "Phenotype")],
    "pde5": [("PDE5", "Protein"), ("phosphodiesterase 5", "Protein"), ("phosphodiesterase type 5", "Protein"), ("PDE5A", "Protein"), ("phosphodiesterase 5A", "Protein")],
    "uchl1": [("UCH-L1", "Gene"), ("ubiquitin C-terminal hydrolase L1", "Gene"), ("ubiquitin carboxyl-terminal hydrolase L1", "Gene"), ("PGP9.5", "Gene")],
    "rfc1": [("replication factor C subunit 1", "Gene"), ("RFC1 repeat expansion", ""), ("RFC1-CANVAS", "Disease")],
    "zikavirus": [("Zika virus", "Disease"), ("ZIKV", "Disease")],
}


PROFILE_QUERY_EXPANSIONS = [
    {
        "pattern": r"\bpde\s*-?\s*5\b|phosphodiesterase\s*5",
        "terms": [("sildenafil", "Chemical"), ("vardenafil", "Chemical"), ("avanafil", "Chemical"), ("mirodenafil", "Chemical"), ("udenafil", "Chemical"), ("tadalafil", "Chemical")],
        "kind": "class_member",
        "reason": "PDE5 inhibitor class members; downstream evidence should separate short- and long-half-life drugs.",
        "weight": 0.88,
    },
    {
        "pattern": r"(zika.*vector|vector.*zika)",
        "terms": [("Zika virus", "Disease"), ("ZIKV", "Disease"), ("mosquito", "Species"), ("Aedes", "Species"), ("Aedes aegypti", "Species"), ("Aedes albopictus", "Species")],
        "kind": "vector_candidate",
        "reason": "Vector questions need organism/species expansions in addition to the pathogen.",
        "weight": 0.86,
    },
    {
        "pattern": r"(glioblastoma|gbm).*(mesenchymal|mes-like|mes like|mes\b).*(transcription|tf)|"
                   r"(transcription|tf).*(mesenchymal|mes-like|mes like|mes\b).*(glioblastoma|gbm)",
        "terms": [("CEBPB", "Gene"), ("CEBPG", "Gene"), ("STAT3", "Gene"), ("TFAP2A", "Gene"), ("FOSL1", "Gene"), ("JUN", "Gene"), ("FOXD1", "Gene"), ("ZEB1", "Gene"), ("SNAI1", "Gene"), ("SNAI2", "Gene"), ("TWIST1", "Gene"), ("TEAD2", "Gene"), ("TEAD4", "Gene"), ("RUNX1", "Gene"), ("NFATC2", "Gene")],
        "kind": "candidate_gene",
        "reason": "MES-like GBM TF questions benefit from candidate TF expansion before KG validation.",
        "weight": 0.84,
    },
    {
        "pattern": r"(cms|consensus molecular subtype).*(colorectal|crc)|"
                   r"(colorectal|crc).*(cms|consensus molecular subtype)",
        "terms": [("consensus molecular subtype", ""), ("CMS1", ""), ("CMS2", ""), ("CMS3", ""), ("CMS4", ""), ("mixed CMS", ""), ("unclassified CMS", ""), ("colorectal cancer", "Disease")],
        "kind": "subtype_term",
        "reason": "CMS questions need subtype labels and unclassified/mixed terminology.",
        "weight": 0.82,
    },
    {
        "pattern": r"rfc1.*(repeat|aaagg|aaggg|canvas)|(?:aaagg|aaggg).*rfc1",
        "terms": [("RFC1", "Gene"), ("AAAGG", ""), ("AAGGG", ""), ("ACAGG", ""), ("CANVAS", "Disease"), ("cerebellar ataxia neuropathy vestibular areflexia syndrome", "Disease")],
        "kind": "motif_context",
        "reason": "Repeat motif questions need exact motif plus comparator motif and disease context.",
        "weight": 0.82,
    },
]


QUERY_INTENT_HINTS = [
    (r"\b(vector|transmit|transmission)\b", "target_type_hint", "Species", "Vector questions usually need Species or organism targets."),
    (r"\b(treatment|therapy|drug|inhibitor|chemotherapy|bevacizumab|folfoxiri)\b", "target_type_hint", "Chemical", "Treatment questions usually need Chemical targets."),
    (r"\b(transcription factor|tf|regulat|driver)\b", "target_type_hint", "Gene,Protein", "TF/driver questions usually need Gene or Protein targets."),
    (r"\b(variant|mutation|loss-of-function|lof|caused by)\b", "target_type_hint", "Disease,Phenotype", "Variant-causality questions usually need Disease or Phenotype targets."),
    (r"\b(marker|biomarker)\b", "relation_hint", "marker_of", "Marker questions should try KG relation marker_of."),
    (r"\b(treatment|therapy|drug|inhibitor)\b", "relation_hint", "pharmaco_effect_on,interact_with", "Treatment questions should try pharmaco_effect_on and interact_with."),
]


def _add_expansion(rows, term, kind, source, weight, reason, entity_type=None):
    term = "" if pd.isna(term) else str(term).strip()
    if not term:
        return
    rows.append({
        "term": term,
        "kind": kind,
        "source": source,
        "weight": float(weight),
        "reason": reason,
        "entity_type": entity_type or "",
    })


def _dedupe_expansion_rows(rows, max_terms):
    best = {}
    for row in rows:
        key = (_normalize_label(row.get("term")), str(row.get("kind", "")))
        if not key[0]:
            continue
        if key not in best or row.get("weight", 0) > best[key].get("weight", 0):
            best[key] = row
    ordered = sorted(best.values(), key=lambda x: (-float(x.get("weight", 0)), str(x.get("term", "")).lower()))
    if max_terms:
        query_rows = [r for r in ordered if r.get("kind") not in {"target_type_hint", "relation_hint"}]
        hint_rows = [r for r in ordered if r.get("kind") in {"target_type_hint", "relation_hint"}]
        ordered = query_rows[:max_terms] + hint_rows
    return pd.DataFrame(ordered, columns=["term", "kind", "source", "weight", "reason", "entity_type"])


def _curated_type_for_term(term):
    key = _normalize_label(term)
    if not key:
        return ""
    for values in CURATED_QUERY_EXPANSIONS.values():
        for item in values:
            value, entity_type = item if isinstance(item, (list, tuple)) else (item, "")
            if _normalize_label(value) == key and entity_type:
                return entity_type
    for item in PROFILE_QUERY_EXPANSIONS:
        for term_item in item["terms"]:
            value, entity_type = term_item if isinstance(term_item, (list, tuple)) else (term_item, "")
            if _normalize_label(value) == key and entity_type:
                return entity_type
    return ""


def query_terms_from_expansion(expanded, entity_type=None, max_terms=None):
    """Return terms suitable for a brick qr call after type-aware filtering."""
    if expanded is None or expanded.empty:
        return []
    blocked = {"target_type_hint", "relation_hint"}
    entity_type_norm = (entity_type or "").casefold()
    terms = []
    seen = set()
    for _, row in expanded.iterrows():
        kind = str(row.get("kind", ""))
        if kind in blocked:
            continue
        row_type = str(row.get("entity_type", "") or "").casefold()
        if entity_type_norm and row_type and row_type != entity_type_norm:
            continue
        term = str(row.get("term", "")).strip()
        key = term.casefold()
        if term and key not in seen:
            seen.add(key)
            terms.append(term)
        if max_terms and len(terms) >= max_terms:
            break
    return terms


QUESTION_LEAD_RE = re.compile(
    r"^(what|which|who|where|when|why|how|does|do|did|is|are|was|were|can|could|please|list|name|identify)\b",
    flags=re.I,
)

PHRASE_LEADING_STOPWORDS = {
    "a", "an", "and", "any", "are", "as", "be", "been", "between", "by", "does", "for", "from", "has",
    "have", "in", "into", "is", "it", "main", "most", "of", "on", "or", "the", "there", "to", "used",
    "what", "which", "with", "prevalence", "cause", "causes", "symptoms", "indications", "adverse",
    "effects", "effect", "role", "vector",
}

PHRASE_TRAILING_STOPWORDS = {
    "and", "are", "as", "between", "by", "for", "from", "in", "into", "is", "of", "on", "or", "the",
    "to", "with",
}

CAPITALIZED_ENTITY_STOPWORDS = {
    "Adverse", "Are", "Can", "Causes", "Does", "Indications", "Is", "Please", "Prevalence", "Symptoms",
    "Vector", "What", "Which",
}

PHRASE_INTERNAL_CUT_WORDS = {
    "associated", "called", "common", "considered", "correlated", "demonstrate", "driving", "found", "more",
    "linked", "performed", "promote", "reduce", "reduces", "restricted", "used",
}

BIOMEDICAL_KEYWORD_RE = re.compile(
    r"\b("
    r"asthma|atrophy|cancer|carcinoma|cell|cells|challenge|classification|complication|condensate|disease|"
    r"disorder|dystrophy|factor|factors|gene|infection|inhibitor|inhibitors|metastasis|method|"
    r"microbiome|mutation|myocardiopathy|neuropathy|oncogenesis|pathway|polymorphism|polymorphisms|procedure|"
    r"pancreatitis|protein|signature|signatures|software|stem|subtype|surfactant|syndrome|therapy|treatment|"
    r"tumor|variant|variants|virus|trypanosomiasis"
    r")\b",
    flags=re.I,
)


def _clean_question_phrase(value):
    if not value:
        return ""
    value = re.sub(r"\s+", " ", str(value)).strip(" \t\r\n,;:.?!()[]{}")
    value = QUESTION_LEAD_RE.sub("", value).strip(" \t\r\n,;:.?!")
    words = re.findall(r"[A-Za-z0-9][A-Za-z0-9+/-]*", value)
    while words and words[0].casefold() in PHRASE_LEADING_STOPWORDS:
        words.pop(0)
    while words and words[-1].casefold() in PHRASE_TRAILING_STOPWORDS:
        words.pop()
    for idx, word in enumerate(words):
        if idx > 0 and word.casefold() in PHRASE_INTERNAL_CUT_WORDS:
            words = words[:idx]
            break
    if not words:
        return ""
    if len(words) > 8:
        words = words[-8:]
        while words and words[0].casefold() in PHRASE_LEADING_STOPWORDS:
            words.pop(0)
    phrase = " ".join(words).strip()
    if len(phrase) < 3 or phrase.casefold() in PHRASE_LEADING_STOPWORDS:
        return ""
    return phrase


def _acronym(words):
    letters = []
    for word in words:
        pieces = re.split(r"[-/]+", word)
        for piece in pieces:
            if piece:
                letters.append(piece[0].casefold())
    return "".join(letters)


def _best_long_form_for_short_form(words, short_form):
    short_norm = re.sub(r"[^A-Za-z0-9]", "", short_form).casefold()
    if not short_norm:
        return ""
    max_words = min(10, len(words))
    best = ""
    for size in range(1, max_words + 1):
        tail = words[-size:]
        candidate = _clean_question_phrase(" ".join(tail))
        if not candidate:
            continue
        candidate_words = re.findall(r"[A-Za-z0-9][A-Za-z0-9+/-]*", candidate)
        acronym = _acronym(candidate_words)
        if acronym == short_norm or short_norm.startswith(acronym):
            best = candidate
            break
    return best


def _extract_parenthetical_pairs(text):
    pairs = []
    if not text:
        return pairs
    for match in re.finditer(r"\(([A-Za-z0-9][A-Za-z0-9\-/]{1,15})\)", text):
        short_form = match.group(1).strip()
        prefix = text[:match.start()]
        words = re.findall(r"[A-Za-z0-9][A-Za-z0-9+/-]*", prefix[-160:])
        long_form = _best_long_form_for_short_form(words, short_form) or _clean_question_phrase(" ".join(words[-8:]))
        if long_form and short_form:
            pairs.append((long_form, short_form))
    return pairs


def _extract_biomedical_tokens(text):
    tokens = []
    if not text:
        return tokens
    for token in re.findall(r"\b[A-Z][A-Z0-9/-]{1,15}\b", text):
        token = token.strip("-/")
        if token not in {"DNA", "RNA"}:
            tokens.append(token)
    for token in re.findall(r"\b[A-Z]{1,8}[0-9][A-Z0-9-]{0,8}\b", text):
        tokens.append(token)
    for token in re.findall(r"\b[A-Z]{2,}[a-z]*[0-9]?[A-Za-z0-9-]*\b", text):
        token = token.strip("-/")
        if len(token) <= 16:
            tokens.append(token)
    return tokens


def _extract_capitalized_entities(text):
    entities = []
    if not text:
        return entities
    for match in re.finditer(r"\b[A-Z][a-z][A-Za-z0-9-]*(?:\s+[A-Z][a-z][A-Za-z0-9-]*){0,4}\b", text):
        value = match.group(0).strip()
        if value in CAPITALIZED_ENTITY_STOPWORDS:
            continue
        words = value.split()
        if all(word in CAPITALIZED_ENTITY_STOPWORDS for word in words):
            continue
        entities.append(value)
    return entities


def _extract_question_phrases(text):
    phrases = []
    if not text:
        return phrases
    cleaned_text = re.sub(r"\([^)]{1,40}\)", " ", text)
    cleaned_text = re.sub(r"\b(AKA|aka)\b.*?(?=[,.;?]|$)", " ", cleaned_text)
    cue_patterns = [
        r"\b(?:prevalence|symptoms|indications|causes?|vector|role|complications?)\s+of\s+([^?.;,]+)",
        r"\b(?:adverse effects?|side effects?)\s+of\s+([^?.;,]+)",
        r"\b(?:treatments?|therapies|surgical treatments?)\s+for\s+([^?.;,]+)",
        r"\b(?:caused by|associated with|linked to|correlated with|performed on|patients with|risk of|use of)\s+([^?.;,]+)",
        r"\bwhat are\s+([^?.;,]+)",
        r"\bin\s+([^?.;,]+(?:cancer|disease|syndrome|infection|metastasis|population|cells?|states?)[^?.;,]*)",
        r"\bfor\s+([^?.;,]+(?:cancer|disease|syndrome|infection|therapy|trials?|cells?)[^?.;,]*)",
    ]
    for pattern in cue_patterns:
        for match in re.finditer(pattern, cleaned_text, flags=re.I):
            raw_phrase = match.group(1)
            pieces = re.split(r"\s+(?:and|or)\s+", raw_phrase, flags=re.I)
            for piece in pieces if len(pieces) > 1 else [raw_phrase]:
                phrase = _clean_question_phrase(piece)
                if phrase:
                    phrases.append(phrase)
    for match in re.finditer(
        r"\b[A-Za-z0-9][A-Za-z0-9+/-]*(?:\s+[A-Za-z0-9][A-Za-z0-9+/-]*){0,6}\b",
        cleaned_text,
    ):
        raw = match.group(0)
        phrase = _clean_question_phrase(raw)
        if not phrase:
            continue
        has_signal = bool(BIOMEDICAL_KEYWORD_RE.search(phrase) or re.search(r"[A-Z]{2,}|\d", phrase))
        if not has_signal:
            continue
        token_count = len(phrase.split())
        if token_count == 1 and not re.search(r"[A-Z]{2,}|\d", phrase):
            continue
        phrases.append(phrase)
    for entity in _extract_capitalized_entities(cleaned_text):
        phrase = _clean_question_phrase(entity)
        if phrase:
            phrases.append(phrase)
    return phrases


def _split_seed_entities(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        raw = []
        for item in value:
            raw.extend(_split_seed_entities(item))
        return raw
    return [x.strip() for x in re.split(r"[,;\n\r\t]+", str(value)) if x.strip()]


def expand_query_df(
    question=None,
    entities=None,
    entity_type=None,
    lexicon=None,
    profile="bioasq",
    max_terms=30,
):
    """Create transparent query expansion terms for KG search."""
    text = " ".join(x for x in [question or "", " ".join(_split_seed_entities(entities))] if x).strip()
    lowered = text.casefold()
    rows = []

    for entity in _split_seed_entities(entities):
        seed_type = _curated_type_for_term(entity) or entity_type
        _add_expansion(rows, entity, "seed", "user", 1.0, "Original query seed.", seed_type)

    for long_form, short_form in _extract_parenthetical_pairs(text):
        _add_expansion(rows, long_form, "parenthetical_long_form", "question", 0.93, f"Expanded from parenthetical pair with {short_form}.", entity_type)
        _add_expansion(rows, short_form, "parenthetical_short_form", "question", 0.91, f"Expanded from parenthetical pair with {long_form}.", entity_type)

    for token in _extract_biomedical_tokens(text):
        token_type = _curated_type_for_term(token) or entity_type
        _add_expansion(rows, token, "biomedical_token", "question", 0.78, "Uppercase/gene-like biomedical token found in question.", token_type)

    for phrase in _extract_question_phrases(text):
        phrase_type = _curated_type_for_term(phrase) or entity_type
        _add_expansion(rows, phrase, "question_phrase", "question", 0.64, "Biomedical noun phrase fallback extracted from question text.", phrase_type)

    normalized_text = _normalize_label(text)
    for key, terms in CURATED_QUERY_EXPANSIONS.items():
        if key in normalized_text:
            for item in terms:
                term, term_type = item if isinstance(item, (list, tuple)) else (item, entity_type)
                _add_expansion(rows, term, "curated_alias", "curated", 0.87, f"Curated expansion for {key}.", term_type)

    if profile and profile.lower() != "none":
        for item in PROFILE_QUERY_EXPANSIONS:
            if re.search(item["pattern"], lowered, flags=re.I):
                for term_item in item["terms"]:
                    term, term_type = term_item if isinstance(term_item, (list, tuple)) else (term_item, entity_type)
                    _add_expansion(rows, term, item["kind"], profile, item["weight"], item["reason"], term_type)

        for pattern, kind, term, reason in QUERY_INTENT_HINTS:
            if re.search(pattern, lowered, flags=re.I):
                _add_expansion(rows, term, kind, "intent", 0.5, reason, "")

    records, alias_map = _load_lexicon(lexicon)
    if records:
        for entity in _split_seed_entities(entities) + _extract_biomedical_tokens(text):
            name, entity_id, score, method = _ground_one(entity, entity_type or "", records, alias_map, 0.78)
            if method != "unmatched":
                _add_expansion(rows, name, "lexicon_grounded", "lexicon", 0.95 * score, f"{method}; id={entity_id}.", entity_type)
                for record in alias_map.get(_normalize_label(entity), []):
                    for alias in record.get("aliases", [])[:10]:
                        _add_expansion(rows, alias, "lexicon_alias", "lexicon", 0.9, f"Alias for grounded entity {name}.", record.get("entity_type", ""))

    return _dedupe_expansion_rows(rows, max_terms=max_terms)


def cmd_expand_query(args):
    entities = _split_seed_entities(args.entity)
    if args.entity_file:
        entities.extend(Path(args.entity_file).read_text(encoding="utf-8").splitlines())
    expanded = expand_query_df(
        question=args.question,
        entities=entities,
        entity_type=args.entity_type,
        lexicon=args.lexicon,
        profile=args.profile,
        max_terms=args.max_terms,
    )
    if args.format == "json":
        payload = {
            "query_terms": query_terms_from_expansion(expanded, args.entity_type, args.max_terms),
            "hints": [r.to_dict() for _, r in expanded.iterrows() if r["kind"] in {"target_type_hint", "relation_hint"}],
            "terms": expanded.to_dict(orient="records"),
        }
        text = json.dumps(payload, ensure_ascii=False, indent=2)
        _write_or_print(text, args.output)
    elif args.format == "txt":
        text = "\n".join(query_terms_from_expansion(expanded, args.entity_type, args.max_terms))
        _write_or_print(text, args.output)
    elif args.output:
        expanded.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(expanded.to_csv(index=False))


def _ground_one(label, entity_type, records, alias_map, min_score):
    label = "" if pd.isna(label) else str(label).strip()
    entity_type = "" if pd.isna(entity_type) else str(entity_type).strip()
    if not label:
        return "", "", 0.0, "missing"
    if not records:
        return label, label, 1.0, "self"

    key = _normalize_label(label)
    candidates = alias_map.get(key, [])
    if entity_type:
        typed = [x for x in candidates if not x["entity_type"] or x["entity_type"].casefold() == entity_type.casefold()]
        candidates = typed or candidates
    if candidates:
        best = candidates[0]
        return best["canonical_name"], best["canonical_id"], 1.0, "exact_alias"

    best = None
    best_score = 0.0
    for record in records:
        if entity_type and record["entity_type"] and record["entity_type"].casefold() != entity_type.casefold():
            continue
        score = SequenceMatcher(None, key, _normalize_label(record["canonical_name"])).ratio()
        if score > best_score:
            best = record
            best_score = score
    if best and best_score >= min_score:
        return best["canonical_name"], best["canonical_id"], float(best_score), "fuzzy_name"
    return label, label, 0.5 if records else 1.0, "unmatched"


def cmd_ground(args):
    evidence = _load_evidence_units(args.input)
    records, alias_map = _load_lexicon(args.lexicon)
    for side in ["source", "target"]:
        names = []
        ids = []
        scores = []
        methods = []
        type_col = f"{side}_type"
        for _, row in evidence.iterrows():
            name, entity_id, score, method = _ground_one(
                row.get(side),
                row.get(type_col, ""),
                records,
                alias_map,
                args.min_score,
            )
            names.append(name)
            ids.append(entity_id)
            scores.append(score)
            methods.append(method)
        evidence[f"{side}_canonical_name"] = names
        evidence[f"{side}_canonical_id"] = ids
        evidence[f"{side}_grounding_score"] = scores
        evidence[f"{side}_grounding_method"] = methods

    if args.output:
        evidence.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(evidence.to_csv(index=False))


RELATION_WEIGHTS = {
    "marker_of": 1.0,
    "produce": 1.0,
    "develope_into": 1.0,
    "pharmaco_effect_on": 1.0,
    "causally_related_to": 0.95,
    "regulated": 0.9,
    "interact_with": 0.85,
    "participate_in": 0.75,
    "functionally_related_to": 0.75,
    "is_a": 0.55,
}


def _relation_support_weight(row):
    relation = str(row.get("relation", "")).strip()
    if _edge_is_co_mentioned(row):
        # Co-mentioned/co-occurrence edges are candidate evidence only.  Even if
        # metadata is rich, they should not rank like direct curated assertions.
        detail_score = _to_float(row.get("edge_detail_score"), None)
        if detail_score is None:
            _, detail_score, _ = _co_mentioned_edge_review(row)
        return round(min(0.35, 0.10 + 0.25 * max(0.0, min(1.0, detail_score))), 6)
    if relation in MAPPING_RELATIONS:
        return 0.25
    if relation in NON_ASSERTIVE_RELATIONS:
        return 0.0
    return RELATION_WEIGHTS.get(relation, 0.5)

def _to_float(value, default=None):
    try:
        if pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default


def _bounded(value, default=0.0):
    if value is None or pd.isna(value):
        return default
    if value > 1.0 and value <= 100.0:
        value = value / 100.0
    return max(0.0, min(1.0, float(value)))


NEGATIVE_EVIDENCE_PATTERNS = [
    r"\bno evidence\b",
    r"\bnot associated\b",
    r"\bnot related\b",
    r"\black of association\b",
    r"\black of evidence\b",
    r"\babsence of\b",
    r"\bfailed to\b",
    r"\bfails to\b",
    r"\bnegative evidence\b",
]

EFFECT_DIRECTION_PATTERNS = {
    "inhibitory": [
        r"\binhibit",
        r"\bsuppress",
        r"\bdecreas",
        r"\bdownregulat",
        r"\brepress",
        r"\bantagon",
    ],
    "activating": [
        r"\bactivat",
        r"\bincreas",
        r"\bupregulat",
        r"\binduc",
        r"\bagonist",
    ],
}


def infer_polarity(row):
    for col in ["evidence_polarity", "polarity", "directionality"]:
        if col in row and pd.notna(row.get(col)):
            text = str(row.get(col)).casefold()
            if text in {"negative", "neg", "-1", "false", "contradict"}:
                return "negative"
            if text in {"positive", "pos", "1", "true", "support"}:
                return "positive"
    haystack = " ".join(str(row.get(col, "")) for col in ["relation", "evidence_text", "condition"]).casefold()
    if any(re.search(pattern, haystack) for pattern in NEGATIVE_EVIDENCE_PATTERNS):
        return "negative"
    return "positive"


def infer_effect_direction(row):
    for col in ["effect_direction", "direction"]:
        if col in row and pd.notna(row.get(col)):
            text = str(row.get(col)).strip().casefold()
            if text in {"inhibitory", "activating", "neutral", "unknown"}:
                return text
    haystack = " ".join(str(row.get(col, "")) for col in ["relation", "evidence_text", "condition"]).casefold()
    for direction, patterns in EFFECT_DIRECTION_PATTERNS.items():
        if any(re.search(pattern, haystack) for pattern in patterns):
            return direction
    return "unknown"


def score_evidence_df(evidence):
    evidence = _ensure_edge_support_columns(evidence)
    max_info = pd.to_numeric(evidence.get("info_source_count", pd.Series(dtype=float)), errors="coerce").max()
    max_info = max_info if pd.notna(max_info) and max_info > 0 else 1.0

    scores = []
    confidence_scores = []
    info_scores = []
    relation_scores = []
    grounding_scores = []
    snippet_scores = []
    polarities = []
    effect_directions = []
    score_caps = []

    for _, row in evidence.iterrows():
        confidence = _bounded(_to_float(row.get("confidence"), 0.5), 0.5)
        info = _to_float(row.get("info_source_count"), 0.0) or 0.0
        info_score = math.log1p(max(info, 0.0)) / math.log1p(max_info)
        relation_score = _relation_support_weight(row)
        source_grounding = _bounded(_to_float(row.get("source_grounding_score"), 1.0), 1.0)
        target_grounding = _bounded(_to_float(row.get("target_grounding_score"), 1.0), 1.0)
        grounding_score = (source_grounding + target_grounding) / 2.0
        snippet = str(row.get("evidence_text", "")).strip()
        snippet_score = 1.0 if snippet and snippet.casefold() != "nan" else 0.7
        polarity = infer_polarity(row)
        effect_direction = infer_effect_direction(row)
        raw_score = (
            0.30 * confidence
            + 0.25 * info_score
            + 0.20 * relation_score
            + 0.15 * grounding_score
            + 0.10 * snippet_score
        )
        strength = str(row.get("evidence_strength", ""))
        if strength.startswith("weak_co_mentioned") or _edge_is_co_mentioned(row):
            cap = 0.35
        elif strength in {"mapping_only", "not_biological_support"}:
            cap = 0.25
        elif strength.startswith("weak"):
            cap = 0.50
        else:
            cap = 1.0
        score = min(raw_score, cap)
        confidence_scores.append(confidence)
        info_scores.append(info_score)
        relation_scores.append(relation_score)
        grounding_scores.append(grounding_score)
        snippet_scores.append(snippet_score)
        polarities.append(polarity)
        effect_directions.append(effect_direction)
        score_caps.append(cap)
        scores.append(round(score, 6))

    evidence["confidence_score"] = confidence_scores
    evidence["info_source_score"] = info_scores
    evidence["relation_weight"] = relation_scores
    evidence["grounding_score"] = grounding_scores
    evidence["snippet_score"] = snippet_scores
    evidence["evidence_polarity"] = polarities
    evidence["effect_direction"] = effect_directions
    evidence["score_cap"] = score_caps
    evidence["evidence_score"] = scores
    return evidence

def _aggregate_scored(scored, mode):
    if mode == "none":
        return scored
    if mode == "target":
        group_cols = ["target"]
    elif mode == "source-target":
        group_cols = ["source", "target"]
    elif mode == "triple":
        group_cols = ["source", "relation", "target"]
    else:
        raise ValueError(f"Unsupported aggregate mode: {mode}")

    rows = []
    for key, group in scored.groupby(group_cols, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        row = dict(zip(group_cols, key))
        positive = group[group["evidence_polarity"] == "positive"]
        negative = group[group["evidence_polarity"] == "negative"]
        weak_mask = group.get("manual_review_required", pd.Series(False, index=group.index)).astype(bool)
        co_mask = group.apply(_edge_is_co_mentioned, axis=1)
        row["evidence_count"] = len(group)
        row["positive_count"] = len(positive)
        row["negative_count"] = len(negative)
        row["weak_or_review_count"] = int(weak_mask.sum())
        row["co_mentioned_count"] = int(co_mask.sum())
        row["support_score"] = positive["evidence_score"].sum()
        row["negative_score"] = negative["evidence_score"].sum()
        row["net_support_score"] = row["support_score"] - row["negative_score"]
        row["mean_evidence_score"] = group["evidence_score"].mean()
        row["sources"] = ", ".join(sorted(set(group["source"].dropna().astype(str)))[:25])
        row["relations"] = ", ".join(sorted(set(group["relation"].dropna().astype(str)))[:25])
        if "evidence_strength" in group:
            row["evidence_strengths"] = ", ".join(sorted(set(group["evidence_strength"].dropna().astype(str)))[:25])
        if "source_quality" in group:
            row["source_qualities"] = ", ".join(sorted(set(group["source_quality"].dropna().astype(str)))[:25])
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["net_support_score", "evidence_count"], ascending=[False, False])

def cmd_score(args):
    evidence = _load_evidence_units(args.input)
    scored = score_evidence_df(evidence)
    result = _aggregate_scored(scored, args.aggregate)
    if args.output:
        result.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(result.to_csv(index=False))


def cmd_conflicts(args):
    evidence = _load_evidence_units(args.input)
    if "evidence_score" not in evidence.columns or "evidence_polarity" not in evidence.columns:
        evidence = score_evidence_df(evidence)

    group_cols = ["target"] if args.by == "target" else ["source", "target"]
    if args.by == "triple":
        group_cols = ["source", "relation", "target"]

    rows = []
    for key, group in evidence.groupby(group_cols, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        positive = group[group["evidence_polarity"] == "positive"]
        negative = group[group["evidence_polarity"] == "negative"]
        has_conflict = len(positive) > 0 and len(negative) > 0
        if args.negative_only and len(negative) == 0:
            continue
        if not args.negative_only and not has_conflict:
            continue
        row = dict(zip(group_cols, key))
        row["positive_count"] = len(positive)
        row["negative_count"] = len(negative)
        row["positive_score"] = positive["evidence_score"].sum() if "evidence_score" in positive else len(positive)
        row["negative_score"] = negative["evidence_score"].sum() if "evidence_score" in negative else len(negative)
        row["conflict_score"] = min(row["positive_score"], row["negative_score"])
        row["relations"] = ", ".join(sorted(set(group["relation"].dropna().astype(str))))
        snippets = negative["evidence_text"].dropna().astype(str).head(args.max_snippets)
        row["negative_examples"] = " | ".join([x[:200] for x in snippets if x and x.casefold() != "nan"])
        rows.append(row)

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values("conflict_score", ascending=False)
    if args.output:
        result.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(result.to_csv(index=False))


def _edge_data(graph, left, right):
    data = graph.get_edge_data(left, right) or graph.get_edge_data(right, left) or {}
    if isinstance(data, dict) and "relation" not in data and data:
        first = next(iter(data.values()))
        if isinstance(first, dict):
            data = first
    return data or {}


def _edge_score_lookup(evidence_path):
    if not evidence_path:
        return {}
    evidence = _load_evidence_units(evidence_path)
    if "evidence_score" not in evidence.columns:
        evidence = score_evidence_df(evidence)
    lookup = {}
    for _, row in evidence.iterrows():
        score = _to_float(row.get("evidence_score"), 0.5) or 0.5
        relation = str(row.get("relation", ""))
        for key in [
            (_normalize_label(row.get("source")), _normalize_label(row.get("target")), relation),
            (_normalize_label(row.get("target")), _normalize_label(row.get("source")), relation),
            (_normalize_label(row.get("source")), _normalize_label(row.get("target")), ""),
            (_normalize_label(row.get("target")), _normalize_label(row.get("source")), ""),
        ]:
            lookup[key] = max(score, lookup.get(key, 0.0))
    return lookup


def _score_graph_edge(left, right, data, lookup):
    relation = str(data.get("relation", ""))
    key = (_normalize_label(left), _normalize_label(right), relation)
    fallback = (_normalize_label(left), _normalize_label(right), "")
    if key in lookup:
        return lookup[key]
    if fallback in lookup:
        return lookup[fallback]
    attr_score = _to_float(data.get("evidence_score"), None)
    if attr_score is not None:
        return _bounded(attr_score, 0.5)
    weight = _to_float(data.get("weight"), None)
    if weight is not None:
        return max(0.1, min(1.0, weight / 10.0 if weight > 1 else weight))
    row = {"source": left, "target": right, **data}
    return _relation_support_weight(row)


def cmd_path_score(args):
    graph = nx.read_graphml(args.graph)
    if args.source not in graph:
        raise ValueError(f"source not found in graph: {args.source}")
    if args.target not in graph:
        raise ValueError(f"target not found in graph: {args.target}")
    lookup = _edge_score_lookup(args.evidence)

    rows = []
    for path in nx.all_simple_paths(graph.to_undirected(), args.source, args.target, cutoff=args.max_depth):
        edge_scores = []
        relations = []
        for left, right in zip(path[:-1], path[1:]):
            data = _edge_data(graph, left, right)
            relations.append(str(data.get("relation", "unknown")))
            edge_scores.append(_score_graph_edge(left, right, data, lookup))
        if not edge_scores:
            continue
        mean_score = sum(edge_scores) / len(edge_scores)
        path_score = mean_score / math.sqrt(len(edge_scores))
        rows.append({
            "source": args.source,
            "target": args.target,
            "length": len(edge_scores),
            "path_score": round(path_score, 6),
            "mean_edge_score": round(mean_score, 6),
            "path": " -> ".join(path),
            "relations": " -> ".join(relations),
            "edge_scores": " -> ".join(f"{x:.3f}" for x in edge_scores),
        })
    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values(["path_score", "mean_edge_score"], ascending=False)
    if args.output:
        result.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(result.to_csv(index=False))


def _store_connect(db_path):
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    return sqlite3.connect(db_path)


def _sanitize_sql_identifier(name):
    return re.sub(r"[^a-zA-Z0-9_]", "_", str(name))


def _ensure_store_table(conn, df):
    existing = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table' AND name='evidence'", conn)
    safe_cols = {_sanitize_sql_identifier(col): col for col in df.columns}
    if existing.empty:
        empty = pd.DataFrame(columns=list(safe_cols.keys()))
        empty.to_sql("evidence", conn, if_exists="replace", index=False)
    pragma = pd.read_sql_query("PRAGMA table_info(evidence)", conn)
    existing_cols = set(pragma["name"]) if not pragma.empty else set()
    for safe_col in safe_cols:
        if safe_col not in existing_cols:
            conn.execute(f'ALTER TABLE evidence ADD COLUMN "{safe_col}" TEXT')
    conn.commit()
    return safe_cols


def cmd_store_init(args):
    with _store_connect(args.db) as conn:
        conn.execute(
            "CREATE TABLE IF NOT EXISTS evidence ("
            "evidence_id TEXT, source TEXT, relation TEXT, target TEXT, "
            "source_type TEXT, target_type TEXT, evidence_score REAL, "
            "evidence_polarity TEXT, provenance_file TEXT)"
        )
        conn.execute("CREATE INDEX IF NOT EXISTS idx_evidence_target ON evidence(target)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_evidence_source ON evidence(source)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_evidence_relation ON evidence(relation)")
        conn.commit()
    print(f"Initialized {args.db}")


def cmd_store_add(args):
    evidence = _load_evidence_units(args.input)
    with _store_connect(args.db) as conn:
        safe_cols = _ensure_store_table(conn, evidence)
        stored = evidence.rename(columns={original: safe for safe, original in safe_cols.items()})
        stored.to_sql("evidence", conn, if_exists="append", index=False)
    print(f"Added {len(evidence)} rows to {args.db}")


def _store_where(args, params):
    clauses = []
    for col in ["source", "target", "relation"]:
        value = getattr(args, col, None)
        if value:
            clauses.append(f"{col} LIKE ?")
            params.append(f"%{value}%")
    return " WHERE " + " AND ".join(clauses) if clauses else ""


def cmd_store_query(args):
    params = []
    where = _store_where(args, params)
    limit = max(args.limit, 1)
    query = f"SELECT * FROM evidence{where} LIMIT {limit}"
    with _store_connect(args.db) as conn:
        result = pd.read_sql_query(query, conn, params=params)
    if args.output:
        result.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(result.to_csv(index=False))


def cmd_store_export(args):
    with _store_connect(args.db) as conn:
        result = pd.read_sql_query("SELECT * FROM evidence", conn)
    if args.output:
        result.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(result.to_csv(index=False))


def cmd_store_stats(args):
    with _store_connect(args.db) as conn:
        total = pd.read_sql_query("SELECT COUNT(*) AS n FROM evidence", conn)["n"].iloc[0]
        stats = {"rows": int(total)}
        for col in ["source", "target", "relation"]:
            try:
                stats[f"distinct_{col}"] = int(pd.read_sql_query(f"SELECT COUNT(DISTINCT {col}) AS n FROM evidence", conn)["n"].iloc[0])
            except Exception:
                stats[f"distinct_{col}"] = None
    print(json.dumps(stats, indent=2, ensure_ascii=False))


def _prediction_column(df, preferred=None):
    if preferred:
        return preferred
    return _first_existing(df, ["path.2.name", "target", "prediction", "name", "entity"]) or df.columns[0]


def _score_column(df, preferred=None):
    if preferred:
        return preferred
    return _first_existing(df, ["net_support_score", "support_score", "evidence_score", "path.2.info_source_count"])


def cmd_eval(args):
    pred = pd.read_csv(args.predicted)
    ref = pd.read_csv(args.reference)
    pred_col = _prediction_column(pred, args.predicted_col)
    ref_col = _prediction_column(ref, args.reference_col)
    score_col = _score_column(pred, args.score_col)

    ranked = pred.copy()
    if score_col and score_col in ranked.columns:
        ranked = ranked.sort_values(score_col, ascending=False)
    elif "path.2.rank_voting" in ranked.columns:
        ranked = ranked.sort_values("path.2.rank_voting", ascending=True)

    predictions = [str(x) for x in ranked[pred_col].dropna()]
    gold = set(str(x) for x in ref[ref_col].dropna())
    if not gold:
        raise ValueError("Reference set is empty.")

    evidence_targets = set()
    evidence_scores = {}
    if args.evidence:
        evidence = _load_evidence_units(args.evidence)
        if "evidence_score" not in evidence.columns:
            evidence = score_evidence_df(evidence)
        evidence_targets = set(str(x) for x in evidence["target"].dropna())
        tmp = evidence.copy()
        tmp["evidence_score"] = pd.to_numeric(tmp["evidence_score"], errors="coerce")
        evidence_scores = tmp.groupby("target")["evidence_score"].mean().to_dict()

    rows = []
    for k in args.k:
        top = predictions[:k]
        hits = [x for x in top if x in gold]
        row = {
            "k": k,
            "predictions": len(top),
            "reference_size": len(gold),
            "hits": len(hits),
            "precision_at_k": len(hits) / len(top) if top else 0.0,
            "recall_at_k": len(hits) / len(gold),
            "hit_at_k": 1.0 if hits else 0.0,
        }
        if evidence_targets:
            covered = [x for x in top if x in evidence_targets]
            row["evidence_coverage_at_k"] = len(covered) / len(top) if top else 0.0
            scored_values = [evidence_scores[x] for x in top if x in evidence_scores and pd.notna(evidence_scores[x])]
            row["mean_evidence_score_at_k"] = sum(scored_values) / len(scored_values) if scored_values else None
        rows.append(row)

    reciprocal_rank = 0.0
    precisions = []
    hit_count = 0
    for idx, pred_name in enumerate(predictions, start=1):
        if pred_name in gold:
            hit_count += 1
            precisions.append(hit_count / idx)
            if reciprocal_rank == 0.0:
                reciprocal_rank = 1.0 / idx
    summary = {
        "mrr": reciprocal_rank,
        "average_precision": sum(precisions) / len(gold) if precisions else 0.0,
    }
    metrics = pd.DataFrame(rows)
    for key, value in summary.items():
        metrics[key] = value

    if args.format == "json":
        text = json.dumps({"metrics": metrics.to_dict(orient="records")}, indent=2, ensure_ascii=False)
        _write_or_print(text, args.output)
    elif args.output:
        metrics.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(metrics.to_csv(index=False))


def cmd_evidence(args):
    df = pd.read_csv(args.input)
    df = _filter_df(df, args)
    if df.empty:
        raise ValueError("No evidence rows after filtering.")

    evidence = normalize_evidence_df(df, provenance_file=args.input)

    group_cols = [col for col in ["target", "target_type"] if col in evidence.columns]
    rows = []
    for key, group in evidence.groupby(group_cols or [evidence.index], dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        target = key[0] if group_cols else "all"
        target_type = key[1] if len(key) > 1 else None
        sources = sorted(set(group.get("source", pd.Series(dtype=str)).dropna().astype(str)))
        relations = sorted(set(group.get("relation", pd.Series(dtype=str)).dropna().astype(str)))
        support_classes = sorted(set(group.get("edge_support_class", pd.Series(dtype=str)).dropna().astype(str)))
        info_count = pd.to_numeric(group.get("info_source_count", pd.Series(dtype=float)), errors="coerce").sum()
        confidence = pd.to_numeric(group.get("confidence", pd.Series(dtype=float)), errors="coerce").mean()
        edge_detail = pd.to_numeric(group.get("edge_detail_score", pd.Series(dtype=float)), errors="coerce").mean()
        co_mentioned_count = (group.get("relation", pd.Series(dtype=str)).astype(str) == "co_mentioned").sum()
        rows.append({
            "target": target,
            "target_type": target_type,
            "match_count": len(sources),
            "info_source_count": info_count,
            "mean_confidence": confidence,
            "mean_edge_detail_score": edge_detail,
            "co_mentioned_count": co_mentioned_count,
            "relations": ", ".join(relations),
            "edge_support_classes": ", ".join(support_classes),
            "sources": ", ".join(sources[:args.max_sources]),
        })

    summary = pd.DataFrame(rows).sort_values(
        by=[c for c in ["match_count", "info_source_count"] if c in rows[0]],
        ascending=False,
    )

    if args.format == "csv":
        if args.output:
            summary.to_csv(args.output, index=False)
            print(f"Wrote {args.output}")
        else:
            print(summary.to_csv(index=False))
        return

    lines = ["# brick Evidence Summary", ""]
    for _, row in summary.head(args.top_n).iterrows():
        title = row["target"]
        if pd.notna(row.get("target_type")):
            title += f" ({row['target_type']})"
        lines.append(f"## {title}")
        lines.append(f"- Matched sources: {row['match_count']}")
        if pd.notna(row.get("info_source_count")):
            lines.append(f"- Information sources: {row['info_source_count']}")
        if pd.notna(row.get("mean_confidence")):
            lines.append(f"- Mean confidence: {row['mean_confidence']:.3f}")
        if pd.notna(row.get("mean_edge_detail_score")):
            lines.append(f"- Mean edge detail score: {row['mean_edge_detail_score']:.3f}")
        if row.get("co_mentioned_count"):
            lines.append(f"- Co-mentioned rows: {int(row['co_mentioned_count'])}")
        if row.get("relations"):
            lines.append(f"- Relations: {row['relations']}")
        if row.get("edge_support_classes"):
            lines.append(f"- Edge support classes: {row['edge_support_classes']}")
        if row.get("sources"):
            lines.append(f"- Source examples: {row['sources']}")
        target_group = evidence[evidence["target"].astype(str) == str(row["target"])]
        snippets = target_group[target_group["evidence_text"].apply(_has_text)].head(args.max_snippets)
        if len(snippets):
            lines.append("- Edge evidence:")
            for _, edge in snippets.iterrows():
                details = []
                for label, col in [
                    ("relation", "relation"),
                    ("support", "edge_support_class"),
                    ("source", "info_source"),
                    ("condition", "condition"),
                    ("confidence", "confidence"),
                ]:
                    value = edge.get(col)
                    if _has_text(value):
                        details.append(f"{label}={value}")
                prefix = f"  - [{' ; '.join(details)}]" if details else "  -"
                lines.append(f"{prefix} {str(edge.get('evidence_text'))[:300]}")
        notes = target_group[target_group["relation"].astype(str) == "co_mentioned"]["edge_review_note"].dropna().astype(str).unique()
        for note in notes[:3]:
            lines.append(f"- Co-mentioned review: {note}")
        lines.append("")

    _write_or_print("\n".join(lines), args.output)


def cmd_hypotheses(args):
    df = pd.read_csv(args.input)
    if args.target_type and "path.2.type" in df:
        df = df[df["path.2.type"] == args.target_type]
    if df.empty:
        raise ValueError("No hypothesis rows after filtering.")

    result = df.copy()
    rank_col = _first_existing(result, ["path.2.rank_voting", "path.2.enrich_pvalue", "path.2.enrich_fdr"])
    if rank_col:
        result = result.sort_values(rank_col, ascending=True)

    reason_bits = []
    for _, row in result.iterrows():
        bits = []
        if "path.2.match_count" in result:
            bits.append(f"matches={row.get('path.2.match_count')}")
        if "path.2.info_source_count" in result:
            bits.append(f"evidence={row.get('path.2.info_source_count')}")
        if "path.2.enrich_pvalue" in result:
            bits.append(f"p={row.get('path.2.enrich_pvalue')}")
        if "path.2.enrich_fdr" in result:
            bits.append(f"fdr={row.get('path.2.enrich_fdr')}")
        reason_bits.append("; ".join(bits))
    result["reasoning_summary"] = reason_bits

    cols = [
        c for c in [
            "path.2.name",
            "path.2.type",
            rank_col,
            "path.2.match_count",
            "path.2.info_source_count",
            "path.2.enrich_pvalue",
            "path.2.enrich_fdr",
            "reasoning_summary",
        ]
        if c and c in result.columns
    ]
    result = result[cols].head(args.top_n)
    if args.output:
        result.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(result.to_string(index=False))


def cmd_paths(args):
    graph = nx.read_graphml(args.graph)
    if args.source not in graph:
        raise ValueError(f"source not found in graph: {args.source}")
    if args.target not in graph:
        raise ValueError(f"target not found in graph: {args.target}")

    paths = []
    for path in nx.all_simple_paths(graph.to_undirected(), args.source, args.target, cutoff=args.max_depth):
        edges = []
        for left, right in zip(path[:-1], path[1:]):
            data = graph.get_edge_data(left, right) or graph.get_edge_data(right, left) or {}
            if isinstance(data, dict) and "relation" not in data and data:
                first = next(iter(data.values()))
                if isinstance(first, dict):
                    data = first
            edges.append(str(data.get("relation", "unknown")))
        paths.append({
            "source": args.source,
            "target": args.target,
            "length": len(path) - 1,
            "path": " -> ".join(path),
            "relations": " -> ".join(edges),
        })

    result = pd.DataFrame(paths).sort_values("length") if paths else pd.DataFrame(columns=["source", "target", "length", "path", "relations"])
    if args.output:
        result.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(result.to_string(index=False))


def cmd_compare(args):
    pred = pd.read_csv(args.predicted)
    ref = pd.read_csv(args.reference)
    pred_col = args.predicted_col or ("path.2.name" if "path.2.name" in pred else pred.columns[0])
    ref_col = args.reference_col or ("path.2.name" if "path.2.name" in ref else ref.columns[0])
    ref_values = set(ref[ref_col].dropna().astype(str))
    result = pred.copy()
    result["is_supported_by_reference"] = result[pred_col].astype(str).isin(ref_values)
    result["reference_size"] = len(ref_values)
    result["overlap_count"] = int(result["is_supported_by_reference"].sum())

    if args.output:
        result.to_csv(args.output, index=False)
        print(f"Wrote {args.output}")
    else:
        print(result.to_string(index=False))


# ---------------------------------------------------------------------------
# literature: inspect the publications behind co-mentioned edges
# ---------------------------------------------------------------------------

def _parse_list_field(value):
    """Parse KG list-valued CSV fields like "['PMID:123', 'PMID:456']"."""
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "[]"}:
        return []
    try:
        parsed = ast.literal_eval(text)
        if isinstance(parsed, (list, tuple)):
            return [str(v) for v in parsed]
    except Exception:
        pass
    return [text]


def collect_co_mentioned_literature(df):
    """Extract per-PMID evidence entries from co_mentioned query rows.

    ``path.1.info_source`` (PMID list) and ``path.1.original_relation``
    (source sentences) are aligned 1:1 in KG query output. Returns a
    deduplicated, order-preserving list of dicts with keys
    pmid / sentence / source / target.
    """
    entries = []
    rel_col = _first_existing(df, ["path.1.relation", "relation"])
    src_col = _first_existing(df, ["path.0.name", "source"])
    tgt_col = _first_existing(df, ["path.2.name", "target"])
    if not rel_col:
        return entries
    for _, row in df.iterrows():
        if not _is_co_mentioned_relation(row.get(rel_col)):
            continue
        pmids = [p for p in _parse_list_field(row.get("path.1.info_source"))
                 if "pmid" in p.lower()]
        sentences = _parse_list_field(row.get("path.1.original_relation"))
        for i, pmid in enumerate(pmids):
            entries.append({
                "pmid": pmid if pmid.upper().startswith("PMID:") else f"PMID:{pmid}",
                "sentence": sentences[i] if i < len(sentences) else "",
                "source": str(row.get(src_col, "")) if src_col else "",
                "target": str(row.get(tgt_col, "")) if tgt_col else "",
            })
    seen = set()
    uniq = []
    for e in entries:
        key = (e["pmid"], e["sentence"])
        if key in seen:
            continue
        seen.add(key)
        uniq.append(e)
    return uniq


def _pubmed_summaries(pmids):
    """Fetch title/journal/pubdate for PMIDs via NCBI E-utilities (one call).

    Returns {pmid_number: {"title":..., "journal":..., "pubdate":...}};
    empty dict on any failure (metadata is best-effort)."""
    if not pmids:
        return {}
    try:
        import requests
    except Exception:
        return {}
    try:
        resp = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "pubmed", "retmode": "json", "id": ",".join(pmids)},
            timeout=30)
        data = resp.json().get("result", {})
    except Exception:
        return {}
    out = {}
    for pmid in pmids:
        rec = data.get(pmid) or {}
        out[pmid] = {
            "title": rec.get("title", ""),
            "journal": rec.get("fulljournalname", ""),
            "pubdate": rec.get("pubdate", ""),
        }
    return out


def cmd_literature(args):
    try:
        from querygraph import query_relation_with_fallback
    except Exception as e:
        raise ValueError(f"querygraph module unavailable: {e}")
    df = query_relation_with_fallback(
        source_entity_set=[args.source],
        target_entity_set=[args.target],
        source_entity_type=args.st,
        target_entity_type=args.tt,
        multi_hop=1,
        return_type="dataframe",
    )
    if df is None or df.empty:
        raise ValueError(f"No KG edges between {args.source!r} and {args.target!r}")
    entries = collect_co_mentioned_literature(df)
    if not entries:
        raise ValueError(
            f"No co-mentioned literature edges between {args.source!r} and "
            f"{args.target!r} ({len(df)} KG rows returned; none carry PMIDs)")

    truncated = False
    if args.limit and len(entries) > args.limit:
        entries = entries[: args.limit]
        truncated = True

    meta = {}
    if not args.no_pubmed:
        meta = _pubmed_summaries(
            [e["pmid"].replace("PMID:", "") for e in entries])

    if args.format == "csv":
        rows = []
        for e in entries:
            m = meta.get(e["pmid"].replace("PMID:", ""), {})
            rows.append({**e, "title": m.get("title", ""),
                         "journal": m.get("journal", ""),
                         "pubdate": m.get("pubdate", "")})
        result = pd.DataFrame(rows)
        if args.output:
            result.to_csv(args.output, index=False)
            print(f"Wrote {args.output} ({len(result)} entries)")
        else:
            print(result.to_string(index=False))
        return

    lines = [
        f"# Literature behind co-mentioned edges: {args.source} — {args.target}",
        "",
        f"KG query returned {len(df)} rows; {len(entries)} publication–sentence pairs"
        + (f" (truncated to --limit {args.limit})" if truncated else "") + ".",
        "",
    ]
    for e in entries:
        m = meta.get(e["pmid"].replace("PMID:", ""), {})
        lines.append(f"## {e['pmid']}")
        if m.get("title"):
            lines.append(f"- Title: {m['title']}")
            lines.append(f"- Journal: {m.get('journal', '')} ({m.get('pubdate', '')})")
        if e["sentence"]:
            lines.append(f"- Sentence: \"{e['sentence']}\"")
        lines.append("")
    report = "\n".join(lines)
    if args.output:
        Path(args.output).write_text(report)
        print(f"Wrote {args.output} ({len(entries)} entries)")
    else:
        print(report)


def build_parser():
    parser = argparse.ArgumentParser(description="brick reason - evidence and reasoning helpers")
    subparsers = parser.add_subparsers(dest="command")

    normalize = subparsers.add_parser("normalize", help="Normalize brick qr CSV into canonical evidence units")
    normalize.add_argument("-i", "--input", required=True)
    normalize.add_argument("--source")
    normalize.add_argument("--target")
    normalize.add_argument("--relation")
    normalize.add_argument("--min-confidence", type=float)
    normalize.add_argument("--provenance-file")
    normalize.add_argument("--deduplicate", action="store_true")
    normalize.add_argument("--format", choices=["csv", "jsonl"], default="csv")
    normalize.add_argument("-o", "--output")
    normalize.set_defaults(func=cmd_normalize)

    evidence = subparsers.add_parser("evidence", help="Summarize evidence from brick qr CSV")
    evidence.add_argument("-i", "--input", required=True)
    evidence.add_argument("--source")
    evidence.add_argument("--target")
    evidence.add_argument("--relation")
    evidence.add_argument("--top-n", type=int, default=10)
    evidence.add_argument("--max-sources", type=int, default=20)
    evidence.add_argument("--max-snippets", type=int, default=5)
    evidence.add_argument("--format", choices=["md", "csv"], default="md")
    evidence.add_argument("-o", "--output")
    evidence.set_defaults(func=cmd_evidence)

    hypo = subparsers.add_parser("hypotheses", help="Create hypothesis shortlist from ranked/enriched CSV")
    hypo.add_argument("-i", "--input", required=True)
    hypo.add_argument("--target-type")
    hypo.add_argument("--top-n", type=int, default=20)
    hypo.add_argument("-o", "--output")
    hypo.set_defaults(func=cmd_hypotheses)

    explain = subparsers.add_parser("explain", help="Generate an interpretability report from evidence units")
    explain.add_argument("-e", "--evidence", required=True)
    explain.add_argument("--hypotheses")
    explain.add_argument("--target")
    explain.add_argument("--top-n", type=int, default=10)
    explain.add_argument("--max-sources", type=int, default=20)
    explain.add_argument("--max-snippets", type=int, default=5)
    explain.add_argument("-o", "--output")
    explain.set_defaults(func=cmd_explain)

    ground = subparsers.add_parser("ground", help="Ground source and target entities to canonical IDs/names")
    ground.add_argument("-i", "--input", required=True)
    ground.add_argument("--lexicon", help="CSV with canonical_name/name, optional canonical_id/entity_id, entity_type, aliases")
    ground.add_argument("--min-score", type=float, default=0.84)
    ground.add_argument("-o", "--output")
    ground.set_defaults(func=cmd_ground)

    expand = subparsers.add_parser("expand-query", help="Expand question/entity terms before brick KG querying")
    expand.add_argument("-q", "--question", help="Natural-language question or query context")
    expand.add_argument("-e", "--entity", action="append", help="Seed entity or comma-separated seed entities")
    expand.add_argument("--entity-file", help="Optional file with seed entities")
    expand.add_argument("--entity-type", help="Expected entity type for lexicon matching")
    expand.add_argument("--lexicon", help="CSV with canonical_name/name, optional canonical_id/entity_id, entity_type, aliases")
    expand.add_argument("--profile", default="bioasq", help="Expansion profile: bioasq or none")
    expand.add_argument("--max-terms", type=int, default=30)
    expand.add_argument("--format", choices=["csv", "json", "txt"], default="csv")
    expand.add_argument("-o", "--output")
    expand.set_defaults(func=cmd_expand_query)

    packet = subparsers.add_parser("packet", help="Build an answer-boundary evidence packet from KG evidence")
    packet.add_argument("-i", "--input", required=True)
    packet.add_argument("--expected-answer-type")
    packet.add_argument("--scope-term", action="append", help="Question scope term; may be repeated")
    packet.add_argument("--scope-match", choices=["any", "all"], default="any")
    packet.add_argument("--candidate-term", action="append", help="Candidate anchor term; may be repeated")
    packet.add_argument("--candidate-side", choices=["source", "target", "either"], default="target")
    packet.add_argument("--weak-evidence-policy", choices=["context_only", "allow_weak_candidates"], default="context_only")
    packet.add_argument("--max-candidate-rows", type=int)
    packet.add_argument("--kg-schema")
    packet.add_argument("--max-rows", type=int, default=12, help="Maximum strong rows for prompt rendering")
    packet.add_argument("--format", choices=["json", "csv", "prompt"], default="json")
    packet.add_argument("-o", "--output")
    packet.set_defaults(func=cmd_packet)

    filter_co = subparsers.add_parser(
        "filter-co-mentioned",
        help="Drop co-mentioned query rows whose original_relation snippet does not support the edge",
    )
    filter_co.add_argument("-i", "--input", required=True, help="brick qr CSV to review")
    filter_co.add_argument("-o", "--output", required=True, help="Filtered query CSV for downstream rk/reason steps")
    filter_co.add_argument("--rejected-output", help="Optional CSV of dropped co-mentioned rows")
    filter_co.add_argument("--audit-output", help="Optional CSV with keep/drop decisions for every input row")
    filter_co.add_argument("--source", help="Optional source filter before review")
    filter_co.add_argument("--target", help="Optional target filter before review")
    filter_co.add_argument("--relation", help="Optional relation filter before review")
    filter_co.add_argument(
        "--mode",
        choices=["strict", "permissive"],
        default="strict",
        help="strict requires snippet review; permissive only drops unresolved/metadata-only co-mentioned rows",
    )
    filter_co.add_argument(
        "--min-detail-score",
        type=float,
        default=0.55,
        help="Minimum co-mentioned detail score required to keep a row in strict mode",
    )
    filter_co.add_argument(
        "--keep-bare-co-occurrence",
        action="store_true",
        help="Keep co-mentioned rows even when original_relation lacks an explicit relationship cue",
    )
    filter_co.add_argument(
        "--allow-missing-endpoint",
        action="store_true",
        help="Do not require both endpoints to appear in original_relation",
    )
    filter_co.add_argument(
        "--include-audit-columns",
        action="store_true",
        help="Keep co_mentioned_* review columns in the main filtered output",
    )
    filter_co.set_defaults(func=cmd_filter_co_mentioned)

    score = subparsers.add_parser("score", help="Score evidence units and optionally aggregate support")
    score.add_argument("-i", "--input", required=True)
    score.add_argument("--aggregate", choices=["none", "target", "source-target", "triple"], default="target")
    score.add_argument("-o", "--output")
    score.set_defaults(func=cmd_score)

    conflicts = subparsers.add_parser("conflicts", help="Detect conflict or negative evidence")
    conflicts.add_argument("-i", "--input", required=True)
    conflicts.add_argument("--by", choices=["target", "source-target", "triple"], default="source-target")
    conflicts.add_argument("--negative-only", action="store_true")
    conflicts.add_argument("--max-snippets", type=int, default=3)
    conflicts.add_argument("-o", "--output")
    conflicts.set_defaults(func=cmd_conflicts)

    paths = subparsers.add_parser("paths", help="Enumerate graph paths between two entities")
    paths.add_argument("-g", "--graph", required=True)
    paths.add_argument("--source", required=True)
    paths.add_argument("--target", required=True)
    paths.add_argument("--max-depth", type=int, default=4)
    paths.add_argument("-o", "--output")
    paths.set_defaults(func=cmd_paths)

    path_score = subparsers.add_parser("path-score", help="Score mechanism paths between two graph nodes")
    path_score.add_argument("-g", "--graph", required=True)
    path_score.add_argument("--source", required=True)
    path_score.add_argument("--target", required=True)
    path_score.add_argument("-e", "--evidence")
    path_score.add_argument("--max-depth", type=int, default=4)
    path_score.add_argument("-o", "--output")
    path_score.set_defaults(func=cmd_path_score)

    compare = subparsers.add_parser("compare", help="Compare predictions with reference entities")
    compare.add_argument("-p", "--predicted", required=True)
    compare.add_argument("-r", "--reference", required=True)
    compare.add_argument("--predicted-col")
    compare.add_argument("--reference-col")
    compare.add_argument("-o", "--output")
    compare.set_defaults(func=cmd_compare)

    store = subparsers.add_parser("store", help="Manage a local SQLite evidence store")
    store_sub = store.add_subparsers(dest="store_command", required=True)

    store_init = store_sub.add_parser("init", help="Initialize an evidence store")
    store_init.add_argument("--db", required=True)
    store_init.set_defaults(func=cmd_store_init)

    store_add = store_sub.add_parser("add", help="Add evidence units to a store")
    store_add.add_argument("--db", required=True)
    store_add.add_argument("-i", "--input", required=True)
    store_add.set_defaults(func=cmd_store_add)

    store_query = store_sub.add_parser("query", help="Query an evidence store")
    store_query.add_argument("--db", required=True)
    store_query.add_argument("--source")
    store_query.add_argument("--target")
    store_query.add_argument("--relation")
    store_query.add_argument("--limit", type=int, default=50)
    store_query.add_argument("-o", "--output")
    store_query.set_defaults(func=cmd_store_query)

    store_export = store_sub.add_parser("export", help="Export all evidence rows from a store")
    store_export.add_argument("--db", required=True)
    store_export.add_argument("-o", "--output")
    store_export.set_defaults(func=cmd_store_export)

    store_stats = store_sub.add_parser("stats", help="Summarize an evidence store")
    store_stats.add_argument("--db", required=True)
    store_stats.set_defaults(func=cmd_store_stats)

    evaluate = subparsers.add_parser("eval", help="Evaluate predictions against a reference set")
    evaluate.add_argument("-p", "--predicted", required=True)
    evaluate.add_argument("-r", "--reference", required=True)
    evaluate.add_argument("-e", "--evidence")
    evaluate.add_argument("--predicted-col")
    evaluate.add_argument("--reference-col")
    evaluate.add_argument("--score-col")
    evaluate.add_argument("--k", type=int, nargs="+", default=[1, 5, 10, 20])
    evaluate.add_argument("--format", choices=["csv", "json"], default="csv")
    evaluate.add_argument("-o", "--output")
    evaluate.set_defaults(func=cmd_eval)

    lit = subparsers.add_parser(
        "literature",
        help="Inspect publications behind co-mentioned edges between two entities")
    lit.add_argument("-s", "--source", required=True, help="Source entity (standard KG name)")
    lit.add_argument("-t", "--target", required=True, help="Target entity (standard KG name)")
    lit.add_argument("--st", "--source-type", dest="st", help="Source entity type (e.g. Gene)")
    lit.add_argument("--tt", "--target-type", dest="tt", help="Target entity type (e.g. Disease)")
    lit.add_argument("--limit", type=int, default=50,
                     help="Max publication-sentence pairs to report (default 50; 0 = all)")
    lit.add_argument("--no-pubmed", action="store_true",
                     help="Skip PubMed metadata lookup (title/journal/date)")
    lit.add_argument("--format", choices=["md", "csv"], default="md")
    lit.add_argument("-o", "--output")
    lit.set_defaults(func=cmd_literature)

    return parser


def _cli_main():
    parser = build_parser()
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit(1)
    try:
        args.func(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _cli_main()
