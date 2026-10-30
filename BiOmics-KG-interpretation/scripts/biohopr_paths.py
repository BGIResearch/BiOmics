"""BioHopR PathAnswerSet utilities for BiOmics-BRICK.

This module keeps BioHopR graph traversal inside the BRICK skill layer.
Benchmark runners should only load records, pass TaskContract objects here,
and write/evaluate the rendered answer files.
"""

from __future__ import annotations

import json
import math
import re
import time
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any, MutableMapping

import pandas as pd


GOLD_FILES = {
    "1hop_single": "questions_1hop_single_gold.json",
    "2hop_single": "questions_2hop_single_gold.json",
    "1hop_full": "questions_1hop_full_gold.json",
    "2hop_full": "questions_2hop_full_gold.json",
}

WEAK_RELATIONS = {"co_mentioned", "mentioned", "no_relation", ""}
STRONG_RELATIONS = {
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
GENERIC_DISEASE_TERMS = {
    "disease",
    "diseases",
    "syndrome",
    "disorder",
    "cancer",
    "neoplasm",
    "neoplasms",
    "tumor",
    "tumors",
    "infection",
    "inflammation",
    "genetic diseases, inborn",
    "genetic disease, inborn",
    "genetic diseases",
    "inborn genetic diseases",
}
HUB_LIKE_PREFIXES = ("hla-", "hla ", "cytochrome p450", "immunoglobulin")
TYPE_FALLBACKS = {
    "Gene": ["Gene", "Protein"],
    "Protein": ["Protein", "Gene"],
    "Chemical": ["Chemical"],
    "Disease": ["Disease", "Phenotype"],
    "Phenotype": ["Phenotype", "Disease"],
    "Pathway": ["Pathway", "Process", "Function"],
    "Process": ["Process", "Function", "Pathway"],
    "Function": ["Function", "Process", "Pathway"],
}
CURATED_INFO_SOURCE_PREFIXES = {
    "AGRKB",
    "AGRKB_DISEASE",
    "BTO",
    "CHEBI",
    "CL",
    "CLINICALTRAILS",
    "CLINICALTRIALS",
    "CLINVAR",
    "CLINVAR_GENE",
    "CLO",
    "DOID",
    "DRUGBANK",
    "DRUGCENTRAL",
    "EFO",
    "FOODB_COMPOUND",
    "FOODB_FOOD",
    "FOODON",
    "GENE2PROTEIN",
    "GENE2PROTEIN_MM",
    "GO",
    "GOA",
    "GOA_MM",
    "HMDB",
    "HPO_ANNO",
    "HPO_D2P_ANNO",
    "HPO_G2P_ANNO",
    "INTACT",
    "IUPHAR",
    "KEGG",
    "KEGG_ANNO",
    "KEGG_ANNO_MM",
    "MEDGEN",
    "MESH",
    "MONDO",
    "NCBI_GENE",
    "NCIT",
    "OMIM",
    "OMIM_ANNO",
    "ORPHANET",
    "PANTHER_PATHWAY",
    "PANGLAODB",
    "PHARMGKB",
    "PRIMEKG",
    "PUBCHEMCID",
    "REACTOME",
    "REACTOME_ANNO",
    "REPODB",
    "RHEA",
    "SIDER",
    "SIGNOR",
    "SMPDB",
    "STITCH",
    "STRINGDB",
    "TRRUST",
    "UNIPROT_HUMAN",
    "UNIPROT_MAPPING",
    "UNIPROT_RODENTS",
}
LOW_SIGNAL_INFO_SOURCE_PREFIXES = {"PMID", "PMCID", "PUBMED", "LITERATURE"}

BIOHOPR_BENCHMARK_POLICIES: dict[str, dict[str, Any]] = {
    "balanced": {
        "description": "Strict extracted-query BioHopR main preset: extract the PathAnswerSet contract from the natural-language question, then use native answers plus score-gated KG path supplements.",
        "kg_limit": 1000,
        "max_answers": 120,
        "kg_answer_strength": "all",
        "path_co_mentioned_policy": "allow",
        "co_mentioned_fallback_min_score": 0.2,
        "answer_selection_mode": "top_k",
        "score_drop": 1.35,
        "min_answers": 3,
        "answer_policy": "native_plus_kg",
        "kg_add_max_1hop": 120,
        "kg_add_max_2hop": 120,
        "kg_add_min_score_1hop": 2.0,
        "kg_add_min_score_2hop": 4.0,
        "kg_add_allow_weak": True,
        "biohopr_source_policy": "none",
        "weak_source_rescue_max_1hop": 0,
        "weak_source_rescue_max_2hop": 0,
        "weak_source_rescue_min_score_1hop": None,
        "weak_source_rescue_min_score_2hop": None,
        "kg_info_source": "",
        "query_source": "extracted",
    },
    "strict_e2e_question_only": {
        "description": "Strict question-only BioHopR end-to-end preset for main results. Structured queries are extracted from question text; gold answers, gold hops, official bridge entities, and official query/target types are not used during inference.",
        "kg_limit": 1000,
        "max_answers": 120,
        "kg_answer_strength": "all",
        "path_co_mentioned_policy": "allow",
        "co_mentioned_fallback_min_score": 0.2,
        "answer_selection_mode": "top_k",
        "score_drop": 1.35,
        "min_answers": 3,
        "answer_policy": "native_plus_kg",
        "kg_add_max_1hop": 120,
        "kg_add_max_2hop": 120,
        "kg_add_min_score_1hop": 2.0,
        "kg_add_min_score_2hop": 4.0,
        "kg_add_allow_weak": True,
        "biohopr_source_policy": "none",
        "weak_source_rescue_max_1hop": 0,
        "weak_source_rescue_max_2hop": 0,
        "weak_source_rescue_min_score_1hop": None,
        "weak_source_rescue_min_score_2hop": None,
        "kg_info_source": "",
        "query_source": "extracted",
    },
    "sourcefiltered": {
        "description": "BioHopR diagnostic preset: question-only extracted contract plus curated source preference with limited 1-hop weak literature rescue.",
        "kg_limit": 1000,
        "max_answers": 120,
        "kg_answer_strength": "all",
        "path_co_mentioned_policy": "allow",
        "co_mentioned_fallback_min_score": 0.2,
        "answer_selection_mode": "top_k",
        "score_drop": 1.35,
        "min_answers": 3,
        "answer_policy": "native_plus_kg",
        "kg_add_max_1hop": 120,
        "kg_add_max_2hop": 120,
        "kg_add_min_score_1hop": 2.0,
        "kg_add_min_score_2hop": 4.0,
        "kg_add_allow_weak": True,
        "biohopr_source_policy": "curated_with_weak_rescue",
        "weak_source_rescue_max_1hop": 3,
        "weak_source_rescue_max_2hop": 0,
        "weak_source_rescue_min_score_1hop": 2.0,
        "weak_source_rescue_min_score_2hop": 4.0,
        "kg_info_source": "",
        "query_source": "extracted",
    },
}


def _with_policy_defaults(policy: dict[str, Any]) -> dict[str, Any]:
    out = dict(policy)
    out.setdefault("entity_linking", "off")
    out.setdefault("max_query_entities", 1)
    out.setdefault("type_fallback", False)
    out.setdefault("discover_bridge_fallback", False)
    out.setdefault("max_bridge_entities", 25)
    out.setdefault("relax_source_type_fallback", False)
    out.setdefault("relax_target_type_fallback", False)
    out.setdefault("kg_add_max_weak_1hop", out.get("kg_add_max_1hop", 0))
    out.setdefault("kg_add_max_weak_2hop", out.get("kg_add_max_2hop", 0))
    out.setdefault("entity_linking_1hop", None)
    out.setdefault("entity_linking_2hop", None)
    out.setdefault("max_query_entities_1hop", None)
    out.setdefault("max_query_entities_2hop", None)
    out.setdefault("type_fallback_1hop", False)
    out.setdefault("type_fallback_2hop", False)
    out.setdefault("discover_bridge_fallback_1hop", False)
    out.setdefault("discover_bridge_fallback_2hop", False)
    out.setdefault("hop_aware", False)
    return out


BIOHOPR_BENCHMARK_POLICIES = {
    key: _with_policy_defaults(value)
    for key, value in BIOHOPR_BENCHMARK_POLICIES.items()
}


def _clone_policy(base: str, description: str, **overrides: Any) -> dict[str, Any]:
    policy = dict(BIOHOPR_BENCHMARK_POLICIES[base])
    policy["description"] = description
    policy.update(overrides)
    return _with_policy_defaults(policy)


BIOHOPR_BENCHMARK_POLICIES["strict_e2e_baseline"] = _clone_policy(
    "strict_e2e_question_only",
    "Strict question-only baseline matching the unoptimized extracted-query main result.",
)
BIOHOPR_BENCHMARK_POLICIES["strict_e2e_entitylink"] = _clone_policy(
    "strict_e2e_question_only",
    "Strict question-only ablation with surface-form entity candidate expansion before KG querying.",
    entity_linking="surface",
    max_query_entities=8,
)
BIOHOPR_BENCHMARK_POLICIES["strict_e2e_fallback"] = _clone_policy(
    "strict_e2e_question_only",
    "Strict question-only ablation with entity expansion plus type/relaxed-type fallback for empty KG path queries.",
    entity_linking="surface",
    max_query_entities=8,
    type_fallback=True,
    discover_bridge_fallback=True,
    max_bridge_entities=25,
)
BIOHOPR_BENCHMARK_POLICIES["strict_e2e_scorecutoff"] = _clone_policy(
    "strict_e2e_question_only",
    "Strict question-only ablation with score-window KG answer selection and capped native-plus-KG merging.",
    max_answers=20,
    answer_selection_mode="score_cutoff",
    min_path_score=2.8,
    score_drop=0.8,
    kg_add_max_1hop=5,
    kg_add_max_2hop=8,
    kg_add_min_score_1hop=3.0,
    kg_add_min_score_2hop=4.8,
    kg_add_allow_weak=False,
    kg_add_max_weak_1hop=0,
    kg_add_max_weak_2hop=0,
)
BIOHOPR_BENCHMARK_POLICIES["strict_e2e_entitylink_fallback"] = _clone_policy(
    "strict_e2e_fallback",
    "Strict question-only ablation combining entity expansion with staged empty-query repair.",
)
BIOHOPR_BENCHMARK_POLICIES["strict_e2e_entitylink_fallback_scorecutoff"] = _clone_policy(
    "strict_e2e_scorecutoff",
    "Strict question-only optimized preset: entity expansion, staged empty-query repair, and capped score-cutoff KG merging.",
    entity_linking="surface",
    max_query_entities=8,
    type_fallback=True,
    discover_bridge_fallback=True,
    max_bridge_entities=25,
)
BIOHOPR_BENCHMARK_POLICIES["strict_e2e_soft_rerank"] = _clone_policy(
    "strict_e2e_fallback",
    (
        "Strict question-only soft rerank preset: keep entity expansion and discovered bridge "
        "fallback for 2-hop recall, then use top-k path ranking plus strong/weak evidence quotas "
        "instead of a hard score cutoff."
    ),
    max_answers=40,
    answer_selection_mode="top_k",
    min_path_score=None,
    score_drop=1.35,
    kg_add_max_1hop=8,
    kg_add_max_2hop=30,
    kg_add_min_score_1hop=2.5,
    kg_add_min_score_2hop=4.0,
    kg_add_allow_weak=True,
    kg_add_max_weak_1hop=0,
    kg_add_max_weak_2hop=6,
)
BIOHOPR_BENCHMARK_POLICIES["strict_e2e_hop_aware"] = _clone_policy(
    "strict_e2e_question_only",
    (
        "Strict question-only hop-aware BRICK/KG preset: use the same TaskContract -> KG "
        "traversal -> path scoring -> native/KG merge pipeline for every BioHopR item, while "
        "routing 1-hop to conservative direct retrieval and 2-hop to entity-linking plus "
        "discovered-bridge fallback."
    ),
    hop_aware=True,
    max_answers=40,
    answer_selection_mode="top_k",
    min_path_score=None,
    score_drop=1.35,
    entity_linking="off",
    max_query_entities=1,
    type_fallback=False,
    discover_bridge_fallback=False,
    entity_linking_1hop="off",
    max_query_entities_1hop=1,
    type_fallback_1hop=False,
    discover_bridge_fallback_1hop=False,
    entity_linking_2hop="surface",
    max_query_entities_2hop=8,
    type_fallback_2hop=True,
    discover_bridge_fallback_2hop=True,
    max_bridge_entities=25,
    kg_add_max_1hop=12,
    kg_add_max_2hop=30,
    kg_add_min_score_1hop=2.0,
    kg_add_min_score_2hop=4.0,
    kg_add_allow_weak=True,
    kg_add_max_weak_1hop=0,
    kg_add_max_weak_2hop=6,
)


@dataclass(frozen=True)
class PathAnswerOptions:
    max_answers: int = 20
    kg_answer_strength: str = "strong_path"
    path_co_mentioned_policy: str = "strong_first_fallback"
    co_mentioned_fallback_min_score: float = 0.2
    selection_mode: str = "top_k"
    min_path_score: float | None = None
    score_drop: float = 1.35
    min_answers: int = 3
    source_policy: str = "none"
    weak_source_rescue_max: int = 0
    weak_source_rescue_min_score: float | None = None
    weak_source_rescue_score_bonus: float = 0.0


def biohopr_benchmark_policy(name: str = "balanced") -> dict[str, Any]:
    key = norm(name).replace("-", "_")
    aliases = {
        "d_main_v2": "balanced",
        "nativepluskg": "balanced",
        "native_plus_kg": "balanced",
        "baseline": "strict_e2e_baseline",
        "strict_e2e": "strict_e2e_question_only",
        "question_only": "strict_e2e_question_only",
        "extracted": "strict_e2e_question_only",
        "entitylink": "strict_e2e_entitylink",
        "entity_link": "strict_e2e_entitylink",
        "fallback": "strict_e2e_fallback",
        "scorecut": "strict_e2e_scorecutoff",
        "score_cutoff": "strict_e2e_scorecutoff",
        "soft": "strict_e2e_hop_aware",
        "soft_rerank": "strict_e2e_hop_aware",
        "rerank": "strict_e2e_hop_aware",
        "hop_aware": "strict_e2e_hop_aware",
        "optimized": "strict_e2e_hop_aware",
        "optimized_x": "strict_e2e_hop_aware",
        "source_filtered": "sourcefiltered",
    }
    key = aliases.get(key, key)
    if key not in BIOHOPR_BENCHMARK_POLICIES:
        known = ", ".join(sorted(BIOHOPR_BENCHMARK_POLICIES))
        raise ValueError(f"Unknown BioHopR benchmark policy {name!r}. Known policies: {known}")
    return dict(BIOHOPR_BENCHMARK_POLICIES[key])


def biohopr_path_options_from_policy(
    policy: dict[str, Any] | str = "balanced",
    *,
    split: str = "1hop_single",
) -> PathAnswerOptions:
    values = biohopr_benchmark_policy(policy) if isinstance(policy, str) else dict(policy)
    is_1hop = str(split).startswith("1hop")
    return PathAnswerOptions(
        max_answers=int(values.get("max_answers", 120)),
        kg_answer_strength=str(values.get("kg_answer_strength", "all")),
        path_co_mentioned_policy=str(values.get("path_co_mentioned_policy", "allow")),
        co_mentioned_fallback_min_score=float(values.get("co_mentioned_fallback_min_score", 0.2)),
        selection_mode=str(values.get("answer_selection_mode", "top_k")),
        min_path_score=values.get("min_path_score"),
        score_drop=float(values.get("score_drop", 1.35)),
        min_answers=int(values.get("min_answers", 3)),
        source_policy=str(values.get("biohopr_source_policy", "none")),
        weak_source_rescue_max=int(
            values.get("weak_source_rescue_max_1hop" if is_1hop else "weak_source_rescue_max_2hop", 0)
            or 0
        ),
        weak_source_rescue_min_score=values.get(
            "weak_source_rescue_min_score_1hop" if is_1hop else "weak_source_rescue_min_score_2hop"
        ),
    )


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def record_id(row: dict[str, Any]) -> int:
    value = row.get("id") if row.get("id") is not None else row.get("sample_id")
    return int(value)


def norm(value: Any) -> str:
    return " ".join(str(value or "").strip().casefold().split())


def split_csv(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        out: list[str] = []
        for item in value:
            out.extend(split_csv(item))
        return out
    text = str(value).strip()
    if not text:
        return []
    return [part.strip() for part in text.split(",") if part.strip()]


def _dedupe_keep_order(values: list[Any], *, max_items: int | None = None) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value or "").strip()
        if not text:
            continue
        key = norm(text)
        if key in seen:
            continue
        out.append(text)
        seen.add(key)
        if max_items is not None and len(out) >= max_items:
            break
    return out


def _dedupe_surface_keep_order(values: list[Any], *, max_items: int | None = None) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value or "").strip()
        if not text:
            continue
        key = " ".join(text.split())
        if key in seen:
            continue
        out.append(text)
        seen.add(key)
        if max_items is not None and len(out) >= max_items:
            break
    return out


def _compact_entity_surface(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip().strip("\"'`")).strip()


def _strip_common_chemical_suffix(text: str) -> str:
    suffixes = [
        " hydrochloride",
        " citrate",
        " sodium",
        " potassium",
        " sulfate",
        " sulphate",
        " phosphate",
        " acetate",
        " maleate",
        " mesylate",
        " tartrate",
        " hydrobromide",
    ]
    lowered = text.casefold()
    for suffix in suffixes:
        if lowered.endswith(suffix):
            return text[: -len(suffix)].strip()
    return text


def biohopr_entity_candidates(
    sources: list[str] | str,
    source_type: str | None = None,
    *,
    mode: str = "off",
    max_candidates: int = 1,
) -> list[str]:
    """Generate conservative source-name candidates without reading gold metadata."""
    source_list = [sources] if isinstance(sources, str) else list(sources or [])
    seeds = [_compact_entity_surface(str(item)) for item in source_list if str(item or "").strip()]
    if not seeds:
        return []
    if mode in {"", "off", "none"} or max_candidates <= 1:
        return _dedupe_keep_order(seeds, max_items=max(1, max_candidates))

    stype = norm(source_type).replace(" ", "_")
    candidates: list[str] = []
    for seed in seeds:
        candidates.append(seed)
        compact = _compact_entity_surface(seed)
        candidates.append(compact)
        if compact != compact.strip("."):
            candidates.append(compact.strip("."))

        if stype in {"gene", "protein", "gene/protein"}:
            candidates.append(compact.upper())
            candidates.append(compact.replace(" ", "-"))
            candidates.append(compact.replace("-", " "))
            no_sep = re.sub(r"[-\s]+", "", compact)
            if len(no_sep) >= 2:
                candidates.append(no_sep)
            interleukin = re.fullmatch(r"interleukin[-\s]*(\d+[a-z]?)", compact, flags=re.I)
            if interleukin:
                number = interleukin.group(1).upper()
                candidates.extend([f"IL{number}", f"IL-{number}", f"IL {number}"])
            il = re.fullmatch(r"IL[-\s]*(\d+[A-Za-z]?)", compact, flags=re.I)
            if il:
                number = il.group(1)
                candidates.extend([f"IL{number}", f"IL-{number}", f"interleukin {number}"])
            tnf = re.match(r"^(TNF)[-\s]*(alpha|beta|a|b)?$", compact, flags=re.I)
            if tnf:
                candidates.extend(["TNF", "TNF alpha", "TNF-alpha", "tumor necrosis factor"])
            if "-" in compact and re.fullmatch(r"[A-Za-z0-9-]+", compact):
                prefix = compact.split("-", 1)[0]
                if len(prefix) >= 2:
                    candidates.append(prefix)

        elif stype in {"chemical", "drug"}:
            candidates.extend([compact.lower(), compact.upper(), compact.title()])
            stripped = _strip_common_chemical_suffix(compact)
            if stripped != compact:
                candidates.append(stripped)

        elif stype in {"disease", "phenotype", "effect/phenotype"}:
            candidates.extend([compact.lower(), compact.title()])
            if compact.casefold().endswith(" disease"):
                candidates.append(compact[:-8].strip())
            if compact.casefold().endswith(" syndrome"):
                candidates.append(compact[:-9].strip())

    return _dedupe_surface_keep_order(candidates, max_items=max(1, max_candidates))


def biohopr_type_candidates(entity_type: str | None, *, allow_fallback: bool = False) -> str | list[str] | None:
    if not entity_type:
        return None
    canonical = str(entity_type).strip()
    if not allow_fallback:
        return canonical
    return _dedupe_keep_order(TYPE_FALLBACKS.get(canonical, [canonical]))


def _accepted_type_set(entity_type: str | None, *, allow_fallback: bool = False) -> set[str]:
    candidates = biohopr_type_candidates(entity_type, allow_fallback=allow_fallback)
    if candidates is None:
        return set()
    values = candidates if isinstance(candidates, list) else [candidates]
    return {norm(value) for value in values if norm(value)}


def _parse_list_like(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item).strip() for item in value if str(item or "").strip()]
    text = str(value).strip()
    if not text:
        return []
    if text.startswith("[") and text.endswith("]"):
        try:
            parsed = json.loads(text.replace("'", '"'))
            if isinstance(parsed, list):
                return [str(item).strip() for item in parsed if str(item or "").strip()]
        except Exception:
            pass
    return [part.strip(" []\"'") for part in text.split(",") if part.strip(" []\"'")]


def normalize_info_source_token(value: Any) -> str:
    token = str(value or "").strip().strip("[]\"'")
    if not token:
        return ""
    upper = token.upper()
    if upper.startswith(("PMID:", "PMCID:")):
        return upper.split(":", 1)[0]
    if ":" in upper:
        return upper.split(":", 1)[0]
    return re.sub(r"[^A-Z0-9_+-]+", "_", upper).strip("_")


def row_info_sources(row: pd.Series, hops: int) -> set[str]:
    sources: set[str] = set()
    for idx in range(1, 2 * hops, 2):
        col = f"path.{idx}.info_source"
        if col not in row.index:
            continue
        for raw in _parse_list_like(row.get(col)):
            token = normalize_info_source_token(raw)
            if token:
                sources.add(token)
    return sources


def path_info_sources(path: dict[str, Any]) -> set[str]:
    sources: set[str] = set()
    for key, value in path.items():
        if not key.endswith(".info_source"):
            continue
        for raw in _parse_list_like(value):
            token = normalize_info_source_token(raw)
            if token:
                sources.add(token)
    return sources


def has_curated_info_source(sources: set[str]) -> bool:
    return bool(sources & CURATED_INFO_SOURCE_PREFIXES)


def has_low_signal_info_source(sources: set[str]) -> bool:
    return bool(sources & LOW_SIGNAL_INFO_SOURCE_PREFIXES)


def is_generic_answer_label(value: Any) -> bool:
    answer = norm(value)
    return bool(answer and answer in GENERIC_DISEASE_TERMS)


def has_value(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, (list, tuple, set, dict)):
        return bool(value)
    try:
        return not bool(pd.isna(value))
    except Exception:
        return True


def cache_key(**kwargs: Any) -> str:
    return json.dumps(kwargs, ensure_ascii=False, sort_keys=True, default=str)


def relation_columns(hops: int) -> list[str]:
    return [f"path.{idx}.relation" for idx in range(1, 2 * hops, 2)]


def endpoint_name_col(hops: int) -> str:
    return f"path.{2 * hops}.name"


def endpoint_type_col(hops: int) -> str:
    return f"path.{2 * hops}.type"


def load_gold_hop_lookup(gold_dir: Path, split: str) -> dict[int, dict[str, Any]]:
    path = gold_dir / GOLD_FILES[split]
    if not path.exists():
        return {}
    rows = read_json(path)
    lookup: dict[int, dict[str, Any]] = {}
    for row in rows:
        try:
            lookup[int(row["id"])] = row
        except Exception:
            continue
    return lookup


def gold_hop_contract_fields(split: str, gold: dict[str, Any] | None) -> dict[str, Any]:
    """Map BioHopR official hop metadata to query fields without gold answers."""
    if not gold:
        return {"structured_query_source": "gold_hops", "gold_hop_missing": True}
    common = {
        "target_type": gold.get("target_type"),
        "relation_hop1": gold.get("relation_hop1"),
        "relation_hop2": gold.get("relation_hop2"),
        "structured_query_source": "gold_hops",
        "gold_hop_missing": False,
    }
    if split.startswith("1hop"):
        return {
            **common,
            "query_entity": gold.get("hop1"),
            "query_type": gold.get("hop1_type"),
            "hops": 1,
        }
    return {
        **common,
        "query_entity": gold.get("hop2"),
        "query_type": gold.get("hop2_type"),
        "bridge_entity": gold.get("hop1"),
        "bridge_type": gold.get("hop1_type"),
        "hops": 2,
    }


def enrich_records_with_gold_hops(
    split: str,
    records: list[dict[str, Any]],
    gold_dir: Path,
) -> list[dict[str, Any]]:
    gold_by_id = load_gold_hop_lookup(gold_dir, split)
    enriched: list[dict[str, Any]] = []
    for row in records:
        gold = gold_by_id.get(record_id(row))
        enriched.append({**row, **gold_hop_contract_fields(split, gold)})
    return enriched


def _cache_context(lock: Any):
    return lock if lock is not None else nullcontext()


def query_kg_cached(
    *,
    source: list[str] | str,
    source_type: str | None,
    target_type: str | None,
    hops: int,
    target: list[str] | str | None = None,
    directed: bool = False,
    limit: int | None = None,
    relation_info_source: list[str] | str | None = None,
    cache: MutableMapping[str, pd.DataFrame] | None = None,
    cache_lock: Any = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    from querygraph import query_relation

    source_list = [source] if isinstance(source, str) else list(source or [])
    source_list = [str(item).strip() for item in source_list if str(item or "").strip()]
    if target is None:
        target_list = None
    else:
        target_list = [target] if isinstance(target, str) else list(target or [])
        target_list = [str(item).strip() for item in target_list if str(item or "").strip()]
    if not source_list:
        return pd.DataFrame(), {"cache_hit": False, "error": "missing_source", "elapsed": 0.0}

    source_type_values = source_type if isinstance(source_type, (list, tuple, set)) else [source_type]
    target_type_values = target_type if isinstance(target_type, (list, tuple, set)) else [target_type]
    source_type_values = list(source_type_values) or [None]
    target_type_values = list(target_type_values) or [None]
    if len(source_type_values) == 1:
        source_type = source_type_values[0]
    if len(target_type_values) == 1:
        target_type = target_type_values[0]
    if len(source_type_values) > 1 or len(target_type_values) > 1:
        frames: list[pd.DataFrame] = []
        metas: list[dict[str, Any]] = []
        for source_type_value in source_type_values:
            for target_type_value in target_type_values:
                df_part, meta_part = query_kg_cached(
                    source=source_list,
                    source_type=source_type_value,
                    target_type=target_type_value,
                    hops=hops,
                    target=target_list,
                    directed=directed,
                    limit=limit,
                    relation_info_source=relation_info_source,
                    cache=cache,
                    cache_lock=cache_lock,
                )
                if not df_part.empty:
                    frames.append(df_part)
                metas.append(meta_part)
        elapsed = round(sum(float(meta.get("elapsed") or 0.0) for meta in metas), 2)
        errors = [str(meta.get("error")) for meta in metas if meta.get("error")]
        combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        return combined, {
            "cache_hit": bool(metas and all(meta.get("cache_hit") for meta in metas)),
            "error": "; ".join(_dedupe_keep_order(errors)),
            "elapsed": elapsed,
            "fanout_calls": len(metas),
        }

    key = cache_key(
        source=source_list,
        source_type=source_type,
        target=target_list,
        target_type=target_type,
        hops=hops,
        directed=directed,
        limit=limit,
        relation_info_source=split_csv(relation_info_source),
    )
    if cache is not None:
        with _cache_context(cache_lock):
            cached = cache.get(key)
        if cached is not None:
            return cached.copy(), {"cache_hit": True, "error": "", "elapsed": 0.0}

    start = time.time()
    try:
        df = query_relation(
            source_entity_set=source_list,
            target_entity_set=target_list,
            source_entity_type=source_type,
            target_entity_type=target_type,
            multi_hop=hops,
            directed=directed,
            relation_info_source=relation_info_source,
            limit_number=limit,
            return_type="dataframe",
        )
        if cache is not None:
            with _cache_context(cache_lock):
                cache[key] = df.copy()
        return df, {"cache_hit": False, "error": "", "elapsed": round(time.time() - start, 2)}
    except Exception as exc:
        return pd.DataFrame(), {"cache_hit": False, "error": repr(exc), "elapsed": round(time.time() - start, 2)}


def relation_priority(relation: Any) -> float:
    rel = norm(relation).replace(" ", "_")
    if rel in STRONG_RELATIONS:
        return 3.0
    if rel in WEAK_RELATIONS:
        return 0.2
    return 1.0


def numeric_value(value: Any, default: float = 0.0) -> float:
    try:
        if value is None or pd.isna(value):
            return default
    except Exception:
        pass
    try:
        return float(value)
    except Exception:
        return default


def literature_support_score(row: pd.Series, hops: int) -> float:
    """Estimate how useful a weak literature trace is for BioHopR answer sets."""
    support = 0.0
    for col in relation_columns(hops):
        rel = norm(row.get(col, "")).replace(" ", "_")
        if rel not in {"co_mentioned", "mentioned"}:
            continue
        info_col = col.replace(".relation", ".info_source_length")
        count = numeric_value(row.get(info_col), 0.0)
        if count > 0:
            support += min(math.log1p(count) * 0.65, 2.2)

        edge_idx = col.split(".")[1]
        snippets = " ".join(
            norm(row.get(f"path.{edge_idx}.{field}", ""))
            for field in ["original_relation", "condition", "info_source"]
            if f"path.{edge_idx}.{field}" in row.index
        )
        if snippets:
            answer = norm(row.get(endpoint_name_col(hops), ""))
            source = norm(row.get("path.0.name", ""))
            if answer and answer in snippets:
                support += 0.35
            if source and source in snippets:
                support += 0.35
    return min(support, 3.0)


def path_template_bonus(row: pd.Series, hops: int, contract: Any) -> float:
    policy = dict(contract.evidence_policy or {})
    if not policy.get("path_template_rerank", True):
        return 0.0
    target_type = norm(policy.get("target_type") or contract.expected_answer_type)
    intermediate_type = norm(policy.get("intermediate_type"))
    bonus = 0.0
    for col in relation_columns(hops):
        rel = norm(row.get(col, "")).replace(" ", "_")
        if target_type == "disease" and rel in {"disease_relation", "causally_related_to", "is_a"}:
            bonus += 0.45
        elif target_type == "chemical" and rel in {"pharmaco_effect_on", "chemical_relation", "chemicals_relation"}:
            bonus += 0.45
        elif target_type in {"gene", "protein"} and rel in {"regulated", "interact_with", "functionally_related_to"}:
            bonus += 0.35
        elif intermediate_type == "gene" and rel in {"regulated", "interact_with", "functionally_related_to"}:
            bonus += 0.25
    return min(bonus, 1.2)


def generic_endpoint_penalty(row: pd.Series, hops: int, contract: Any) -> float:
    policy = dict(contract.evidence_policy or {})
    if not policy.get("generic_disease_penalty", True):
        return 0.0
    answer = norm(row.get(endpoint_name_col(hops), ""))
    if not answer:
        return 0.0
    if answer in GENERIC_DISEASE_TERMS:
        return 1.25
    if len(answer.split()) == 1 and any(term in answer for term in {"cancer", "disease", "syndrome"}):
        return 0.7
    return 0.0


def hub_like_penalty(row: pd.Series, hops: int, contract: Any) -> float:
    policy = dict(contract.evidence_policy or {})
    if not policy.get("hub_penalty", True):
        return 0.0
    penalty = 0.0
    for idx in range(0, 2 * hops + 1, 2):
        name = norm(row.get(f"path.{idx}.name", ""))
        if any(name.startswith(prefix) for prefix in HUB_LIKE_PREFIXES):
            penalty += 0.25
    return min(penalty, 0.8)


def row_has_strong_path(row: pd.Series, hops: int) -> bool:
    for col in relation_columns(hops):
        rel = norm(row.get(col, "")).replace(" ", "_")
        if rel in STRONG_RELATIONS:
            return True
    return False


def row_has_weak_edge(row: pd.Series, hops: int) -> bool:
    for col in relation_columns(hops):
        rel = norm(row.get(col, "")).replace(" ", "_")
        if rel in WEAK_RELATIONS:
            return True
    return False


def row_score(row: pd.Series, hops: int, contract: Any | None = None) -> float:
    score = 0.0
    has_weak = row_has_weak_edge(row, hops)
    has_strong = row_has_strong_path(row, hops)
    for col in relation_columns(hops):
        score += relation_priority(row.get(col, ""))
        info_col = col.replace(".relation", ".info_source_length")
        try:
            score += min(math.log1p(float(row.get(info_col) or 0)), 3.0) * 0.05
        except Exception:
            pass
    if has_weak:
        score += literature_support_score(row, hops)
    if contract is not None:
        score += path_template_bonus(row, hops, contract)
        score -= generic_endpoint_penalty(row, hops, contract)
        score -= hub_like_penalty(row, hops, contract)
    if has_weak and not has_strong:
        score -= 0.85
    elif has_weak:
        score -= 0.35
    return score


def trim_ranked_answers(
    ranked: list[dict[str, Any]],
    options: PathAnswerOptions,
) -> list[dict[str, Any]]:
    if options.source_policy == "curated_with_weak_rescue" and options.weak_source_rescue_max >= 0:
        filtered: list[dict[str, Any]] = []
        weak_kept = 0
        for item in ranked:
            if item.get("source_role") == "weak_literature_rescue":
                if weak_kept >= options.weak_source_rescue_max:
                    continue
                weak_kept += 1
            filtered.append(item)
        ranked = filtered
    if options.selection_mode != "score_cutoff" or not ranked:
        return ranked[: options.max_answers]
    best = float(ranked[0].get("score") or 0.0)
    threshold = best - options.score_drop
    if options.min_path_score is not None:
        threshold = max(threshold, options.min_path_score)

    kept = [item for item in ranked if float(item.get("score") or 0.0) >= threshold]
    if len(kept) < options.min_answers:
        kept = ranked[: options.min_answers]
    return kept[: options.max_answers]


def _best_bridge_validation_row(df: pd.DataFrame) -> pd.Series | None:
    if df.empty:
        return None
    scored = [(row_score(row, 1), idx) for idx, row in df.iterrows()]
    _, best_idx = max(scored, key=lambda item: item[0])
    return df.loc[best_idx]


def _edge_value(row: pd.Series | None, edge_idx: int, field: str) -> Any:
    if row is None:
        return ""
    return row.get(f"path.{edge_idx}.{field}", "")


def compose_structured_bridge_rows(
    *,
    contract: Any,
    bridge_validation_df: pd.DataFrame,
    answer_df: pd.DataFrame,
) -> pd.DataFrame:
    if answer_df.empty:
        return pd.DataFrame()
    policy = dict(contract.evidence_policy or {})
    metadata = dict(contract.metadata or {})
    source = (contract.source_mentions or [""])[0]
    source_type = policy.get("source_type")
    bridge = str(metadata.get("bridge_entity") or "").strip()
    bridge_type = policy.get("intermediate_type") or metadata.get("bridge_type")
    bridge_row = _best_bridge_validation_row(bridge_validation_df)
    bridge_validated = bridge_row is not None
    composed: list[dict[str, Any]] = []

    for _, answer_row in answer_df.iterrows():
        out: dict[str, Any] = {
            "path.0.name": source,
            "path.0.type": source_type,
            "path.2.name": bridge,
            "path.2.type": bridge_type,
            "path.4.name": answer_row.get("path.2.name", ""),
            "path.4.type": answer_row.get("path.2.type", ""),
            "query.structured_bridge": True,
            "query.bridge_validated": bridge_validated,
        }
        if bridge_validated:
            for field in [
                "relation",
                "info_source",
                "info_source_length",
                "original_relation",
                "condition",
                "relation_confidence",
            ]:
                out[f"path.1.{field}"] = _edge_value(bridge_row, 1, field)
        else:
            out["path.1.relation"] = "provided_bridge"
            out["path.1.original_relation"] = "BioHopR official hop bridge supplied; KG validation returned no row."
        for field in [
            "relation",
            "info_source",
            "info_source_length",
            "original_relation",
            "condition",
            "relation_confidence",
        ]:
            out[f"path.3.{field}"] = answer_row.get(f"path.1.{field}", "")
        composed.append(out)
    return pd.DataFrame(composed)


def query_structured_bridge_path(
    contract: Any,
    *,
    kg_limit: int = 400,
    relation_info_source: list[str] | str | None = None,
    cache: MutableMapping[str, pd.DataFrame] | None = None,
    cache_lock: Any = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    policy = dict(contract.evidence_policy or {})
    metadata = dict(contract.metadata or {})
    source = list(contract.source_mentions or [])
    source_type = policy.get("source_type")
    bridge = str(metadata.get("bridge_entity") or "").strip()
    bridge_type = policy.get("intermediate_type")
    target_type = policy.get("target_type") or contract.expected_answer_type
    start = time.time()

    if not bridge or not bridge_type:
        df, meta = query_kg_cached(
            source=source,
            source_type=source_type,
            target_type=target_type,
            hops=2,
            directed=False,
            limit=kg_limit,
            relation_info_source=relation_info_source,
            cache=cache,
            cache_lock=cache_lock,
        )
        meta.update({"mode": "gold_hops_missing_bridge_fallback_2hop", "calls": 1})
        return df, meta

    bridge_df, bridge_meta = query_kg_cached(
        source=source,
        source_type=source_type,
        target=[bridge],
        target_type=bridge_type,
        hops=1,
        directed=False,
        limit=kg_limit,
        relation_info_source=relation_info_source,
        cache=cache,
        cache_lock=cache_lock,
    )
    answer_df, answer_meta = query_kg_cached(
        source=[bridge],
        source_type=bridge_type,
        target_type=target_type,
        hops=1,
        directed=False,
        limit=kg_limit,
        relation_info_source=relation_info_source,
        cache=cache,
        cache_lock=cache_lock,
    )
    df = compose_structured_bridge_rows(
        contract=contract,
        bridge_validation_df=bridge_df,
        answer_df=answer_df,
    )
    errors = [meta.get("error", "") for meta in (bridge_meta, answer_meta) if meta.get("error")]
    return df, {
        "cache_hit": bool(bridge_meta.get("cache_hit") and answer_meta.get("cache_hit")),
        "error": "; ".join(errors),
        "elapsed": round(time.time() - start, 2),
        "mode": "gold_hops_bridge",
        "calls": 2,
        "bridge_entity": bridge,
        "bridge_type": bridge_type,
        "bridge_validation_rows": int(len(bridge_df)),
        "bridge_answer_rows": int(len(answer_df)),
        "bridge_validation_error": bridge_meta.get("error", ""),
        "bridge_answer_error": answer_meta.get("error", ""),
    }


def _filter_endpoint_type(
    df: pd.DataFrame,
    *,
    hops: int,
    entity_type: str | None,
    allow_type_fallback: bool,
) -> pd.DataFrame:
    if df.empty or not entity_type:
        return df
    col = endpoint_type_col(hops)
    if col not in df.columns:
        return df
    accepted = _accepted_type_set(entity_type, allow_fallback=allow_type_fallback)
    return df[df[col].fillna("").astype(str).map(norm).isin(accepted)]


def _rank_bridge_names(
    bridge_df: pd.DataFrame,
    *,
    max_bridge_entities: int,
) -> list[str]:
    if bridge_df.empty or "path.2.name" not in bridge_df.columns:
        return []
    scored: dict[str, tuple[float, str]] = {}
    for _, row in bridge_df.iterrows():
        bridge = str(row.get("path.2.name") or "").strip()
        if not bridge:
            continue
        key = norm(bridge)
        score = row_score(row, 1, None)
        current = scored.get(key)
        if current is None or score > current[0]:
            scored[key] = (score, bridge)
    ranked = sorted(scored.values(), key=lambda item: (-item[0], norm(item[1])))
    return [bridge for _score, bridge in ranked[: max(1, max_bridge_entities)]]


def compose_discovered_bridge_rows(
    *,
    contract: Any,
    bridge_df: pd.DataFrame,
    answer_df: pd.DataFrame,
) -> pd.DataFrame:
    if bridge_df.empty or answer_df.empty:
        return pd.DataFrame()
    policy = dict(contract.evidence_policy or {})
    source = (contract.source_mentions or [""])[0]
    source_type = policy.get("source_type")
    bridge_type = policy.get("intermediate_type")
    best_bridge_rows: dict[str, tuple[float, pd.Series]] = {}
    for _, row in bridge_df.iterrows():
        bridge = str(row.get("path.2.name") or "").strip()
        if not bridge:
            continue
        key = norm(bridge)
        score = row_score(row, 1, None)
        current = best_bridge_rows.get(key)
        if current is None or score > current[0]:
            best_bridge_rows[key] = (score, row)

    composed: list[dict[str, Any]] = []
    for _, answer_row in answer_df.iterrows():
        bridge = str(answer_row.get("path.0.name") or "").strip()
        bridge_key = norm(bridge)
        bridge_match = best_bridge_rows.get(bridge_key)
        if bridge_match is None:
            continue
        bridge_row = bridge_match[1]
        out: dict[str, Any] = {
            "path.0.name": source,
            "path.0.type": source_type,
            "path.2.name": bridge,
            "path.2.type": answer_row.get("path.0.type", "") or bridge_row.get("path.2.type", "") or bridge_type,
            "path.4.name": answer_row.get("path.2.name", ""),
            "path.4.type": answer_row.get("path.2.type", ""),
            "query.discovered_bridge": True,
        }
        for field in [
            "relation",
            "info_source",
            "info_source_length",
            "original_relation",
            "condition",
            "relation_confidence",
        ]:
            out[f"path.1.{field}"] = bridge_row.get(f"path.1.{field}", "")
            out[f"path.3.{field}"] = answer_row.get(f"path.1.{field}", "")
        composed.append(out)
    return pd.DataFrame(composed)


def query_discovered_bridge_path(
    contract: Any,
    *,
    source: list[str],
    source_type: str | None,
    target_type: str | None,
    kg_limit: int,
    relation_info_source: list[str] | str | None,
    cache: MutableMapping[str, pd.DataFrame] | None,
    cache_lock: Any,
    source_candidates: list[str],
    allow_type_fallback: bool,
    max_bridge_entities: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    policy = dict(contract.evidence_policy or {})
    intermediate_type = policy.get("intermediate_type")
    if not intermediate_type:
        return pd.DataFrame(), {"error": "missing_intermediate_type", "elapsed": 0.0, "calls": 0}
    start = time.time()
    bridge_df, bridge_meta = query_kg_cached(
        source=source_candidates or source,
        source_type=biohopr_type_candidates(source_type, allow_fallback=allow_type_fallback),
        target_type=biohopr_type_candidates(intermediate_type, allow_fallback=allow_type_fallback),
        hops=1,
        directed=False,
        limit=kg_limit,
        relation_info_source=relation_info_source,
        cache=cache,
        cache_lock=cache_lock,
    )
    bridge_filtered = _filter_endpoint_type(
        bridge_df,
        hops=1,
        entity_type=intermediate_type,
        allow_type_fallback=allow_type_fallback,
    )
    bridge_names = _rank_bridge_names(
        bridge_filtered,
        max_bridge_entities=max_bridge_entities,
    )
    if not bridge_names:
        return pd.DataFrame(), {
            "error": bridge_meta.get("error", ""),
            "elapsed": round(time.time() - start, 2),
            "calls": 1,
            "bridge_probe_rows": int(len(bridge_df)),
            "bridge_probe_filtered_rows": int(len(bridge_filtered)),
            "bridge_candidates": [],
        }
    answer_df, answer_meta = query_kg_cached(
        source=bridge_names,
        source_type=biohopr_type_candidates(intermediate_type, allow_fallback=allow_type_fallback),
        target_type=biohopr_type_candidates(target_type, allow_fallback=allow_type_fallback),
        hops=1,
        directed=False,
        limit=kg_limit,
        relation_info_source=relation_info_source,
        cache=cache,
        cache_lock=cache_lock,
    )
    answer_filtered = _filter_endpoint_type(
        answer_df,
        hops=1,
        entity_type=target_type,
        allow_type_fallback=allow_type_fallback,
    )
    composed = compose_discovered_bridge_rows(
        contract=contract,
        bridge_df=bridge_filtered,
        answer_df=answer_filtered,
    )
    errors = [meta.get("error", "") for meta in (bridge_meta, answer_meta) if meta.get("error")]
    return composed, {
        "error": "; ".join(_dedupe_keep_order(errors)),
        "elapsed": round(time.time() - start, 2),
        "calls": 2,
        "bridge_probe_rows": int(len(bridge_df)),
        "bridge_probe_filtered_rows": int(len(bridge_filtered)),
        "bridge_answer_rows": int(len(answer_df)),
        "bridge_answer_filtered_rows": int(len(answer_filtered)),
        "bridge_candidates": bridge_names[:10],
    }


def query_path_rows(
    contract: Any,
    *,
    query_source: str = "question",
    kg_limit: int = 400,
    relation_info_source: list[str] | str | None = None,
    cache: MutableMapping[str, pd.DataFrame] | None = None,
    cache_lock: Any = None,
    entity_linking: str = "off",
    max_query_entities: int = 1,
    type_fallback: bool = False,
    discover_bridge_fallback: bool = False,
    max_bridge_entities: int = 25,
    relax_source_type_fallback: bool = False,
    relax_target_type_fallback: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    policy = dict(contract.evidence_policy or {})
    source = list(contract.source_mentions or [])
    source_type = policy.get("source_type")
    target_type = policy.get("target_type") or contract.expected_answer_type
    hops = int(policy.get("max_hops") or 1)
    if query_source == "gold_hops" and hops == 2:
        return query_structured_bridge_path(
            contract,
            kg_limit=kg_limit,
            relation_info_source=relation_info_source,
            cache=cache,
            cache_lock=cache_lock,
        )
    source_candidates = biohopr_entity_candidates(
        source,
        source_type,
        mode=entity_linking,
        max_candidates=max_query_entities,
    )
    strict_source = _dedupe_keep_order(source, max_items=max(1, max_query_entities))
    query_mode = "extracted_contract_path" if query_source == "extracted" else "question_contract_path"
    stages: list[dict[str, Any]] = [
        {
            "stage": "initial",
            "source": strict_source,
            "source_type": source_type,
            "target_type": target_type,
            "allow_type_fallback": False,
        }
    ]
    if source_candidates and source_candidates != strict_source:
        stages.append(
            {
                "stage": "entity_surface",
                "source": source_candidates,
                "source_type": source_type,
                "target_type": target_type,
                "allow_type_fallback": False,
            }
        )
    if type_fallback:
        stages.append(
            {
                "stage": "type_expand",
                "source": source_candidates or strict_source,
                "source_type": biohopr_type_candidates(source_type, allow_fallback=True),
                "target_type": biohopr_type_candidates(target_type, allow_fallback=True),
                "allow_type_fallback": True,
            }
        )
    if relax_source_type_fallback:
        stages.append(
            {
                "stage": "relax_source_type",
                "source": source_candidates or strict_source,
                "source_type": None,
                "target_type": biohopr_type_candidates(target_type, allow_fallback=type_fallback),
                "allow_type_fallback": type_fallback,
            }
        )
    if relax_target_type_fallback:
        stages.append(
            {
                "stage": "relax_target_type",
                "source": source_candidates or strict_source,
                "source_type": biohopr_type_candidates(source_type, allow_fallback=type_fallback),
                "target_type": None,
                "allow_type_fallback": type_fallback,
            }
        )

    start = time.time()
    selected_df = pd.DataFrame()
    selected_stage = stages[0]["stage"]
    selected_meta: dict[str, Any] = {}
    first_nonempty_df: pd.DataFrame | None = None
    first_nonempty_stage = ""
    stage_results: list[dict[str, Any]] = []
    errors: list[str] = []

    for stage in stages:
        df_stage, meta_stage = query_kg_cached(
            source=stage["source"],
            source_type=stage["source_type"],
            target_type=stage["target_type"],
            hops=hops,
            directed=False,
            limit=kg_limit,
            relation_info_source=relation_info_source,
            cache=cache,
            cache_lock=cache_lock,
        )
        if not df_stage.empty:
            df_stage = df_stage.copy()
            df_stage["query.stage"] = stage["stage"]
        filtered_stage = filter_path_rows(
            df_stage,
            contract,
            allow_type_fallback=bool(stage.get("allow_type_fallback")),
        )
        stage_result = {
            "stage": stage["stage"],
            "rows": int(len(df_stage)),
            "filtered_rows": int(len(filtered_stage)),
            "source_count": len(stage["source"]),
            "source_type": stage["source_type"],
            "target_type": stage["target_type"],
            "allow_type_fallback": bool(stage.get("allow_type_fallback")),
            "cache_hit": bool(meta_stage.get("cache_hit")),
            "error": meta_stage.get("error", ""),
            "elapsed": meta_stage.get("elapsed", 0.0),
        }
        stage_results.append(stage_result)
        if meta_stage.get("error"):
            errors.append(str(meta_stage.get("error")))
        if first_nonempty_df is None and not df_stage.empty:
            first_nonempty_df = df_stage
            first_nonempty_stage = stage["stage"]
        if not filtered_stage.empty:
            selected_df = df_stage
            selected_stage = stage["stage"]
            selected_meta = meta_stage
            break

    if selected_df.empty and discover_bridge_fallback and hops == 2:
        bridge_df, bridge_meta = query_discovered_bridge_path(
            contract,
            source=source,
            source_type=source_type,
            target_type=target_type,
            kg_limit=kg_limit,
            relation_info_source=relation_info_source,
            cache=cache,
            cache_lock=cache_lock,
            source_candidates=source_candidates or strict_source,
            allow_type_fallback=type_fallback,
            max_bridge_entities=max_bridge_entities,
        )
        filtered_bridge_df = filter_path_rows(
            bridge_df,
            contract,
            allow_type_fallback=type_fallback,
        )
        bridge_stage_result = {
            "stage": "discovered_bridge",
            "rows": int(len(bridge_df)),
            "filtered_rows": int(len(filtered_bridge_df)),
            "source_count": len(source_candidates or strict_source),
            "source_type": biohopr_type_candidates(source_type, allow_fallback=type_fallback),
            "target_type": biohopr_type_candidates(target_type, allow_fallback=type_fallback),
            "allow_type_fallback": bool(type_fallback),
            "cache_hit": bool(bridge_meta.get("cache_hit")),
            "error": bridge_meta.get("error", ""),
            "elapsed": bridge_meta.get("elapsed", 0.0),
            "bridge_probe_rows": bridge_meta.get("bridge_probe_rows", 0),
            "bridge_probe_filtered_rows": bridge_meta.get("bridge_probe_filtered_rows", 0),
            "bridge_answer_rows": bridge_meta.get("bridge_answer_rows", 0),
            "bridge_answer_filtered_rows": bridge_meta.get("bridge_answer_filtered_rows", 0),
            "bridge_candidates": bridge_meta.get("bridge_candidates", []),
        }
        stage_results.append(bridge_stage_result)
        if bridge_meta.get("error"):
            errors.append(str(bridge_meta.get("error")))
        if not bridge_df.empty and first_nonempty_df is None:
            first_nonempty_df = bridge_df
            first_nonempty_stage = "discovered_bridge"
        if not filtered_bridge_df.empty:
            selected_df = bridge_df
            selected_stage = "discovered_bridge"
            selected_meta = bridge_meta

    if selected_df.empty and first_nonempty_df is not None:
        selected_df = first_nonempty_df
        selected_stage = first_nonempty_stage
    elif selected_df.empty:
        selected_meta = stage_results[-1] if stage_results else {}

    query_mode = "extracted_contract_path" if query_source == "extracted" else "question_contract_path"
    error_text = "; ".join(_dedupe_keep_order(errors)) if errors else str(selected_meta.get("error", "") or "")
    meta = {
        "cache_hit": bool(selected_meta.get("cache_hit", False)),
        "error": error_text,
        "elapsed": round(time.time() - start, 2),
        "mode": f"{query_mode}:{selected_stage}" if selected_stage != "initial" else query_mode,
        "calls": len(stage_results),
        "selected_stage": selected_stage,
        "fallback_level": selected_stage,
        "fallback_used_query_repair": selected_stage != "initial",
        "stage_results": stage_results,
        "source_candidates_preview": (source_candidates or strict_source)[:10],
        "source_candidate_count": len(source_candidates or strict_source),
        "entity_linking": entity_linking,
        "type_fallback": bool(type_fallback),
        "discover_bridge_fallback": bool(discover_bridge_fallback),
        "max_bridge_entities": int(max_bridge_entities),
        "relax_source_type_fallback": bool(relax_source_type_fallback),
        "relax_target_type_fallback": bool(relax_target_type_fallback),
    }
    return selected_df, meta


def filter_path_rows(
    df: pd.DataFrame,
    contract: Any,
    *,
    allow_type_fallback: bool = False,
) -> pd.DataFrame:
    if df.empty:
        return df
    policy = dict(contract.evidence_policy or {})
    hops = int(policy.get("max_hops") or 1)
    target_type = policy.get("target_type") or contract.expected_answer_type
    intermediate_type = policy.get("intermediate_type")
    out = df
    end_type = endpoint_type_col(hops)
    if target_type and end_type in out.columns:
        target_types = _accepted_type_set(target_type, allow_fallback=allow_type_fallback)
        out = out[out[end_type].fillna("").astype(str).map(norm).isin(target_types)]
    if hops == 2 and intermediate_type and "path.2.type" in out.columns:
        intermediate_types = _accepted_type_set(intermediate_type, allow_fallback=allow_type_fallback)
        out = out[out["path.2.type"].fillna("").astype(str).map(norm).isin(intermediate_types)]
    return out


def weak_path_has_context_support(row: pd.Series, hops: int, contract: Any) -> bool:
    source_terms = [norm(value) for value in (contract.source_mentions or []) if norm(value)]
    answer = norm(row.get(endpoint_name_col(hops), ""))
    texts: list[str] = []
    for idx in range(1, 2 * hops, 2):
        for field in ["original_relation", "condition", "info_source"]:
            col = f"path.{idx}.{field}"
            if col in row.index:
                texts.append(norm(row.get(col, "")))
    evidence_text = " ".join(texts)
    if evidence_text and answer and any(source in evidence_text for source in source_terms):
        if answer in evidence_text:
            return True
    return bool(answer and source_terms)


def row_source_role(
    row: pd.Series,
    hops: int,
    options: PathAnswerOptions,
    *,
    score: float,
    weak_rescue_used: int,
) -> tuple[bool, str]:
    """Apply BioHopR source policy after path scoring but before final answer use."""
    policy = options.source_policy
    if policy in {"", "none", "all"}:
        return True, "unfiltered"

    sources = row_info_sources(row, hops)
    has_curated = has_curated_info_source(sources)
    has_low_signal = has_low_signal_info_source(sources)
    has_weak = row_has_weak_edge(row, hops)
    has_strong = row_has_strong_path(row, hops)
    answer = row.get(endpoint_name_col(hops), "")

    if has_curated and has_strong:
        return True, "curated"

    if policy == "curated_only":
        return False, "dropped_non_curated"

    if policy == "curated_with_weak_rescue":
        if (
            has_weak
            and not has_strong
            and has_low_signal
            and score >= (
                options.weak_source_rescue_min_score
                if options.weak_source_rescue_min_score is not None
                else options.co_mentioned_fallback_min_score + options.weak_source_rescue_score_bonus
            )
            and not is_generic_answer_label(answer)
        ):
            return True, "weak_literature_rescue"
        return False, "dropped_non_curated"

    return True, "unfiltered"


def collect_path_answers(
    df: pd.DataFrame,
    contract: Any,
    options: PathAnswerOptions,
    *,
    allow_weak_fallback: bool = False,
) -> tuple[list[str], list[dict[str, Any]]]:
    if df.empty:
        return [], []
    hops = int((contract.evidence_policy or {}).get("max_hops") or 1)
    answer_col = endpoint_name_col(hops)
    if answer_col not in df.columns:
        return [], []
    scored: dict[str, dict[str, Any]] = {}
    for _, row in df.iterrows():
        has_strong = row_has_strong_path(row, hops)
        has_weak = row_has_weak_edge(row, hops)
        if options.kg_answer_strength == "strong_path" and not has_strong and not allow_weak_fallback:
            continue
        if has_weak and not has_strong and not allow_weak_fallback:
            if options.path_co_mentioned_policy in {"fallback_trace", "strong_first_fallback"}:
                continue
        if allow_weak_fallback:
            if not has_weak or has_strong:
                continue
            if not weak_path_has_context_support(row, hops, contract):
                continue
        answer = str(row.get(answer_col) or "").strip()
        if not answer:
            continue
        key = norm(answer)
        score = row_score(row, hops, contract)
        if allow_weak_fallback:
            score += 2.5
            if score < options.co_mentioned_fallback_min_score:
                continue
        source_allowed, source_role = row_source_role(
            row,
            hops,
            options,
            score=score,
            weak_rescue_used=0,
        )
        if not source_allowed:
            continue
        current = scored.get(key)
        if current is None or score > current["score"]:
            path = {
                "answer": answer,
                "score": round(score, 4),
                "source_policy": options.source_policy,
                "source_role": source_role,
                "info_source_tokens": sorted(row_info_sources(row, hops))[:20],
            }
            if allow_weak_fallback:
                path["fallback_used_co_mentioned"] = True
                path["evidence_strength"] = "weak_co_mentioned_fallback"
            for idx in range(0, 2 * hops + 1):
                for field in ["name", "type", "relation", "info_source", "original_relation"]:
                    col = f"path.{idx}.{field}"
                    if col in df.columns:
                        value = row.get(col)
                        if has_value(value):
                            path[col] = str(value)[:500]
            scored[key] = path
    ranked = sorted(scored.values(), key=lambda item: (-float(item["score"]), norm(item["answer"])))
    trimmed = trim_ranked_answers(ranked, options)
    return [item["answer"] for item in trimmed], trimmed


def extract_answer_set(
    df: pd.DataFrame,
    contract: Any,
    options: PathAnswerOptions | None = None,
) -> tuple[list[str], list[dict[str, Any]]]:
    options = options or PathAnswerOptions()
    answers, paths = collect_path_answers(df, contract, options, allow_weak_fallback=False)
    if answers or options.path_co_mentioned_policy != "strong_first_fallback":
        return answers, paths
    return collect_path_answers(df, contract, options, allow_weak_fallback=True)


def path_relation_summary(
    df: pd.DataFrame,
    hops: int,
    *,
    path_co_mentioned_policy: str = "strong_first_fallback",
) -> dict[str, Any]:
    counts: dict[str, int] = {}
    weak_rows = 0
    strong_rows = 0
    if df.empty:
        return {"rows": 0, "weak_rows": 0, "strong_path_rows": 0, "relations": {}}
    rel_cols = [col for col in relation_columns(hops) if col in df.columns]
    for _, row in df.iterrows():
        row_rels = [str(row.get(col) or "") for col in rel_cols]
        if row_rels and all(norm(rel).replace(" ", "_") in WEAK_RELATIONS for rel in row_rels):
            weak_rows += 1
        if any(norm(rel).replace(" ", "_") in STRONG_RELATIONS for rel in row_rels):
            strong_rows += 1
        for rel in row_rels:
            key = rel or ""
            counts[key] = counts.get(key, 0) + 1
    return {
        "rows": int(len(df)),
        "weak_rows": weak_rows,
        "strong_path_rows": strong_rows,
        "strong_path_found": strong_rows > 0,
        "fallback_used_co_mentioned": weak_rows > 0 and strong_rows == 0,
        "co_mentioned_main_path_policy": path_co_mentioned_policy,
        "relations": counts,
    }
