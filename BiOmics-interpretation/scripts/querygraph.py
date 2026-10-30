"""Agent_brick query module - HTTP API client for brick KG API service.

This module provides query_cypher and query_relation functions that call
the remote brick KG API service instead of direct Neo4j connection.
"""

import argparse
import json
import os
import re
import shlex
import sys
from typing import Callable, List, Optional, Union
import pandas as pd
import requests

try:
    from kgpolicy import (
        KG_POLICY_CHOICES,
        canonicalize_info_sources,
        expand_info_source_filter,
        load_kg_schema,
        resolve_relation_aliases,
        select_kg_policy,
    )
except Exception:  # pragma: no cover - keeps qr usable if policy module is unavailable
    KG_POLICY_CHOICES = ("manual", "auto", "source-scoped", "direct", "multihop", "off")
    canonicalize_info_sources = None
    expand_info_source_filter = None
    load_kg_schema = None
    resolve_relation_aliases = None
    select_kg_policy = None

# brick KG API configuration.
API_BASE_URL = os.getenv("BRICK_API_BASE_URL", "https://biomics.dcs.cloud/api")
API_TIMEOUT = float(os.getenv("BRICK_API_TIMEOUT", "120"))
_TRUST_ENV_PROXY_VALUE = os.getenv(
    "BRICK_TRUST_ENV_PROXY",
    os.getenv("BRICK_API_TRUST_ENV_PROXY", "0"),
)
API_TRUST_ENV_PROXY = _TRUST_ENV_PROXY_VALUE.lower() in {
    "1",
    "true",
    "yes",
    "on",
}
MAX_MULTI_HOP = 3


def _post_api(url: str, payload: dict) -> requests.Response:
    """POST to the brick API without inheriting local proxy env by default."""
    session = requests.Session()
    session.trust_env = API_TRUST_ENV_PROXY
    return session.post(url, json=payload, timeout=API_TIMEOUT)


def _extract_entities_from_text(text: str) -> List[str]:
    """
    Extract entity names from text, handling various formats.

    Supports:
    - Comma-separated: "GeneA, GeneB, GeneC"
    - With parentheses: "GeneA (protein), GeneB (gene)"
    - With descriptions: "GeneA - kinase, GeneB"
    - One per line
    - Mixed formats

    Args:
        text: Raw text containing entity names.

    Returns:
        List of cleaned entity names.
    """
    if not text or not isinstance(text, str):
        return []

    # Split by common delimiters (comma, newline, tab)
    # First try to split by newline for one-per-line format
    lines = text.replace('\r\n', '\n').replace('\r', '\n').split('\n')

    entities = []
    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Split line by comma if it contains commas
        parts = [p.strip() for p in line.split(',')]

        for part in parts:
            if not part:
                continue
            part = part.lstrip("\ufeff")

            # Extract entity name by removing common annotations
            # Pattern 1: "Entity (description)" -> "Entity"
            # Pattern 2: "Entity - description" -> "Entity"
            # Pattern 3: "Entity [id]" -> "Entity"
            # Pattern 4: "1. Entity" -> "Entity" (numbered list)

            # Remove leading numbers (e.g., "1. GeneA" or "1) GeneA")
            part = re.sub(r'^\s*[\d]+[.\)]\s*', '', part)

            # Remove content in parentheses, brackets, braces
            part = re.sub(r'\s*[\(\[\{].*?[\)\]\}]', '', part)

            # Remove prose annotations only when the separator is surrounded by
            # whitespace. Colons without surrounding whitespace are part of
            # canonical identifiers such as HGVS, GO:0008150, and MONDO:0005148.
            part = re.sub(r'\s+[\u2013\u2014:;]\s+.*$', '', part)

            # Clean up whitespace
            entity = part.strip()

            # Skip if empty or too short (likely not a valid entity)
            if len(entity) < 2:
                continue

            entities.append(entity)

    # Remove duplicates while preserving order
    seen = set()
    unique_entities = []
    for e in entities:
        key = e.lower()
        if key not in seen:
            seen.add(key)
            unique_entities.append(e)

    return unique_entities


def _extract_entities_from_lines(text: str) -> List[str]:
    """Load exact entity literals from newline-delimited text.

    Each non-empty line is one entity. Only surrounding whitespace and an
    optional UTF-8 BOM are removed; internal spaces, commas, colons,
    parentheses, and HGVS punctuation are preserved verbatim.
    """
    if not text or not isinstance(text, str):
        return []

    entities = []
    seen = set()
    for raw_line in text.replace('\r\n', '\n').replace('\r', '\n').split('\n'):
        entity = raw_line.strip().lstrip("\ufeff")
        if not entity:
            continue
        key = entity.casefold()
        if key in seen:
            continue
        seen.add(key)
        entities.append(entity)
    return entities


def _load_entities_from_file(filepath: str) -> List[str]:
    """
    Load and extract entity names from a file.

    Args:
        filepath: Path to the file containing entity names.

    Returns:
        List of cleaned entity names.

    Raises:
        FileNotFoundError: If file doesn't exist.
        ValueError: If no valid entities found.
    """
    try:
        with open(filepath, 'r', encoding='utf-8-sig') as f:
            content = f.read()
    except UnicodeDecodeError:
        # Try with different encoding
        with open(filepath, 'r', encoding='latin-1') as f:
            content = f.read()

    entities = _extract_entities_from_lines(content)

    if not entities:
        raise ValueError(f"No valid entities found in file: {filepath}")

    return entities


def _parse_comma_list(value: Optional[str]) -> Optional[List[str]]:
    """Parse a comma-separated CLI value into a clean list."""
    if not value:
        return None
    items = [item.strip() for item in value.split(",") if item.strip()]
    return items or None


def looks_like_canonical_entity(text: str) -> bool:
    if not text:
        return False
    value = text.strip()
    id_patterns = [
        r"^rs\d+$",
        r"^GO:\d+$",
        r"^CHEBI:\d+$",
        r"^HP:\d+$",
        r"^DOID:\d+$",
        r"^MONDO:\d+$",
        r"^PMID:\d+$",
        r"^[A-Z0-9]{2,10}$",
    ]
    if any(re.match(pattern, value, re.I) for pattern in id_patterns):
        return True
    lower = value.lower()
    natural_language_markers = {
        "which", "what", "how", "new", "version", "challenge", "trial",
        "cohort", "compared", "proportion", "function", "classification",
        "subtype", "profiles", "outcome", "severity",
    }
    return len(value.split()) <= 3 and not any(marker in lower for marker in natural_language_markers)


def decide_target_mode(target_entities: Optional[List[str]], target_mode: str) -> tuple[Optional[List[str]], str]:
    if not target_entities:
        return None, "none"
    if target_mode == "soft":
        return None, "soft"
    if target_mode == "strict":
        return target_entities, "strict"
    if all(looks_like_canonical_entity(entity) for entity in target_entities):
        return target_entities, "auto_strict"
    return None, "auto_soft"


def _expand_or_relations(items: Optional[List[str]]) -> Optional[List[str]]:
    """Split ``|``-joined OR relation tokens into individual relation names.

    ``is_a|develop_into`` is a server-side Cypher OR pattern; each component
    must be validated/alias-resolved on its own, otherwise the whole token is
    treated as one unknown relation and silently dropped.
    """
    if not items:
        return items
    expanded: List[str] = []
    for item in items:
        parts = [part.strip() for part in str(item).split("|") if part.strip()]
        expanded.extend(parts or [item])
    return expanded or None


def _rejoin_or_groups(originals: Optional[List[str]], resolved: List[str]) -> Optional[List[str]]:
    """Re-join validated relations back into ``|`` OR groups.

    ``originals`` keeps the user's ``|`` grouping; ``resolved`` is the flat
    validated list in the same order. Aliases that map one name to several
    replacements are joined with ``|`` (they are equivalent alternatives).
    """
    if not originals:
        return resolved or None
    groups: List[str] = []
    cursor = 0
    for token in originals:
        width = len([part for part in str(token).split("|") if part.strip()]) or 1
        members = resolved[cursor:cursor + width]
        cursor += width
        if members:
            groups.append("|".join(members))
    groups.extend(resolved[cursor:])
    return groups or None


def resolve_query_relations(
    relations: Optional[Union[str, List[str]]],
    schema_file: Optional[str] = None,
) -> tuple[Optional[List[str]], List[str], dict]:
    if relations is None:
        return None, [], {}
    raw_relations = relations
    if isinstance(raw_relations, str):
        raw_relations = [raw_relations]
    originals = list(raw_relations)
    flat = _expand_or_relations(originals) or []
    if resolve_relation_aliases is None:
        return flat or None, [], {}
    resolved, unknown, aliases = resolve_relation_aliases(flat, schema_file=schema_file)
    if unknown:
        unknown_set = {str(item).casefold() for item in unknown}
        resolved = [rel for rel in resolved if str(rel).casefold() not in unknown_set]
    return (_rejoin_or_groups(originals, resolved) or None), unknown, aliases


def _as_entity_list(value: Optional[Union[str, List[str]]]) -> Optional[List[str]]:
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    return list(value)


def build_fallback_query_plan(
    source_entity_set: Union[str, List[str]],
    target_entity_set: Optional[Union[str, List[str]]] = None,
    relation: Optional[Union[str, List[str]]] = None,
    source_entity_type: Optional[str] = None,
    target_entity_type: Optional[str] = None,
    target_mode: str = "auto",
    kg_schema: Optional[str] = None,
    enable_fallback: bool = True,
) -> dict:
    """Build schema-aware fallback stages without calling the KG API."""
    source_entities = _as_entity_list(source_entity_set) or []
    target_entities = _as_entity_list(target_entity_set)
    target_entities_for_query, target_query_mode = decide_target_mode(target_entities, target_mode)
    resolved_relation, unknown_relations, relation_aliases = resolve_query_relations(relation, schema_file=kg_schema)

    stages = [
        {
            "stage": "initial",
            "source_entity_set": source_entities,
            "target_entity_set": target_entities_for_query,
            "relation": resolved_relation,
            "source_entity_type": source_entity_type,
            "target_entity_type": target_entity_type,
        }
    ]
    if enable_fallback and resolved_relation is not None:
        stages.append({
            "stage": "fallback_no_relation",
            "source_entity_set": source_entities,
            "target_entity_set": target_entities_for_query,
            "relation": None,
            "source_entity_type": source_entity_type,
            "target_entity_type": target_entity_type,
        })
    if enable_fallback and target_entities_for_query is not None and target_entity_type:
        stages.append({
            "stage": "fallback_relax_target_type",
            "source_entity_set": source_entities,
            "target_entity_set": target_entities_for_query,
            "relation": resolved_relation,
            "source_entity_type": source_entity_type,
            "target_entity_type": None,
        })
    if enable_fallback and target_entities_for_query is not None:
        stages.append({
            "stage": "fallback_drop_target",
            "source_entity_set": source_entities,
            "target_entity_set": None,
            "relation": resolved_relation,
            "source_entity_type": source_entity_type,
            "target_entity_type": target_entity_type,
        })
    if enable_fallback and target_entities_for_query is not None and target_entity_type and resolved_relation is not None:
        stages.append({
            "stage": "fallback_drop_target_no_relation",
            "source_entity_set": source_entities,
            "target_entity_set": None,
            "relation": None,
            "source_entity_type": source_entity_type,
            "target_entity_type": target_entity_type,
        })

    preview_fields = ["stage", "target_entity_set", "relation", "source_entity_type", "target_entity_type"]
    return {
        "source_entities_preview": source_entities[:10],
        "target_entities_preview": (target_entities or [])[:10],
        "target_entities_for_query_preview": (target_entities_for_query or [])[:10],
        "target_query_mode": target_query_mode,
        "relations": resolved_relation,
        "unknown_relations": unknown_relations,
        "relation_aliases": relation_aliases,
        "fallback_enabled": enable_fallback,
        "stages": [{key: stage.get(key) for key in preview_fields} for stage in stages],
        "_stages": stages,
    }


def query_relation_with_fallback(
    source_entity_set: Union[str, List[str]],
    target_entity_set: Optional[Union[str, List[str]]] = None,
    relation: Optional[Union[str, List[str]]] = None,
    multi_hop: Union[int, tuple] = 1,
    directed: bool = False,
    source_entity_type: Optional[str] = None,
    target_entity_type: Optional[str] = None,
    relation_info_source: Optional[Union[str, List[str]]] = None,
    query_attribution: str = "name",
    limit_number: Optional[int] = None,
    return_type: str = "dataframe",
    batch_size: int = 100,
    kg_schema: Optional[str] = None,
    target_mode: str = "auto",
    enable_fallback: bool = True,
    return_diagnostics: bool = False,
    stop_condition: str = "nonempty",
    candidate_evaluator: Optional[Callable[[pd.DataFrame], bool]] = None,
) -> Union[pd.DataFrame, tuple[pd.DataFrame, dict]]:
    """Run schema-aware KG retrieval with conservative fallback stages.

    By default, the first non-empty stage is selected. MCQ/list-style callers
    can use ``stop_condition="candidate_eligible"`` with a packet-level
    evaluator so weak/background rows do not prematurely stop fallback.
    Empty results and errors are recorded as diagnostics; empty KG is still
    not treated as negative evidence.
    """
    if stop_condition not in {"nonempty", "candidate_eligible", "all_stages"}:
        raise ValueError("stop_condition must be nonempty, candidate_eligible, or all_stages")
    plan = build_fallback_query_plan(
        source_entity_set=source_entity_set,
        target_entity_set=target_entity_set,
        relation=relation,
        source_entity_type=source_entity_type,
        target_entity_type=target_entity_type,
        target_mode=target_mode,
        kg_schema=kg_schema,
        enable_fallback=enable_fallback,
    )
    stage_results = []
    selected_df = pd.DataFrame()
    selected_meta = None
    nonempty_frames = []
    stopped_by = None

    multi_hop = _normalize_multi_hop(multi_hop)
    for stage in plan["_stages"]:
        try:
            df_stage = query_relation(
                source_entity_set=stage["source_entity_set"],
                target_entity_set=stage["target_entity_set"],
                relation=stage["relation"],
                multi_hop=multi_hop,
                directed=directed,
                source_entity_type=stage["source_entity_type"],
                target_entity_type=stage["target_entity_type"],
                relation_info_source=relation_info_source,
                query_attribution=query_attribution,
                limit_number=limit_number,
                return_type=return_type,
                batch_size=batch_size,
                kg_schema=kg_schema,
            )
            error = None
        except Exception as exc:
            df_stage = pd.DataFrame()
            error = repr(exc)

        if not df_stage.empty:
            df_stage = df_stage.copy()
            df_stage["query.stage"] = stage["stage"]
            df_stage["query.target_query_mode"] = plan["target_query_mode"]
            nonempty_frames.append(df_stage)

        meta = {
            "stage": stage["stage"],
            "row_count": int(len(df_stage)),
            "relation": stage["relation"],
            "target_present": bool(stage["target_entity_set"]),
            "source_type": stage["source_entity_type"],
            "target_type": stage["target_entity_type"],
            "error": error,
        }
        stage_results.append(meta)
        if selected_meta is None or (selected_meta.get("error") and not error):
            selected_df = df_stage
            selected_meta = meta
        if not df_stage.empty:
            combined_df = pd.concat(nonempty_frames, ignore_index=True) if len(nonempty_frames) > 1 else df_stage
            selected_df = combined_df
            selected_meta = meta
            if stop_condition == "nonempty":
                stopped_by = "nonempty"
                break
            if stop_condition == "candidate_eligible":
                eligible = bool(candidate_evaluator(combined_df)) if candidate_evaluator is not None else False
                meta["candidate_eligible"] = eligible
                if eligible:
                    stopped_by = "candidate_eligible"
                    break

    selected_meta = selected_meta or {"stage": None, "row_count": 0, "error": None}
    if stopped_by is None:
        stopped_by = "exhausted_stages"
        if stop_condition == "candidate_eligible" and candidate_evaluator is None and nonempty_frames:
            stopped_by = "exhausted_stages_no_candidate_evaluator"
    selected_nonempty_stages = [item.get("stage") for item in stage_results if item.get("row_count", 0) > 0]
    diagnostics = {
        "selected_stage": selected_meta.get("stage"),
        "selected_nonempty_stages": selected_nonempty_stages,
        "stop_condition": stop_condition,
        "stopped_by": stopped_by,
        "stage_results": stage_results,
        "target_query_mode": plan["target_query_mode"],
        "relation_aliases": plan["relation_aliases"],
        "unknown_relations": plan["unknown_relations"],
        "kg_query_calls": len(stage_results),
        "kg_nonempty_calls": sum(1 for item in stage_results if item.get("row_count", 0) > 0),
        "kg_error_calls": sum(1 for item in stage_results if item.get("error")),
        "initial_nonempty": any(item.get("stage") == "initial" and item.get("row_count", 0) > 0 for item in stage_results),
        "fallback_nonempty": any(str(item.get("stage", "")).startswith("fallback_") and item.get("row_count", 0) > 0 for item in stage_results),
        "fallback_drop_target_nonempty": any(item.get("stage") == "fallback_drop_target" and item.get("row_count", 0) > 0 for item in stage_results),
        "fallback_no_relation_nonempty": any(item.get("stage") == "fallback_no_relation" and item.get("row_count", 0) > 0 for item in stage_results),
        "fallback_enabled": enable_fallback,
    }
    return (selected_df, diagnostics) if return_diagnostics else selected_df


def _normalize_multi_hop(multi_hop: Union[int, List[int], tuple]) -> Union[int, List[int]]:
    """Validate and normalize KG hop count, capped at 3 hops."""
    is_range = isinstance(multi_hop, (list, tuple))
    raw_values = list(multi_hop) if is_range else [multi_hop]
    if not raw_values:
        raise ValueError("multi_hop must contain at least one hop value")

    try:
        hops = [int(value) for value in raw_values]
    except (TypeError, ValueError) as exc:
        raise ValueError("multi_hop must be an integer or a list/range of integers") from exc

    invalid = [hop for hop in hops if hop < 1 or hop > MAX_MULTI_HOP]
    if invalid:
        raise ValueError(f"multi_hop values must be between 1 and {MAX_MULTI_HOP}: {invalid}")

    return hops if is_range else hops[0]


def _path_column_sort_key(column: str):
    match = re.match(r"^path\.(\d+)\.(.+)$", column)
    if not match:
        return (10_000, column)
    field_order = {
        "id": 0,
        "name": 1,
        "type": 2,
        "relation": 3,
        "info_source": 4,
        "info_source_length": 5,
        "relation_confidence": 6,
        "condition": 7,
        "original_relation": 8,
    }
    return (int(match.group(1)), field_order.get(match.group(2), 99), column)


def _path_key_columns(df: pd.DataFrame) -> List[str]:
    """Choose path columns that preserve distinct one-hop or multi-hop paths."""
    candidates = [
        col for col in df.columns
        if re.match(r"^path\.\d+\.(?:id|name|relation|info_source|relation_confidence|original_relation|condition)$", col)
    ]
    return sorted(candidates, key=_path_column_sort_key)


def _drop_duplicate_paths(df: pd.DataFrame, key_columns: List[str]) -> pd.DataFrame:
    """Deduplicate paths without mutating list/dict-valued evidence columns."""
    if not key_columns:
        return df

    def make_hashable(value):
        try:
            hash(value)
            return value
        except TypeError:
            return json.dumps(value, sort_keys=True, ensure_ascii=False, default=str)

    keys = df[key_columns].copy()
    for column in key_columns:
        keys[column] = keys[column].map(make_hashable)
    return df.loc[~keys.duplicated()].copy()


def _normalize_entity_type_filter(value):
    """Serialize one or more API entity types into the supported pipe syntax."""
    if value is None or isinstance(value, str):
        return value
    items = []
    seen = set()
    for item in value:
        item = str(item).strip()
        if not item or item in seen:
            continue
        seen.add(item)
        items.append(item)
    return "|".join(items) or None


def _preferred_output_columns(df: pd.DataFrame, hops: int = 1) -> List[str]:
    """Return readable path columns for one-hop and optional multi-hop results."""
    preferred = []
    max_path_index = 2 * max(1, min(int(hops or 1), MAX_MULTI_HOP))
    for idx in range(max_path_index + 1):
        if idx % 2 == 0:
            fields = ("name", "type", "id")
        else:
            fields = ("relation", "info_source", "info_source_length", "relation_confidence", "condition", "original_relation")
        for field in fields:
            col = f"path.{idx}.{field}"
            if col in df.columns and col not in preferred:
                preferred.append(col)

    if preferred:
        return preferred

    path_cols = sorted(
        [col for col in df.columns if col.startswith("path.")],
        key=_path_column_sort_key,
    )
    return path_cols[:18] or list(df.columns[:12])


def query_cypher(cypher: str, parameters: Optional[dict] = None) -> pd.DataFrame:
    """
    Execute a Cypher query via brick KG API.

    Args:
        cypher: The Cypher query string.
        parameters: Query parameters dictionary.

    Returns:
        pandas.DataFrame: Query results with flattened columns.
    """
    url = f"{API_BASE_URL}/v1/cypher"
    payload = {
        "cypher": cypher,
        "parameters": parameters or {},
        "return_type": "dataframe"
    }

    response = _post_api(url, payload)
    response.raise_for_status()

    result = response.json()
    if not result.get("success"):
        raise RuntimeError(f"API error: {result.get('message', 'Unknown error')}")

    data = result.get("data", {})
    rows = data.get("rows", [])
    columns = data.get("columns", [])

    if not rows:
        return pd.DataFrame(columns=columns)

    return pd.DataFrame(rows, columns=columns)


def query_shortest_path_api(
    source_entities: List[str],
    target_entities: List[str],
    relation: Optional[str] = None,
    multi_hop: Union[int, tuple] = (1, 10),
    directed: bool = False,
    source_entity_type: Optional[str] = None,
    target_entity_type: Optional[str] = None,
    limit_number: Optional[int] = None,
) -> pd.DataFrame:
    """
    Server-side per-pair shortest path via the /v1/shortest_path endpoint.

    One POST per source entity with the full target list, replicating legacy
    ``BRICK.qr.query_shortest_path`` semantics (Neo4j shortestPath per pair),
    avoiding the LIMIT-truncation / timeout trade-off of /v1/relation.
    """
    url = f"{API_BASE_URL}/v1/shortest_path"
    hop_payload = list(multi_hop) if isinstance(multi_hop, tuple) else multi_hop
    frames = []
    seen_pairs = set()  # unordered pairs already queried (shortestPath is undirected)
    for src in source_entities:
        pending = []
        for tgt in target_entities:
            key = frozenset((src, tgt))
            if tgt == src or key in seen_pairs:
                continue
            seen_pairs.add(key)
            pending.append(tgt)
        if not pending:
            continue
        payload = {
            "source_entity_name": src,
            "target_entity_name": pending,
            "relation": relation,
            "multi_hop": hop_payload,
            "directed": directed,
            "source_entity_type": source_entity_type,
            "target_entity_type": target_entity_type,
            "return_type": "dataframe",
        }
        if limit_number is not None:
            payload["limit_number"] = limit_number
        response = _post_api(url, payload)
        response.raise_for_status()
        result = response.json()
        if not result.get("success"):
            raise RuntimeError(f"API error for source '{src}': {result.get('message', 'Unknown error')}")
        data = result.get("data", {})
        rows = data.get("rows", [])
        frames.append(pd.DataFrame(rows, columns=data.get("columns", [])))
        print(f"  shortest_path[{src}]: {len(rows)} row(s)")
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def query_relation(
    source_entity_set: Union[str, List[str]],
    target_entity_set: Optional[Union[str, List[str]]] = None,
    relation: Optional[Union[str, List[str]]] = None,
    multi_hop: Union[int, tuple] = 1,
    directed: bool = False,
    source_entity_type: Optional[str] = None,
    target_entity_type: Optional[str] = None,
    relation_info_source: Optional[Union[str, List[str]]] = None,
    query_attribution: str = "name",
    limit_number: Optional[int] = None,
    return_type: str = "dataframe",
    batch_size: int = 100,
    kg_schema: Optional[str] = None,
) -> pd.DataFrame:
    """
    Query relationships between entities via brick KG API.

    Supports two modes:
    - Neighbor mode (target_entity_set=None): finds all connected entities.
    - Relation mode (target_entity_set provided): validates specific pairs.

    When source_entity_set has more than batch_size entities, the query is
    automatically split into multiple batches to avoid server timeouts.

    Args:
        source_entity_set: Source entity or list of entities.
        target_entity_set: Target entity or list (None for neighbor mode).
        relation: Relationship type(s) to filter.
        multi_hop: Number of path hops to query. Values are capped at 3.
        directed: Whether to query directed paths.
        source_entity_type: Type of source entities (e.g., 'Cell', 'Gene').
        target_entity_type: Type of target entities.
        relation_info_source: Filter by info source.
        query_attribution: Attribute to match ('name' or 'id').
        limit_number: Maximum results to return per batch.
        return_type: Return format ('dataframe' or 'graph').
        batch_size: Maximum number of source entities per batch (default: 100).

    Returns:
        pandas.DataFrame: Query results with relationship data.
    """
    # Normalize source_entity_set to list
    if isinstance(source_entity_set, str):
        source_entity_set = [source_entity_set]
    else:
        source_entity_set = list(source_entity_set)

    # Normalize target_entity_set
    if target_entity_set is not None:
        if isinstance(target_entity_set, str):
            target_entity_set = [target_entity_set]
        else:
            target_entity_set = list(target_entity_set)

    # Resolve relation aliases and drop blocked/unknown relation filters before
    # sending a live KG request. This keeps concept labels such as
    # associated_with from becoming hard API filters.
    relation, _, _ = resolve_query_relations(relation, schema_file=kg_schema)

    # Normalize edge info_source filters to a list.
    if relation_info_source is not None and expand_info_source_filter is not None:
        relation_info_source = expand_info_source_filter(relation_info_source)
    elif relation_info_source is not None and isinstance(relation_info_source, str):
        relation_info_source = _parse_comma_list(relation_info_source) or [relation_info_source]
    elif relation_info_source is not None:
        relation_info_source = list(relation_info_source)

    multi_hop = _normalize_multi_hop(multi_hop)

    # If source entities are within batch size, query directly
    if len(source_entity_set) <= batch_size:
        return _query_relation_single_batch(
            source_entity_set=source_entity_set,
            target_entity_set=target_entity_set,
            relation=relation,
            multi_hop=multi_hop,
            directed=directed,
            source_entity_type=source_entity_type,
            target_entity_type=target_entity_type,
            relation_info_source=relation_info_source,
            query_attribution=query_attribution,
            limit_number=limit_number,
            return_type=return_type
        )

    # Split into batches and query sequentially
    print(f"Source entities ({len(source_entity_set)}) exceed batch size ({batch_size}). Splitting into batches...")

    all_results = []
    total_batches = (len(source_entity_set) + batch_size - 1) // batch_size

    for i in range(0, len(source_entity_set), batch_size):
        batch_num = i // batch_size + 1
        batch = source_entity_set[i:i + batch_size]

        print(f"  Querying batch {batch_num}/{total_batches} ({len(batch)} entities)...")

        try:
            df = _query_relation_single_batch(
                source_entity_set=batch,
                target_entity_set=target_entity_set,
                relation=relation,
                multi_hop=multi_hop,
                directed=directed,
                source_entity_type=source_entity_type,
                target_entity_type=target_entity_type,
                relation_info_source=relation_info_source,
                query_attribution=query_attribution,
                limit_number=limit_number,
                return_type=return_type
            )

            if not df.empty:
                all_results.append(df)
                print(f"    Found {len(df)} relations")
            else:
                print(f"    No results")

        except Exception as e:
            print(f"    Error in batch {batch_num}: {e}")
            # Continue with next batch instead of failing entirely

    if not all_results:
        return pd.DataFrame()

    # Combine all batch results and remove duplicates
    combined = pd.concat(all_results, ignore_index=True)

    # Remove duplicate paths while preserving distinct multi-hop routes.
    available_key_cols = _path_key_columns(combined)
    if available_key_cols:
        combined = combined.drop_duplicates(subset=available_key_cols)

    print(f"\nTotal unique relations from all batches: {len(combined)}")
    return combined


def _query_relation_single_batch(
    source_entity_set: List[str],
    target_entity_set: Optional[List[str]] = None,
    relation: Optional[List[str]] = None,
    multi_hop: Union[int, List[int]] = 1,
    directed: bool = False,
    source_entity_type: Optional[str] = None,
    target_entity_type: Optional[str] = None,
    relation_info_source: Optional[Union[str, List[str]]] = None,
    query_attribution: str = "name",
    limit_number: Optional[int] = None,
    return_type: str = "dataframe"
) -> pd.DataFrame:
    """
    Execute a single batch query to the brick KG API.

    This is the internal implementation that sends one HTTP request.
    Use query_relation() for automatic batching support.
    """
    url = f"{API_BASE_URL}/v1/relation"

    payload = {
        "source_entity_set": source_entity_set,
        "target_entity_set": target_entity_set,
        "relation": relation,
        "multi_hop": multi_hop,
        "directed": directed,
        "source_entity_type": source_entity_type,
        "target_entity_type": target_entity_type,
        # The live /v1/relation endpoint expects the edge provenance filter as
        # "info_source". Sending "relation_info_source" currently triggers a
        # server-side Cypher type error.
        "info_source": relation_info_source,
        "query_attribution": query_attribution,
        "limit_number": limit_number,
        "return_type": return_type
    }

    # Remove None values
    payload = {k: v for k, v in payload.items() if v is not None}

    response = _post_api(url, payload)
    response.raise_for_status()

    result = response.json()
    if not result.get("success"):
        raise RuntimeError(f"API error: {result.get('message', 'Unknown error')}")

    data = result.get("data", {})
    rows = data.get("rows", [])
    columns = data.get("columns", [])

    if not rows:
        return pd.DataFrame(columns=columns)

    return pd.DataFrame(rows, columns=columns)



def _apply_kg_query_policy(args, relations, info_sources):
    """Resolve source/relation/hop policy without overwriting user filters."""
    if info_sources and canonicalize_info_sources is not None:
        info_sources = canonicalize_info_sources(info_sources)
    question = getattr(args, "question", None)
    policy_mode = getattr(args, "kg_policy", None)
    if policy_mode is None:
        policy_mode = "auto" if question else "manual"
    policy_plan = None
    raw_relations = list(relations or [])
    relations_explicitly_dropped = False
    if relations and resolve_relation_aliases is not None:
        # Expand '|' OR tokens (e.g. is_a|develop_into) so each component is
        # validated/alias-resolved on its own, then re-join the validated
        # relations into their original OR groups for the API call.
        flat_relations = _expand_or_relations(relations) or []
        resolved_relations, unknown_relations, relation_aliases = resolve_relation_aliases(
            flat_relations,
            schema_file=getattr(args, "kg_schema", None),
        )
        if unknown_relations:
            unknown_set = {str(item).casefold() for item in unknown_relations}
            resolved_relations = [
                rel for rel in resolved_relations if str(rel).casefold() not in unknown_set
            ]
        relations = _rejoin_or_groups(raw_relations, resolved_relations)
        relations_explicitly_dropped = bool(raw_relations and relation_aliases and not relations)
        policy_plan = {
            "relation_aliases": relation_aliases,
            "unknown_relations": unknown_relations,
            "warnings": [],
        }
        if relation_aliases:
            policy_plan["warnings"].append(f"Resolved relation aliases: {relation_aliases}")
        if unknown_relations:
            policy_plan["warnings"].append(f"Unknown KG relation(s): {unknown_relations}")
            print(
                f"Warning: unknown KG relation(s) dropped from query filter: {unknown_relations}",
                file=sys.stderr,
            )
    if select_kg_policy is None:
        return relations, info_sources, args.hops, policy_plan

    should_route = policy_mode != "manual" or bool(question)
    if not should_route:
        return relations, info_sources, args.hops, policy_plan

    routed_policy = select_kg_policy(
        question=question,
        source_type=args.st,
        target_type=args.tt,
        relation=relations,
        task=getattr(args, "task", None),
        kg_policy=policy_mode,
        policy_file=getattr(args, "policy_file", None),
        schema_file=getattr(args, "kg_schema", None),
        requested_hops=args.hops,
        user_info_sources=info_sources,
    )
    if policy_plan:
        routed_policy["warnings"] = policy_plan.get("warnings", []) + routed_policy.get("warnings", [])
        routed_policy["relation_aliases"] = {**policy_plan.get("relation_aliases", {}), **routed_policy.get("relation_aliases", {})}
        routed_policy["unknown_relations"] = policy_plan.get("unknown_relations", []) or routed_policy.get("unknown_relations", [])
    policy_plan = routed_policy
    if not policy_plan.get("use_kg", True):
        raise ValueError("KG query disabled by --kg-policy off; do not call brick qr for this step.")

    if relations_explicitly_dropped:
        routed_policy["relations"] = []
    if not relations and not relations_explicitly_dropped:
        relations = policy_plan.get("relations") or None
    if not info_sources:
        info_sources = policy_plan.get("info_sources") or None
    hops = int(policy_plan.get("hops") or args.hops or 1)

    if getattr(args, "show_policy", False):
        print("KG query policy:")
        for key in ["kg_policy", "intent", "reason", "evidence_tier", "relations", "info_sources", "hops", "co_mentioned_policy"]:
            if key in policy_plan:
                print(f"  {key}: {policy_plan[key]}")
        for warning in policy_plan.get("warnings", []):
            print(f"  warning: {warning}")
    return relations, info_sources, hops, policy_plan


def _is_broad_query(target_entities, relations, info_sources, target_type, hops):
    return (
        target_entities is None
        and not relations
        and not info_sources
        and not target_type
        and int(hops or 1) > 0
    )


def _broad_query_message():
    return (
        "Broad KG neighbor query: no target entity/type, relation, or info_source filter. "
        "Prefer --tt, -r, --info-source, or --kg-policy auto; pass --allow-broad-query to silence this guardrail."
    )

def _cli_main():
    """CLI entry point for brick qr command."""
    parser = argparse.ArgumentParser(
        description="brick qr - Query entity relations",
        usage='brick qr "entity_name" [options]'
    )

    parser.add_argument(
        "source",
        nargs="*",
        help="Source entity name; adjacent unquoted words are rejoined, or use a legacy comma-separated list",
    )
    parser.add_argument("-i", "--input-file", help="File containing exactly one source entity literal per line")
    parser.add_argument(
        "-t",
        "--target",
        nargs="+",
        help="Target entity name; adjacent unquoted words are rejoined, or use a legacy comma-separated list",
    )
    parser.add_argument("--target-file", help="File containing exactly one target entity literal per line")
    parser.add_argument("-r", "--relation", help="Relation type (e.g., produce, causally_related_to). Comma-separated values are AND filters; use '|' for server-side OR (e.g., 'is_a|develop_into'). Concept aliases such as marker_of are resolved.")
    parser.add_argument("--st", help="Source entity type (e.g., Cell, Gene, Disease)")
    parser.add_argument("--tt", help="Target entity type")
    parser.add_argument(
        "--info-source",
        "--relation-info-source",
        dest="info_source",
        help="Comma-separated edge info_source filter for relationship provenance (e.g., CellMarker,PanglaoDB)",
    )
    parser.add_argument("--question", help="Natural-language question/context used by --kg-policy auto to select sources/relations")
    parser.add_argument("--task", help="Task hint for KG source routing, e.g. celltype, drug, enrichment, trajectory")
    parser.add_argument("--kg-policy", "--source-policy", dest="kg_policy", choices=KG_POLICY_CHOICES + ("manual",), default=None,
                        help="KG query policy: manual uses explicit filters; auto/source-scoped infer source filters; multihop allows limited path search. Defaults to auto when --question is provided, otherwise manual.")
    parser.add_argument("--policy-file", help="Optional JSON KG source policy catalog")
    parser.add_argument("--kg-schema", default=None, help="Path to KG schema JSON. Defaults to references/kg_schema.json.")
    parser.add_argument("--show-schema", action="store_true", help="Show resolved KG schema decisions in dry-run output.")
    parser.add_argument("--target-mode", choices=["auto", "strict", "soft"], default="auto", help="Target filtering mode. auto hard-filters only canonical-looking targets.")
    parser.add_argument("--allow-broad-query", action="store_true",
                        help="Allow an unscoped query with no target type, relation, or info_source filter")
    parser.add_argument("--show-policy", action="store_true", help="Print the resolved KG query policy before querying")
    parser.add_argument("--dry-run-plan", action="store_true",
                        help="Print the resolved KG query policy and API filters without calling the KG API")
    parser.add_argument("--no-fallback", dest="enable_fallback", action="store_false", default=True,
                        help="Disable schema-aware fallback stages after an empty KG query")
    parser.add_argument(
        "--hops",
        "--multi-hop",
        dest="hops",
        type=int,
        default=1,
        metavar="{1,2,3}",
        help="Number of KG path hops to query: 1 direct only, 2/3 optional multi-hop (default: 1). For trajectory/shortest-path use cases pair deeper hops (e.g. --hops 10) with --shortest (server-side per-pair shortestPath; no LIMIT truncation or 60s-timeout risk).",
    )
    parser.add_argument("--directed", action="store_true", help="Respect edge direction during relation/path query")
    parser.add_argument(
        "--shortest",
        action="store_true",
        help="Compute the shortest path per source-target pair on the server via the /v1/shortest_path endpoint (one call per source entity, targets passed as a list; unordered pairs only). Exact legacy BRICK.qr.query_shortest_path semantics; avoids the LIMIT truncation and 60s-timeout issues of the generic /v1/relation endpoint.",
    )
    parser.add_argument("-n", type=int, default=10, help="Number of results (default: 10)")
    parser.add_argument("-o", help="Output CSV file")

    args = parser.parse_args()
    if args.hops < 1 or (args.hops > MAX_MULTI_HOP and not args.shortest):
        parser.error(f"--hops must be between 1 and {MAX_MULTI_HOP} (deeper hops require --shortest)")
    source_text = " ".join(args.source).strip() if args.source else None
    target_text = " ".join(args.target).strip() if args.target else None

    # Validate input: either source or input-file must be provided
    if not source_text and not args.input_file:
        parser.error("Either 'source' argument or '-i/--input-file' must be provided")

    try:
        # Load source entities
        if args.input_file:
            if not args.dry_run_plan:
                print(f"Loading source entities from: {args.input_file}")
            source_entities = _load_entities_from_file(args.input_file)
            if not args.dry_run_plan:
                print(f"Extracted {len(source_entities)} source entities: {source_entities[:5]}{'...' if len(source_entities) > 5 else ''}")
        else:
            source_entities = _extract_entities_from_text(source_text)
            if not source_entities:
                raise ValueError("No valid source entities found")
            if not args.dry_run_plan:
                print(f"Source entities: {source_entities}")

        # Load target entities
        target_entities = None
        if args.target_file:
            if not args.dry_run_plan:
                print(f"Loading target entities from: {args.target_file}")
            target_entities = _load_entities_from_file(args.target_file)
            if not args.dry_run_plan:
                print(f"Extracted {len(target_entities)} target entities: {target_entities[:5]}{'...' if len(target_entities) > 5 else ''}")
        elif target_text:
            target_entities = _extract_entities_from_text(target_text)
            if target_entities and not args.dry_run_plan:
                print(f"Target entities: {target_entities}")

        target_entities_for_query, target_query_mode = decide_target_mode(target_entities, args.target_mode)

        relations = _parse_comma_list(args.relation)
        info_sources = _parse_comma_list(args.info_source)
        relations, info_sources, resolved_hops, policy = _apply_kg_query_policy(args, relations, info_sources)
        fallback_plan = build_fallback_query_plan(
            source_entity_set=source_entities,
            target_entity_set=target_entities,
            relation=relations,
            source_entity_type=args.st,
            target_entity_type=args.tt,
            target_mode=args.target_mode,
            kg_schema=args.kg_schema,
            enable_fallback=args.enable_fallback,
        )

        if _is_broad_query(target_entities_for_query, relations, info_sources, args.tt, resolved_hops):
            message = _broad_query_message()
            if not args.allow_broad_query:
                raise ValueError(message)
            print(f"Warning: {message}")

        if args.dry_run_plan:
            command = ["brick", "qr"]
            if args.input_file:
                command += ["-i", args.input_file]
            elif source_text:
                command.append(source_text)
            if target_text:
                command += ["-t", target_text]
            if args.target_file:
                command += ["--target-file", args.target_file]
            if args.st:
                command += ["--st", args.st]
            if args.tt:
                command += ["--tt", args.tt]
            if relations:
                command += ["-r", ",".join(relations)]
            if info_sources:
                command += ["--info-source", ",".join(info_sources)]
            if resolved_hops != 1:
                command += ["--hops", str(resolved_hops)]
            if args.directed:
                command.append("--directed")
            command += ["-n", str(args.n)]
            if args.o:
                command += ["-o", args.o]
            plan = {
                "source_entities_preview": source_entities[:10],
                "target_entities_preview": (target_entities or [])[:10],
                "target_entities_for_query_preview": (target_entities_for_query or [])[:10],
                "target_query_mode": target_query_mode,
                "source_entity_type": args.st,
                "target_entity_type": args.tt,
                "relations": relations,
                "info_sources": info_sources,
                "hops": resolved_hops,
                "directed": args.directed,
                "limit": args.n,
                "command": command,
                "command_sh": shlex.join(command),
                "recommended_action": "query_kg",
                "policy": policy or {},
                "fallback_enabled": args.enable_fallback,
                "fallback_stages": fallback_plan["stages"],
            }
            if args.show_schema and load_kg_schema is not None:
                schema = load_kg_schema(args.kg_schema)
                plan["schema"] = {
                    "node_types": schema.get("node_types", []),
                    "queryable_edge_types": schema.get("queryable_edge_types", []),
                    "blocked_query_edge_types": schema.get("blocked_query_edge_types", []),
                    "weak_evidence_edge_types": schema.get("weak_evidence_edge_types", []),
                    "relation_aliases": (policy or {}).get("relation_aliases", {}),
                    "unknown_relations": (policy or {}).get("unknown_relations", []),
                }
            print(json.dumps(plan, indent=2, ensure_ascii=True))
            return

        print(f"\nQuerying {len(source_entities)} source entities...")
        if policy:
            print(f"KG policy: {policy.get('kg_policy')} / intent={policy.get('intent')} / reason={policy.get('reason')}")
        if relations:
            print(f"Restricting relation values: {relations}")
        if info_sources:
            print(f"Restricting edge info_source values: {info_sources}")
        if resolved_hops > 1:
            print(f"Using optional multi-hop search: up to {resolved_hops} hops")

        fallback_diagnostics = {}
        if args.shortest:
            df = query_shortest_path_api(
                source_entities=source_entities,
                target_entities=target_entities or source_entities,
                relation=",".join(relations) if relations else None,
                multi_hop=(1, resolved_hops),
                directed=args.directed,
                source_entity_type=args.st,
                target_entity_type=args.tt,
                limit_number=args.n,
            )
            print(f"Server-side shortest_path: {len(df)} row(s) total")
        else:
            df, fallback_diagnostics = query_relation_with_fallback(
                source_entity_set=source_entities,
                target_entity_set=target_entities,
                relation=relations,
                multi_hop=resolved_hops,
                directed=args.directed,
                source_entity_type=args.st,
                target_entity_type=args.tt,
                relation_info_source=info_sources,
                limit_number=args.n,
                kg_schema=args.kg_schema,
                target_mode=args.target_mode,
                enable_fallback=args.enable_fallback,
                return_diagnostics=True,
            )
            print(f"Selected KG query stage: {fallback_diagnostics.get('selected_stage')}")

        if df.empty:
            print("No results found.")
            if fallback_diagnostics.get("kg_error_calls"):
                print("KG query errors:")
                for item in fallback_diagnostics.get("stage_results", []):
                    if item.get("error"):
                        print(f"  {item.get('stage')}: {item.get('error')}")
            return

        print(f"Found {len(df)} relations")

        if args.o:
            df.to_csv(args.o, index=False)
            print(f"Saved to {args.o}")
        else:
            available = _preferred_output_columns(df, resolved_hops)
            print(df[available].to_string())

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _cli_main()
