"""Leakage-resistant drug-discovery helpers for the unified ``brick`` CLI."""

from __future__ import annotations

import argparse
import ast
import hashlib
import inspect
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


DEFAULT_GRACE_SEEDS = (20241204, 20241205, 20241206)
CURATED_DRUG_SOURCES = {
    "drugbank",
    "drugcentral",
    "iuphar",
    "pharmgkb",
    "repodb",
    "clinicaltrials",
    "clinicaltrails",
}
PHARMACO_RELATIONS = {
    "pharmaco_effect_on",
    "treats",
    "indication",
    "contraindication",
    "agonist",
    "antagonist",
    "inhibitor",
}


def _utc_now():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _require_columns(frame, columns, context):
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{context} is missing required columns: {', '.join(missing)}")


def _normalize_label(value):
    return " ".join(str(value).strip().casefold().replace("_", " ").replace("-", " ").split())


def _flatten_values(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    if isinstance(value, (list, tuple, set)):
        result = []
        for item in value:
            result.extend(_flatten_values(item))
        return result
    if isinstance(value, dict):
        result = []
        for item in value.values():
            result.extend(_flatten_values(item))
        return result
    text = str(value).strip()
    if not text:
        return []
    if text[:1] in "[{(" and text[-1:] in "]})":
        try:
            return _flatten_values(ast.literal_eval(text))
        except (SyntaxError, ValueError):
            pass
    for separator in ("|", ";", ","):
        if separator in text:
            return [part.strip() for part in text.split(separator) if part.strip()]
    return [text]


def _scope_label(relations, sources):
    relations = {_normalize_label(value).replace(" ", "_") for value in relations}
    sources = {_normalize_label(value).replace(" ", "") for value in sources}
    curated = bool(sources & {value.replace(" ", "") for value in CURATED_DRUG_SOURCES})
    pharmacologic = bool(relations & PHARMACO_RELATIONS)
    if curated and pharmacologic:
        return "curated_pharmacologic_entity"
    if curated:
        return "curated_drug_database_entity"
    if pharmacologic:
        return "pharmacologic_relation_entity"
    return "chemical_entity_unspecified"


def chemical_scope_table(frame):
    """Tag Chemical entities from evidence without changing candidate eligibility."""
    _require_columns(frame, ["path.2.id", "path.2.name", "path.2.type"], "Chemical evidence")
    chemicals = frame.loc[frame["path.2.type"].astype(str).eq("Chemical")].copy()
    if chemicals.empty:
        raise ValueError("No Chemical rows were found")
    empty_ids = chemicals["path.2.id"].isna() | chemicals["path.2.id"].astype(str).str.strip().eq("")
    if empty_ids.any():
        raise ValueError("Chemical evidence contains empty stable IDs")

    records = []
    for stable_id, group in chemicals.groupby("path.2.id", sort=True, dropna=False):
        names = sorted(set(group["path.2.name"].dropna().astype(str)))
        relations = []
        sources = []
        for column in ("path.1.relation", "path.1.original_relation"):
            if column in group:
                for value in group[column]:
                    relations.extend(_flatten_values(value))
        if "path.1.info_source" in group:
            for value in group["path.1.info_source"]:
                sources.extend(_flatten_values(value))
        records.append(
            {
                "path.2.id": str(stable_id),
                "path.2.name": names[0] if names else "",
                "path.2.type": "Chemical",
                "chemical_scope": _scope_label(relations, sources),
                "formal_candidate_eligible": True,
                "scope_is_diagnostic_only": True,
                "name_aliases": "|".join(names),
                "relations": "|".join(sorted(set(relations))),
                "info_sources": "|".join(sorted(set(sources))),
                "evidence_row_count": int(len(group)),
            }
        )
    return pd.DataFrame(records)


def cmd_scope(args):
    frame = pd.read_csv(args.input)
    result = chemical_scope_table(frame)
    result.to_csv(args.output, index=False)
    print(f"Wrote {args.output} ({len(result)} Chemical entities; diagnostic tags only)")


def cmd_extract(args):
    ranked = pd.read_csv(args.input)
    rank_col = next(
        (column for column in ("consensus_rank", "path.2.rank_voting") if column in ranked.columns),
        None,
    )
    if rank_col is None:
        raise ValueError("Ranked Chemical table is missing consensus_rank or path.2.rank_voting")
    _require_columns(
        ranked,
        ["path.2.id", "path.2.name", "path.2.type"],
        "Ranked Chemical table",
    )
    ranked = ranked.loc[ranked["path.2.type"].astype(str).eq("Chemical")].copy()
    ranked["path.2.id"] = ranked["path.2.id"].astype(str).str.strip()
    if ranked["path.2.id"].eq("").any():
        raise ValueError("Ranked Chemical table contains empty stable IDs")
    if ranked["path.2.id"].duplicated().any():
        raise ValueError("Ranked Chemical table contains duplicate stable IDs")
    ranked = ranked.sort_values(
        [rank_col, "path.2.id", "path.2.name"], kind="mergesort"
    )
    ranked.insert(0, "candidate_rank", range(1, len(ranked) + 1))
    if args.scope:
        scope = pd.read_csv(args.scope)
        _require_columns(scope, ["path.2.id", "chemical_scope"], "Chemical scope table")
        ranked = ranked.merge(
            scope.drop_duplicates("path.2.id"), on="path.2.id", how="left", suffixes=("", "_scope")
        )
    result = ranked if args.top_n <= 0 else ranked.head(args.top_n)
    result.to_csv(args.output, index=False)
    print(f"Wrote {args.output} ({len(result)} candidates)")


def freeze_ranking(ranking_path, output_path):
    ranking_path = Path(ranking_path)
    ranking = pd.read_csv(ranking_path)
    _require_columns(ranking, ["path.2.id", "path.2.name"], "Primary ranking")
    ids = ranking["path.2.id"].fillna("").astype(str).str.strip()
    if ids.eq("").any() or ids.duplicated().any():
        raise ValueError("Primary ranking must contain non-empty, unique stable IDs before freezing")
    digest = _sha256(ranking_path)
    output_path = Path(output_path)
    output_path.write_text(f"{digest}  {ranking_path.name}\n", encoding="utf-8")
    metadata = {
        "ranking": str(ranking_path.resolve()),
        "ranking_sha256": digest,
        "ranking_rows": int(len(ranking)),
        "frozen_at": _utc_now(),
    }
    _write_json(output_path.with_suffix(output_path.suffix + ".json"), metadata)
    return metadata


def cmd_freeze(args):
    metadata = freeze_ranking(args.input, args.output)
    print(f"Frozen {metadata['ranking_rows']} rows: {metadata['ranking_sha256']}")


def _read_freeze_hash(path):
    text = Path(path).read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError("Ranking freeze file is empty")
    return text.split()[0]


def cmd_validate(args):
    ranking_path = Path(args.ranking)
    reference_path = Path(args.reference)
    freeze_path = Path(args.freeze)
    expected_hash = _read_freeze_hash(freeze_path)
    actual_hash = _sha256(ranking_path)
    if expected_hash != actual_hash:
        raise ValueError("Primary ranking hash does not match the frozen SHA-256")
    if not args.allow_preexisting_reference and reference_path.stat().st_mtime_ns < freeze_path.stat().st_mtime_ns:
        raise ValueError(
            "Reference file predates the ranking freeze; use a post-freeze query or explicitly "
            "pass --allow-preexisting-reference for an externally governed evaluation set"
        )

    ranking = pd.read_csv(ranking_path)
    reference = pd.read_csv(reference_path)
    _require_columns(ranking, ["path.2.id", "path.2.name"], "Primary ranking")
    _require_columns(reference, ["path.2.id", "path.2.name"], "Reference table")
    ranking_ids = ranking["path.2.id"].dropna().astype(str).str.strip()
    reference_ids = set(reference["path.2.id"].dropna().astype(str).str.strip()) - {""}
    if ranking_ids.eq("").any() or ranking_ids.duplicated().any():
        raise ValueError("Primary ranking IDs must be non-empty and unique")
    eligible_ids = set(ranking_ids) & reference_ids

    rows = []
    for cutoff in args.cutoffs:
        evaluated = ranking.head(cutoff)
        evaluated_ids = set(evaluated["path.2.id"].astype(str))
        hits = len(evaluated_ids & reference_ids)
        denominator = len(evaluated)
        eligible_count = len(eligible_ids)
        all_count = len(reference_ids)
        rows.append(
            {
                "cutoff": int(cutoff),
                "evaluated_count": denominator,
                "hit_count": hits,
                "precision_denominator": denominator,
                "eligible_reference_count": eligible_count,
                "all_reference_count": all_count,
                "precision": hits / denominator if denominator else np.nan,
                "eligible_reference_recall": hits / eligible_count if eligible_count else np.nan,
                "all_reference_coverage": hits / all_count if all_count else np.nan,
                "primary_identity": "stable_id",
                "ranking_sha256": actual_hash,
            }
        )
    pd.DataFrame(rows).to_csv(args.output, index=False)

    reference_exact_names = set(reference["path.2.name"].dropna().astype(str))
    reference_normalized_names = {_normalize_label(value) for value in reference_exact_names}
    audit = ranking.copy()
    audit["stable_id_hit"] = audit["path.2.id"].astype(str).isin(reference_ids)
    audit["exact_label_hit_diagnostic"] = audit["path.2.name"].astype(str).isin(reference_exact_names)
    audit["normalized_label_hit_diagnostic"] = audit["path.2.name"].map(_normalize_label).isin(reference_normalized_names)
    audit.to_csv(args.identity_audit, index=False)
    print(f"Wrote {args.output} and {args.identity_audit}; primary matches use stable IDs only")


def cmd_filter_enrichment(args):
    frame = pd.read_csv(args.input)
    _require_columns(frame, ["path.2.enrich_pvalue", "path.2.id", "path.2.name"], "Enrichment table")
    pvalues = pd.to_numeric(frame["path.2.enrich_pvalue"], errors="coerce")
    mask = pvalues.lt(args.pvalue_cutoff)
    if "path.2.enrichment_status" in frame:
        mask &= frame["path.2.enrichment_status"].astype(str).eq("exact")
    result = frame.loc[mask].copy().sort_values(
        ["path.2.enrich_pvalue", "path.2.id"], kind="mergesort"
    )
    if result.empty and not args.allow_empty:
        raise ValueError("No exact KG-wide enrichment rows pass the fixed cutoff")
    result.to_csv(args.output, index=False)
    print(f"Wrote {args.output} ({len(result)} strict rows)")


def _canonical_graph_hash(graph):
    node_rows = sorted((str(node), str(data.get("type", ""))) for node, data in graph.nodes(data=True))
    edge_rows = sorted(tuple(sorted((str(left), str(right)))) for left, right in graph.edges())
    payload = json.dumps({"nodes": node_rows, "edges": edge_rows}, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _coerce_identity_list(value):
    """Parse scalar or stringified-list identity fields without reordering."""
    if isinstance(value, list):
        return [str(item) for item in value]
    if isinstance(value, (tuple, set)):
        return [str(item) for item in value]
    if pd.isna(value):
        return []
    if isinstance(value, str):
        text = value.strip()
        if text.startswith("[") and text.endswith("]"):
            try:
                parsed = ast.literal_eval(text)
            except (SyntaxError, ValueError):
                parsed = None
            if isinstance(parsed, (list, tuple, set)):
                return [str(item) for item in parsed]
        return [text]
    return [str(value)]


def _source_identity_tuples(row):
    """Return aligned (name, type, stable_id) tuples for aggregate sources."""
    names = _coerce_identity_list(row["path.0.name"])
    stable_ids = _coerce_identity_list(row["path.0.id"])
    entity_types = _coerce_identity_list(row["path.0.type"])
    if len(entity_types) == 1:
        entity_types = entity_types * len(names)
    if len(names) == len(stable_ids) == len(entity_types):
        return list(zip(names, entity_types, stable_ids))
    if len(names) == 1 and len(entity_types) == 1 and stable_ids:
        return [(names[0], entity_types[0], stable_id) for stable_id in stable_ids]
    if len(stable_ids) == 1 and len(entity_types) == len(names):
        return [(name, entity_type, stable_ids[0]) for name, entity_type in zip(names, entity_types)]
    raise ValueError(
        "Aggregate source names, types, and stable IDs cannot be aligned without guessing: "
        f"{len(names)} names, {len(entity_types)} types, {len(stable_ids)} IDs"
    )


def _qualify_target_labels_for_native_graph(enrichment, omics_node_labels):
    """Make label-keyed target nodes collision-safe while preserving display labels."""
    frame = enrichment.copy()
    for column in ("path.2.id", "path.2.name", "path.2.type"):
        if frame[column].isna().any() or frame[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Graph input contains a missing target identity field: {column}")
    frame["__display_path2_name"] = frame["path.2.name"].astype(str)

    targets = (
        frame[["path.2.type", "path.2.id", "__display_path2_name"]]
        .astype(str)
        .drop_duplicates()
        .sort_values(["__display_path2_name", "path.2.type", "path.2.id"], kind="mergesort")
    )
    source_names = {
        name
        for _, row in frame.iterrows()
        for name, _entity_type, _stable_id in _source_identity_tuples(row)
    }
    reserved_native_labels = set(map(str, omics_node_labels)) | source_names
    target_identity_counts = targets.groupby("__display_path2_name", dropna=False).size()
    colliding_display_labels = {
        str(label) for label, count in target_identity_counts.items() if int(count) > 1
    }
    colliding_display_labels |= set(targets["__display_path2_name"]) & reserved_native_labels

    gene_targets = targets.loc[targets["path.2.type"].eq("Gene")]
    ambiguous_gene_labels = {
        str(label)
        for label, group in gene_targets.groupby("__display_path2_name", dropna=False)
        if group["path.2.id"].nunique(dropna=False) > 1
    }
    if ambiguous_gene_labels:
        examples = ", ".join(sorted(ambiguous_gene_labels)[:10])
        raise ValueError(
            "Gene target display labels map to multiple stable IDs and cannot be split "
            f"without breaking omics integration: {examples}"
        )

    target_mappings = []
    occupied_labels = reserved_native_labels | (
        set(targets["__display_path2_name"]) - colliding_display_labels
    )
    for row in targets.itertuples(index=False, name=None):
        entity_type, stable_id, display_label = map(str, row)
        needs_qualification = entity_type != "Gene" and display_label in colliding_display_labels
        graph_label = display_label
        if needs_qualification:
            graph_label = f"{display_label} [{stable_id}]"
            if graph_label in occupied_labels:
                graph_label = f"{display_label} [{entity_type}:{stable_id}]"
            if graph_label in occupied_labels:
                raise ValueError(
                    "Stable-ID target qualification still collides with an existing graph label: "
                    f"{graph_label}"
                )
            occupied_labels.add(graph_label)
        target_mappings.append(
            {
                "entity_type": entity_type,
                "display_label": display_label,
                "stable_id": stable_id,
                "graph_label": graph_label,
                "qualified": bool(needs_qualification),
            }
        )

    mapping = {
        (row["entity_type"], row["stable_id"], row["display_label"]): row["graph_label"]
        for row in target_mappings
    }
    qualified_input_rows = 0
    for key, graph_label in mapping.items():
        entity_type, stable_id, display_label = key
        if graph_label == display_label:
            continue
        mask = (
            frame["path.2.type"].astype(str).eq(entity_type)
            & frame["path.2.id"].astype(str).eq(stable_id)
            & frame["__display_path2_name"].astype(str).eq(display_label)
        )
        qualified_input_rows += int(mask.sum())
        frame.loc[mask, "path.2.name"] = graph_label

    qualified_entities = [row for row in target_mappings if row["qualified"]]
    audit = {
        "policy": (
            "use '<display label> [<stable ID>]' only as the native graph key for "
            "non-Gene target labels that would otherwise merge distinct entities or omics nodes"
        ),
        "display_labels_changed": False,
        "selection_or_filtering_used": False,
        "colliding_display_label_count": len(colliding_display_labels),
        "colliding_display_labels": sorted(colliding_display_labels),
        "qualified_entity_count": len(qualified_entities),
        "qualified_input_row_count": qualified_input_rows,
        "qualified_entities": qualified_entities,
    }
    return frame, audit


def cmd_build_graph(args):
    import scanpy as sc
    from embedgraph import combine_graph

    enrichment = pd.read_csv(args.input)
    _require_columns(
        enrichment,
        ["path.0.id", "path.0.name", "path.0.type", "path.2.id", "path.2.name", "path.2.type"],
        "Strict enrichment table",
    )
    if enrichment.empty:
        raise ValueError("Strict enrichment table is empty; refusing an all-edge fallback")
    if "path.2.enrich_pvalue" in enrichment:
        pvalues = pd.to_numeric(enrichment["path.2.enrich_pvalue"], errors="coerce")
        if pvalues.isna().any() or pvalues.ge(args.pvalue_cutoff).any():
            raise ValueError("Graph input contains non-exact or above-cutoff enrichment rows")

    genes = [line.strip() for line in Path(args.genes).read_text(encoding="utf-8").splitlines() if line.strip()]
    if not genes:
        raise ValueError("Gene list is empty")
    adata = sc.read_h5ad(args.adata)
    selected = sorted(set(genes) & set(adata.var_names.astype(str)))
    if not selected:
        raise ValueError("None of the selected genes occur in AnnData.var_names")
    omics_node_labels = set(selected) | set(adata.obs_names.astype(str))
    enrichment, qualification_audit = _qualify_target_labels_for_native_graph(
        enrichment, omics_node_labels
    )

    anchor_ids = set(
        enrichment.loc[
            enrichment["path.2.type"].astype(str).eq("Disease")
            & enrichment["__display_path2_name"].astype(str).eq(str(args.disease)),
            "path.2.id",
        ].dropna().astype(str)
    )
    if len(anchor_ids) != 1:
        raise ValueError(
            f"Expected display label {args.disease!r} to resolve to exactly one Disease stable ID; "
            f"found {len(anchor_ids)}"
        )

    graph = combine_graph(adata[:, selected].copy(), enrichment)
    nx.write_graphml(graph, args.output)

    identity = {}
    display_labels = {}
    for _, row in enrichment.iterrows():
        for label, entity_type, stable_id in _source_identity_tuples(row):
            native_type = entity_type if entity_type == "Gene" else f"KG_{entity_type}"
            key = (label, native_type)
            identity.setdefault(key, set()).add(stable_id)
            display_labels.setdefault(key, set()).add(label)

        label = str(row["path.2.name"])
        entity_type = str(row["path.2.type"])
        native_type = entity_type if entity_type == "Gene" else f"KG_{entity_type}"
        key = (label, native_type)
        identity.setdefault(key, set()).add(str(row["path.2.id"]))
        display_labels.setdefault(key, set()).add(str(row["__display_path2_name"]))

    node_rows = []
    for index, (node, data) in enumerate(graph.nodes(data=True)):
        native_type = str(data.get("type", ""))
        stable_ids = sorted(identity.get((str(node), native_type), set()))
        node_display_labels = sorted(display_labels.get((str(node), native_type), {str(node)}))
        node_rows.append(
            {
                "node_order": index,
                "native_label": str(node),
                "display_label": node_display_labels[0] if len(node_display_labels) == 1 else str(node),
                "native_type": native_type,
                "stable_id": stable_ids[0] if len(stable_ids) == 1 else "",
                "stable_id_set": "|".join(stable_ids),
                "stable_id_count": len(stable_ids),
            }
        )
    nodes = pd.DataFrame(node_rows)
    ambiguous = nodes.loc[
        nodes["native_type"].isin(["KG_Chemical", "KG_Disease"])
        & nodes["stable_id_count"].ne(1)
    ]
    if not ambiguous.empty:
        examples = ", ".join(ambiguous["native_label"].head(10).astype(str))
        raise ValueError(f"Chemical/Disease graph nodes are not uniquely grounded: {examples}")
    disease = nodes.loc[
        nodes["native_type"].eq("KG_Disease") & nodes["display_label"].eq(args.disease)
    ]
    if len(disease) != 1:
        raise ValueError(f"Expected exactly one KG_Disease node named {args.disease!r}; found {len(disease)}")
    nodes.to_csv(args.identity_output, index=False)

    audit = {
        "created_at": _utc_now(),
        "graph": str(Path(args.output).resolve()),
        "graph_sha256": _sha256(args.output),
        "canonical_topology_sha256": _canonical_graph_hash(graph),
        "node_count": graph.number_of_nodes(),
        "edge_count": graph.number_of_edges(),
        "connected_components": nx.number_connected_components(graph.to_undirected()),
        "isolated_nodes": int(len(list(nx.isolates(graph)))),
        "selected_gene_count": len(selected),
        "strict_enrichment_rows": len(enrichment),
        "all_edge_fallback_used": False,
        "stable_id_label_qualification": qualification_audit,
        "native_graph_defaults": {
            "directed": False,
            "weighted": False,
            "cell_neighbor": False,
        },
    }
    _write_json(args.audit_output, audit)
    print(f"Wrote {args.output}, {args.identity_output}, and {args.audit_output}")


def _official_grace_model():
    """Import the third-party model without the former local-name collision."""
    from embedgraph import GRACE_model

    origin = Path(inspect.getfile(GRACE_model)).resolve()
    local_dir = Path(__file__).resolve().parent
    if origin.parent == local_dir or origin == local_dir / "brickembed.py":
        raise ImportError("The third-party embedgraph package is shadowed by the brick skill")
    return GRACE_model, origin


def _cosine(left, right):
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator == 0.0:
        return float("nan")
    return float(np.dot(left, right) / denominator)


def cmd_grace_rank(args):
    import torch

    if args.device != "cuda":
        raise ValueError("Formal GRACE drug ranking requires --device cuda")
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable; CPU fallback is forbidden")
    graph = nx.read_graphml(args.graph)
    nodes = pd.read_csv(args.identity)
    _require_columns(
        nodes,
        ["node_order", "native_label", "display_label", "native_type", "stable_id", "stable_id_count"],
        "Graph identity table",
    )
    graph_nodes = [str(node) for node in graph.nodes()]
    identity_nodes = nodes.sort_values("node_order")["native_label"].astype(str).tolist()
    if graph_nodes != identity_nodes:
        raise ValueError("Graph node order does not match the frozen identity table")
    candidates = nodes.loc[nodes["native_type"].eq("KG_Chemical")].copy()
    disease = nodes.loc[
        nodes["native_type"].eq("KG_Disease") & nodes["display_label"].eq(args.disease)
    ]
    if len(disease) != 1:
        raise ValueError(f"Expected exactly one KG_Disease node named {args.disease!r}; found {len(disease)}")
    if candidates.empty:
        raise ValueError("No KG_Chemical nodes occur in the graph")
    if candidates["stable_id_count"].astype(int).ne(1).any() or candidates["stable_id"].fillna("").astype(str).str.strip().eq("").any():
        raise ValueError("Every ranked Chemical must resolve to exactly one stable ID")

    GRACE_model, grace_origin = _official_grace_model()
    # PyG's from_networkx attempts to tensorize every node/edge attribute. KG
    # labels and native types are strings, so keep them in the frozen identity
    # table and pass an order-preserving topology-only graph to GRACE.
    primary_graph = graph.__class__()
    primary_graph.add_nodes_from(graph_nodes)
    primary_graph.add_edges_from((str(left), str(right)) for left, right in graph.edges())
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    seed_rank_rows = []
    loss_rows = []
    runtime_rows = []
    signature = inspect.signature(GRACE_model)
    resolved_defaults = {
        name: parameter.default
        for name, parameter in signature.parameters.items()
        if parameter.default is not inspect.Parameter.empty
    }
    node_hash = hashlib.sha256("\n".join(graph_nodes).encode("utf-8")).hexdigest()

    for seed in args.seeds:
        torch.cuda.reset_peak_memory_stats()
        started = time.perf_counter()
        model = GRACE_model(primary_graph, device="cuda", seed=seed)
        model.fit(epoches=args.epochs)
        elapsed = time.perf_counter() - started
        embedding = np.asarray(model.get_embbedding(), dtype=float)
        if embedding.shape[0] != len(graph_nodes):
            raise RuntimeError("GRACE embedding row count does not match graph node order")
        losses = list(model.losses)
        if len(losses) != args.epochs:
            raise RuntimeError(f"GRACE seed {seed} returned {len(losses)} losses, expected {args.epochs}")

        loss_frame = pd.DataFrame({"seed": seed, "epoch": range(1, args.epochs + 1), "loss": losses})
        loss_path = output_dir / f"grace_loss_seed_{seed}.csv"
        loss_frame.to_csv(loss_path, index=False)
        loss_rows.append(loss_frame)

        embedding_frame = pd.DataFrame(embedding, columns=[f"embedding_{index + 1}" for index in range(embedding.shape[1])])
        embedding_frame.insert(0, "native_label", graph_nodes)
        embedding_frame.insert(0, "node_order", range(len(graph_nodes)))
        embedding_path = output_dir / f"embedding_seed_{seed}_{embedding.shape[1]}d.csv"
        embedding_frame.to_csv(embedding_path, index=False)

        disease_index = int(disease.iloc[0]["node_order"])
        disease_vector = embedding[disease_index]
        rows = []
        for _, candidate in candidates.iterrows():
            node_index = int(candidate["node_order"])
            rows.append(
                {
                    "seed": int(seed),
                    "path.2.id": str(candidate["stable_id"]),
                    "path.2.name": str(candidate["display_label"]),
                    "path.2.type": "Chemical",
                    "native_type": "KG_Chemical",
                    "raw_cosine_similarity": _cosine(embedding[node_index], disease_vector),
                    "graph_degree_diagnostic": int(graph.degree[graph_nodes[node_index]]),
                }
            )
        seed_frame = pd.DataFrame(rows).sort_values(
            ["raw_cosine_similarity", "path.2.id", "path.2.name"],
            ascending=[False, True, True],
            kind="mergesort",
        )
        seed_frame["raw_rank"] = range(1, len(seed_frame) + 1)
        seed_rank_rows.append(seed_frame)
        runtime_rows.append(
            {
                "seed": int(seed),
                "elapsed_seconds": elapsed,
                "final_loss": float(losses[-1]),
                "embedding_sha256": _sha256(embedding_path),
                "loss_sha256": _sha256(loss_path),
                "peak_gpu_memory_bytes": int(torch.cuda.max_memory_allocated()),
            }
        )

    seed_ranks = pd.concat(seed_rank_rows, ignore_index=True)
    seed_ranks.to_csv(output_dir / "chemical_seed_ranks.csv", index=False)
    pd.concat(loss_rows, ignore_index=True).to_csv(output_dir / "primary_grace_loss.csv", index=False)

    consensus = seed_ranks.pivot_table(
        index=["path.2.id", "path.2.name", "path.2.type", "native_type"],
        columns="seed",
        values=["raw_cosine_similarity", "raw_rank"],
        aggfunc="first",
    )
    consensus.columns = [f"{metric}_seed_{seed}" for metric, seed in consensus.columns]
    consensus = consensus.reset_index()
    rank_columns = [f"raw_rank_seed_{seed}" for seed in args.seeds]
    score_columns = [f"raw_cosine_similarity_seed_{seed}" for seed in args.seeds]
    consensus["mean_raw_rank"] = consensus[rank_columns].mean(axis=1)
    consensus["median_raw_rank"] = consensus[rank_columns].median(axis=1)
    consensus["raw_rank_std"] = consensus[rank_columns].std(axis=1, ddof=0)
    consensus["mean_raw_cosine_similarity"] = consensus[score_columns].mean(axis=1)
    consensus = consensus.sort_values(
        ["mean_raw_rank", "path.2.id", "path.2.name"], kind="mergesort"
    ).reset_index(drop=True)
    consensus.insert(0, "consensus_rank", range(1, len(consensus) + 1))
    ranking_path = output_dir / "primary_chemical_ranking.csv"
    consensus.to_csv(ranking_path, index=False)
    freeze = freeze_ranking(ranking_path, output_dir / "ranking_frozen.sha256")

    runtime = {
        "created_at": _utc_now(),
        "method": "three_seed_mean_raw_rank_consensus",
        "validation_used_for_ranking": False,
        "seeds": list(args.seeds),
        "epochs": args.epochs,
        "device": "cuda",
        "cuda_used": True,
        "gpu": torch.cuda.get_device_name(torch.cuda.current_device()),
        "torch_version": torch.__version__,
        "python_version": platform.python_version(),
        "grace_model_origin": str(grace_origin),
        "grace_signature": str(signature),
        "resolved_constructor_defaults": resolved_defaults,
        "ordered_node_list_sha256": node_hash,
        "grace_graph_attributes": "topology_only; labels/types preserved in external identity table",
        "runs": runtime_rows,
        "freeze": freeze,
    }
    _write_json(output_dir / "primary_grace_runtime.json", runtime)
    print(f"Wrote frozen GRACE consensus ranking to {ranking_path}")


def build_parser():
    parser = argparse.ArgumentParser(prog="brick drug", description="Drug-discovery identity, scope, GRACE, and validation helpers")
    sub = parser.add_subparsers(dest="command")

    scope = sub.add_parser("scope", help="Tag Chemical entities without filtering the formal candidate set")
    scope.add_argument("-i", "--input", required=True)
    scope.add_argument("-o", "--output", required=True)
    scope.set_defaults(func=cmd_scope)

    extract = sub.add_parser("extract", help="Extract a stable-ID-first Chemical ranking")
    extract.add_argument("-i", "--input", required=True)
    extract.add_argument("--scope")
    extract.add_argument("--top-n", type=int, default=50)
    extract.add_argument("-o", "--output", required=True)
    extract.set_defaults(func=cmd_extract)

    freeze = sub.add_parser("freeze", help="Freeze a primary ranking before any validation query")
    freeze.add_argument("-i", "--input", required=True)
    freeze.add_argument("-o", "--output", required=True)
    freeze.set_defaults(func=cmd_freeze)

    validate = sub.add_parser("validate", help="Evaluate a frozen ranking using stable IDs")
    validate.add_argument("--ranking", required=True)
    validate.add_argument("--freeze", required=True)
    validate.add_argument("--reference", required=True)
    validate.add_argument("--cutoffs", type=lambda value: [int(x) for x in value.split(",")], default=[10, 25, 50, 100])
    validate.add_argument("--allow-preexisting-reference", action="store_true")
    validate.add_argument("--identity-audit", required=True)
    validate.add_argument("-o", "--output", required=True)
    validate.set_defaults(func=cmd_validate)

    filter_enrich = sub.add_parser("filter-enrichment", help="Keep exact KG-wide rows below a fixed raw-p cutoff")
    filter_enrich.add_argument("-i", "--input", required=True)
    filter_enrich.add_argument("--pvalue-cutoff", type=float, default=0.0001)
    filter_enrich.add_argument("--allow-empty", action="store_true")
    filter_enrich.add_argument("-o", "--output", required=True)
    filter_enrich.set_defaults(func=cmd_filter_enrichment)

    build_graph = sub.add_parser("build-graph", help="Build the native joint graph from strict enrichment")
    build_graph.add_argument("-a", "--adata", required=True)
    build_graph.add_argument("--genes", required=True)
    build_graph.add_argument("-i", "--input", required=True)
    build_graph.add_argument("--disease", default="COVID-19")
    build_graph.add_argument("--pvalue-cutoff", type=float, default=0.0001)
    build_graph.add_argument("--identity-output", required=True)
    build_graph.add_argument("--audit-output", required=True)
    build_graph.add_argument("-o", "--output", required=True)
    build_graph.set_defaults(func=cmd_build_graph)

    grace = sub.add_parser("grace-rank", help="Train fixed-seed GPU GRACE and freeze a raw-cosine consensus")
    grace.add_argument("-g", "--graph", required=True)
    grace.add_argument("--identity", required=True)
    grace.add_argument("--disease", default="COVID-19")
    grace.add_argument("--seeds", type=lambda value: [int(x) for x in value.split(",")], default=list(DEFAULT_GRACE_SEEDS))
    grace.add_argument("--epochs", type=int, default=1000)
    grace.add_argument("--device", choices=["cuda"], default="cuda")
    grace.add_argument("-o", "--output-dir", required=True)
    grace.set_defaults(func=cmd_grace_rank)
    return parser


def _cli_main():
    parser = build_parser()
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        return
    try:
        args.func(args)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    _cli_main()
