import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import re
import sys
import time

import pandas as pd
import requests
from scipy.stats import fisher_exact

from querygraph import (
    _drop_duplicate_paths,
    _normalize_entity_type_filter,
    _path_key_columns,
    _query_relation_single_batch,
    query_cypher,
    query_relation,
)

try:
    from kgpolicy import source_quality_weight
except Exception:  # pragma: no cover
    def source_quality_weight(info_source):
        return 0.55


def _benjamini_hochberg(pvalues):
    """Return Benjamini-Hochberg adjusted p-values in original order."""
    pvalues = pd.Series(pvalues, dtype=float)
    n = len(pvalues)
    if n == 0:
        return []

    ranked = pvalues.sort_values(ascending=False)
    adjusted = []
    running_min = 1.0
    for rank_from_low, (idx, pvalue) in enumerate(ranked.items(), start=1):
        rank = n - rank_from_low + 1
        running_min = min(running_min, pvalue * n / rank)
        adjusted.append((idx, min(running_min, 1.0)))

    return pd.Series(dict(adjusted)).loc[pvalues.index].tolist()


def _relation_column(df):
    """Return the relation column used by brick query output."""
    for col in ('path.1.relation', 'path.1.name', 'path.1'):
        if col in df.columns:
            return col
    raise ValueError(
        "Input DataFrame is missing a relation column. Expected one of "
        "'path.1.relation', 'path.1.name', or 'path.1'."
    )


def infer_terminal_node_prefix(df):
    """Infer the terminal node prefix for one-hop or multi-hop brick qr output."""
    indices = []
    for col in df.columns:
        match = re.match(r"^path\.(\d+)\.(?:name|id|type)$", str(col))
        if match:
            idx = int(match.group(1))
            if idx % 2 == 0 and idx > 0:
                values = df[col].dropna().astype(str)
                if not values.empty and values.str.strip().ne("").any():
                    indices.append(idx)
    return f"path.{max(indices)}" if indices else "path.2"


def _target_prefix_from_arg(df, target_index="auto"):
    if target_index in (None, "", "auto"):
        return infer_terminal_node_prefix(df)
    target_index = str(target_index).strip()
    if target_index.isdigit():
        target_index = f"path.{target_index}"
    if not re.match(r"^path\.[02468]$", target_index):
        raise ValueError("--target-index must be auto or an even path index such as path.2, path.4, or path.6")
    return target_index


def _apply_target_index(df, target_index="auto"):
    """Copy the selected terminal node columns onto path.2.* for legacy ranking code."""
    target_prefix = _target_prefix_from_arg(df, target_index)
    if target_prefix == "path.2":
        return df.copy(), target_prefix
    out = df.copy()
    for suffix in ("id", "name", "type", "info_source_count", "rank_voting"):
        src = f"{target_prefix}.{suffix}"
        dst = f"path.2.{suffix}"
        if src in out.columns:
            out[dst] = out[src]
    if "path.2.name" not in out.columns or out["path.2.name"].dropna().astype(str).str.strip().eq("").all():
        raise ValueError(f"Selected target index {target_prefix} does not contain path node names")
    return out, target_prefix


def _stable_identity_column(df, prefix="path.2"):
    """Return the stable-ID column when it is usable, otherwise the label column."""
    id_col = f"{prefix}.id"
    name_col = f"{prefix}.name"
    if id_col in df.columns:
        values = df[id_col].dropna().astype(str).str.strip()
        if not values.empty and values.ne("").all():
            return id_col
    if name_col in df.columns:
        return name_col
    raise ValueError(f"Input DataFrame is missing {id_col} and {name_col}")


def _query_df2count_df(df, group_col=None,
                       retained_col = ['path.0.name',
                                       'path.0.id',
                                       'path.1.relation', 'path.1.info_source', 'path.1.info_source_length', 'path.1.relation_confidence',
                                       'path.1.original_relation', 'path.1.condition', 'path.2.id', 'path.2.name', 'path.2.type']):
    """
    Groups the given DataFrame by a specified column and processes the retained columns, returning a new DataFrame.

    This function takes a DataFrame, groups it based on a specified column, and processes the retained columns
    as per the requirements. The result is a new DataFrame that contains the processed information.

    Args:
        df (pandas.DataFrame): The input DataFrame to be processed.
        group_col (str, default='path.2.name'): The column name to group by.
        retained_col (list of str, optional): A list of column names to retain and process. Default includes 'path.0.name',
                                              'path.1.relation', 'path.1.info_source_length',
                                              'path.1.relation_confidence', 'path.2.id', 'path.2.name', and 'path.2.type'.

    Returns:
        pandas.DataFrame: A new DataFrame with the processed results after grouping by the specified column
                          and retaining the specified columns.
    """
    group_col = group_col or _stable_identity_column(df)
    group_path = '.'.join(group_col.split('.')[:2])

    retained_col = [x for x in retained_col if x in df]
    df_count = defaultdict(list)
    for x, y in df.groupby(group_col):
        #df_count[group_col].append(x)
        for i in retained_col:
            tmp = list(y[i])
            if group_path in i:
                df_count[i].append(tmp[0])
            else:
                df_count[i].append(tmp)

    df_count = pd.DataFrame(df_count)
    return df_count

def _prune_parent(query_df):
    """
    Prune the parent nodes of the query DataFrame.

    Args:
        query_df (pandas.DataFrame): The input DataFrame to be processed.

    Returns:
        pandas.DataFrame: A new DataFrame with the processed results after pruning the parent nodes.
    """
    query_df = query_df.copy()
    all_possible_nodes = list(query_df['path.2.name'].unique())
    all_possible_types = "|".join(query_df['path.2.type'].unique())
    parent2children = query_relation(all_possible_nodes, target_entity_type=all_possible_types, relation='is_a', directed=True)
    parent_cells = [x for x, y in parent2children.groupby('path.2.name') if y.shape[0] == 1]
    only_child_parent = parent2children.loc[parent2children['path.2.name'].isin(parent_cells)]
    only_child_parent = dict(zip(only_child_parent['path.2.name'], only_child_parent['path.0.name']))
    query_df['path.2.name'] = [ only_child_parent[x] if x in only_child_parent else x for x in query_df['path.2.name'] ]
    return query_df

def _match_count(df, target_group_df=None,  prune_parent=False, **kwargs):
    """
    Groups the given DataFrame by specified columns, processes the retained columns, and returns a new DataFrame
    with statistical results.

    This function allows grouping the input DataFrame by the specified columns, processes the retained columns,
    and computes statistical results based on those columns. It can also integrate additional statistics from
    a pre-grouped target DataFrame when multiple sorting indicators are used.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be processed.
        target_group_df (pandas.DataFrame, optional): A pre-grouped DataFrame with 'target' values.
                                                      This is used when computing additional statistics
                                                      with multiple sorting indicators.
        group_col (str, optional): The column name to group by. The default value is 'path.2.name' and it can be
                                   passed via **kwargs.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        retained_col (list of str, optional): A list of column names to retain and process. This can be passed
                                               via **kwargs.

    Returns:
        pandas.DataFrame: A new DataFrame containing the processed statistics and calculated results.
    """
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df
    source_col = "path.0.id" if "path.0.id" in df_count else "path.0.name"
    df_count['path.2.match_count'] = [len(set(x)) for x in df_count[source_col]]
    return df_count.sort_values(by='path.2.match_count', ascending=False)

def _info_source_count(df, target_group_df=None, prune_parent=False, **kwargs):
    """
    Groups the given DataFrame by specified columns, processes the retained columns, and returns a new DataFrame
    with statistical results, including the count of information sources.

    This function allows grouping the input DataFrame by specified columns, processes the retained columns,
    and computes the total count of information sources based on the 'path.1.info_source_length' column.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be processed.
        target_group_df (pandas.DataFrame, optional): A pre-grouped DataFrame with 'target' values.
                                                      This is used when computing additional statistics
                                                      with multiple sorting indicators.
        group_col (str, optional): The column name to group by. The default value is 'path.2.name' and it can be
                                   passed via **kwargs.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        retained_col (list of str, optional): A list of column names to retain and process. This can be passed
                                               via **kwargs.

    Returns:
        pandas.DataFrame: A new DataFrame containing the processed statistics and calculated results,
                           including the information source count.
    """
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df
    df_count['path.2.info_source_count'] = [sum(x) for x in df_count['path.1.info_source_length']]
    return df_count.sort_values(by='path.2.info_source_count', ascending=False)

def _match_probability(df, target_group_df=None, prune_parent=False, **kwargs):
    """
    Groups the given DataFrame by specified columns, processes the retained columns, and returns a new DataFrame
    with statistical results, including the calculated match probability.

    This function groups the input DataFrame based on specified columns, processes the retained columns, and
    calculates the match probability for each group, providing a new DataFrame with these results.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be processed.
        target_group_df (pandas.DataFrame, optional): A pre-grouped DataFrame with 'target' values.
                                                      This is used when calculating additional statistics
                                                      with multiple sorting indicators.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        group_col (str, optional): The column name to group by. The default value is 'path.2.name' and it can
                                   be passed via **kwargs.
        retained_col (list of str, optional): A list of column names to retain and process. This can be passed
                                               via **kwargs.

    Returns:
        pandas.DataFrame: A new DataFrame containing the processed statistics and match probability results.
    """
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df

    source_col = _stable_identity_column(df, "path.0")
    target_col = _stable_identity_column(df, "path.2")
    adj = defaultdict(dict)
    for (a, b), y in df.groupby([source_col, target_col]):
        adj[a][b] = 1

    adj = pd.DataFrame(adj)
    adj.fillna(0, inplace=True)
    adj = adj / adj.sum()
    target2prob = dict(adj.mean(axis=1))
    output_target_col = target_col if target_col in df_count else "path.2.name"
    df_count['path.2.match_probability'] = [target2prob[x] for x in df_count[output_target_col]]
    return df_count.sort_values(by='path.2.match_probability', ascending=False)

def _info_source_probability(df, target_group_df=None, prune_parent=False, **kwargs):
    """
    Groups the given DataFrame by specified columns, processes the retained columns, and returns a new DataFrame
    with statistical results, including the calculated information source probability.

    This function groups the input DataFrame based on specified columns, processes the retained columns, and
    calculates the information source probability for each group, providing a new DataFrame with these results.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be processed.
        target_group_df (pandas.DataFrame, optional): A pre-grouped DataFrame with 'target' values.
                                                      This is used when calculating additional statistics
                                                      with multiple sorting indicators.
        group_col (str, optional): The column name to group by. The default value is 'path.2.name' and it can
                                   be passed via **kwargs.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        retained_col (list of str, optional): A list of column names to retain and process. This can be passed
                                               via **kwargs.

    Returns:
        pandas.DataFrame: A new DataFrame containing the processed statistics and information source probability results.
    """
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df

    source_col = _stable_identity_column(df, "path.0")
    target_col = _stable_identity_column(df, "path.2")
    adj = defaultdict(dict)
    for (a, b), y in df.groupby([source_col, target_col]):
        adj[a][b] = y['path.1.info_source_length'].sum()

    adj = pd.DataFrame(adj)
    adj.fillna(0, inplace=True)
    adj = adj / adj.sum()
    target2prob = dict(adj.mean(axis=1))
    output_target_col = target_col if target_col in df_count else "path.2.name"
    df_count['path.2.info_source_probability'] = [target2prob[x] for x in df_count[output_target_col]]
    return df_count.sort_values(by='path.2.info_source_probability', ascending=False)


def _missing_cypher_endpoint(exc):
    """Return True only when HTTP explicitly reports no Cypher endpoint."""
    response = getattr(exc, "response", None)
    return getattr(response, "status_code", None) in {404, 405, 501}


def _positive_int(value, default, env_name):
    if value is None:
        value = os.getenv(env_name, default)
    try:
        value = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{env_name} must be a positive integer") from exc
    if value < 1:
        raise ValueError(f"{env_name} must be a positive integer")
    return value


def _nonnegative_int(value, default, env_name):
    if value is None:
        value = os.getenv(env_name, default)
    try:
        value = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{env_name} must be a non-negative integer") from exc
    if value < 0:
        raise ValueError(f"{env_name} must be a non-negative integer")
    return value


def _retryable_relation_exception(exc):
    """Retry only transient transport or server-capacity failures."""
    if isinstance(exc, (requests.Timeout, requests.ConnectionError)):
        return True
    response = getattr(exc, "response", None)
    return getattr(response, "status_code", None) in {
        408, 425, 429, 500, 502, 503, 504,
    }


def _deduplicated_background_rows(frame):
    """Count exact returned paths per source target after transport deduplication."""
    if frame.empty:
        return {}, frame
    frame = _drop_duplicate_paths(frame, _path_key_columns(frame))
    if "path.0.id" not in frame.columns:
        raise RuntimeError("KG relation response is missing path.0.id")
    frame = frame.copy()
    frame["path.0.id"] = frame["path.0.id"].astype(str)
    counts = frame.groupby("path.0.id", dropna=False).size().astype(int).to_dict()
    return counts, frame


def _background_failure_threshold(match_count, source_set_length, pvalue_cutoff):
    """Return the first KG background count that cannot pass the fixed cutoff.

    The enrichment table used by this module is monotone in the background
    count for fixed ``match_count`` and ``source_set_length``.  Finding this
    threshold before the live request lets the relation endpoint stop after it
    has returned enough independent KG-wide evidence to prove that a target
    cannot enter the strict graph.  No query-subset rows are used.
    """
    if pvalue_cutoff is None:
        return None
    match_count = int(match_count)
    source_set_length = int(source_set_length)
    cutoff = float(pvalue_cutoff)
    low = 0
    high = max(1, match_count)
    while fisher_exact(
        [[match_count, high], [source_set_length, 20000]],
        alternative="greater",
    ).pvalue < cutoff:
        low = high
        high *= 2
        if high > 2_000_000:
            raise RuntimeError(
                "Could not bound the enrichment failure threshold below 2,000,000 rows"
            )
    while low + 1 < high:
        middle = (low + high) // 2
        pvalue = fisher_exact(
            [[match_count, middle], [source_set_length, 20000]],
            alternative="greater",
        ).pvalue
        if pvalue >= cutoff:
            high = middle
        else:
            low = middle
    return int(high)


def _relation_background_counts(
    df,
    targets,
    relation_col,
    batch_size=None,
    max_workers=None,
    max_retries=None,
    pvalue_cutoff=None,
    source_set_length=None,
    server_row_cap=None,
):
    """Obtain exact KG-wide target backgrounds through reverse relation queries.

    This is a server-backed reconstruction of the legacy ``count(m)`` query,
    not an estimate from the submitted source subset. Only target identity,
    source/target types, and relation labels already present in the live KG
    response are used.
    """
    required = {
        "path.2.id",
        "path.2.name",
        "path.2.type",
        "path.2.match_count",
    }
    missing = required.difference(targets.columns)
    if missing:
        raise ValueError(
            "Relation-backed enrichment requires target identity columns: "
            + ", ".join(sorted(missing))
        )

    source_types = _normalize_entity_type_filter(
        [
            str(value).strip()
            for value in df["path.0.type"].dropna().unique()
            if str(value).strip()
        ]
    )
    if not source_types:
        raise ValueError("Cannot compute a KG background without a source entity type")
    target_relation_lookup = {}
    for (target_type, target_id), group in df.groupby(
        ["path.2.type", "path.2.id"],
        sort=False,
        dropna=False,
    ):
        relation_signature = tuple(
            sorted(
                {
                    str(value).strip()
                    for value in group[relation_col].dropna()
                    if str(value).strip()
                }
            )
        )
        if relation_signature:
            target_relation_lookup[(str(target_type).strip(), str(target_id).strip())] = relation_signature
    if not target_relation_lookup:
        raise ValueError("Cannot compute a KG background without relationship labels")

    targets = targets[
        ["path.2.id", "path.2.name", "path.2.type", "path.2.match_count"]
    ].copy()
    for column in ("path.2.id", "path.2.name", "path.2.type"):
        targets[column] = targets[column].astype(str).str.strip()
    if targets[["path.2.id", "path.2.name", "path.2.type"]].eq("").any(axis=None):
        raise ValueError("Relation-backed enrichment requires non-empty target IDs, names, and types")
    targets = targets.drop_duplicates(subset=["path.2.id", "path.2.type"])
    targets["_relation_signature"] = [
        target_relation_lookup.get((target_type, target_id), tuple())
        for target_type, target_id in zip(
            targets["path.2.type"],
            targets["path.2.id"],
        )
    ]
    missing_relations = targets.loc[
        targets["_relation_signature"].map(len).eq(0),
        ["path.2.type", "path.2.id"],
    ]
    if not missing_relations.empty:
        raise RuntimeError(
            "Primary KG evidence is missing relation labels for evaluated targets: "
            + "; ".join(
                f"{row['path.2.type']}:{row['path.2.id']}"
                for _, row in missing_relations.head(10).iterrows()
            )
        )

    if targets.empty:
        return {}, {}, {}

    # A bounded KG-wide request must be target-local: a shared server LIMIT
    # across several targets could be consumed by one hub and silently hide the
    # other targets.  ``batch_size`` remains accepted for API compatibility and
    # provenance, but the exact/censored implementation intentionally submits
    # one stable target per request.
    _positive_int(batch_size, 10, "BRICK_ENRICH_RELATION_BATCH_SIZE")
    max_workers = _positive_int(
        max_workers,
        4,
        "BRICK_ENRICH_RELATION_WORKERS",
    )
    max_retries = _nonnegative_int(
        max_retries,
        2,
        "BRICK_ENRICH_RELATION_RETRIES",
    )
    server_row_cap = _positive_int(
        server_row_cap,
        10000,
        "BRICK_ENRICH_RELATION_SERVER_ROW_CAP",
    )
    if pvalue_cutoff is not None and source_set_length is None:
        raise ValueError("source_set_length is required for cutoff-aware KG-wide screening")

    jobs = [
        (
            str(record["path.2.type"]),
            tuple(record["_relation_signature"]),
            record,
        )
        for record in targets.sort_values(
            ["path.2.type", "path.2.id"], kind="mergesort"
        ).to_dict("records")
    ]

    def query_with_transient_retries(**query_kwargs):
        for attempt in range(max_retries + 1):
            try:
                return _query_relation_single_batch(**query_kwargs), attempt
            except Exception as exc:
                if attempt >= max_retries or not _retryable_relation_exception(exc):
                    raise
                time.sleep(min(2 ** attempt, 4))
        raise AssertionError("unreachable")

    def fetch(target_type, relation_types, record):
        target_id = record["path.2.id"]
        key = (target_type, target_id)
        match_count = int(record["path.2.match_count"])
        failure_threshold = _background_failure_threshold(
            match_count,
            source_set_length,
            pvalue_cutoff,
        )
        requested_limit = (
            server_row_cap + 1
            if failure_threshold is None
            else min(int(failure_threshold), server_row_cap + 1)
        )
        method_parts = ["relation_id_target_relations_kg_wide"]

        while True:
            frame, id_transport_retries = query_with_transient_retries(
                source_entity_set=[target_id],
                target_entity_set=None,
                relation=relation_types,
                multi_hop=1,
                directed=False,
                source_entity_type=target_type,
                target_entity_type=source_types,
                relation_info_source=None,
                query_attribution="id",
                limit_number=requested_limit,
                return_type="dataframe",
            )
            raw_count = int(len(frame))
            counts_by_id, deduplicated = _deduplicated_background_rows(frame)
            background_count = int(counts_by_id.get(target_id, 0))
            if id_transport_retries:
                method_parts.append(f"transport_retry_{id_transport_retries}")

            if background_count < match_count and raw_count < requested_limit:
                retry, label_transport_retries = query_with_transient_retries(
                    source_entity_set=[record["path.2.name"]],
                    target_entity_set=None,
                    relation=relation_types,
                    multi_hop=1,
                    directed=False,
                    source_entity_type=target_type,
                    target_entity_type=source_types,
                    relation_info_source=None,
                    query_attribution="name",
                    limit_number=requested_limit,
                    return_type="dataframe",
                )
                retry_raw_count = int(len(retry))
                retry_counts, retry_deduplicated = _deduplicated_background_rows(retry)
                if target_id in retry_counts:
                    background_count = int(retry_counts[target_id])
                elif not retry_deduplicated.empty:
                    source_ids = retry_deduplicated["path.0.id"].dropna().astype(str).unique()
                    source_names = retry_deduplicated.get(
                        "path.0.name", pd.Series(dtype=str)
                    )
                    exact_name = source_names.astype(str).str.casefold().eq(
                        record["path.2.name"].casefold()
                    )
                    if len(source_ids) == 1 and exact_name.any():
                        background_count = int(len(retry_deduplicated))
                raw_count = retry_raw_count
                method_parts[0] = "relation_label_retry_target_relations_kg_wide"
                if label_transport_retries:
                    method_parts.append(
                        f"label_transport_retry_{label_transport_retries}"
                    )

            if failure_threshold is not None and background_count >= failure_threshold:
                method_parts.append(f"failure_bound_{failure_threshold}")
                return key, background_count, "_".join(method_parts), (
                    "screened_kg_wide_failure_bound"
                )

            # The live service currently caps relation results at 10,000 rows
            # even when a larger limit is requested.  Treat a cap hit as a
            # structural hub exclusion, never as an exact count.
            if raw_count >= server_row_cap:
                method_parts.append(f"server_cap_{server_row_cap}")
                return key, background_count, "_".join(method_parts), (
                    "excluded_kg_wide_saturated_hub"
                )

            if raw_count < requested_limit:
                if background_count < match_count:
                    raise RuntimeError(
                        f"KG relation background is incomplete for {key}: "
                        f"background={background_count}, match={match_count}"
                    )
                method_parts.append(f"complete_below_limit_{requested_limit}")
                return key, background_count, "_".join(method_parts), "exact"

            # The server filled the requested window but exact transport
            # duplicates kept the unique path count below the decision bound.
            # Grow the window deterministically without changing query scope.
            next_limit = min(server_row_cap + 1, max(requested_limit + 1, requested_limit * 2))
            if next_limit == requested_limit:
                raise RuntimeError(
                    f"Could not resolve KG-wide background under the server cap for {key}"
                )
            requested_limit = next_limit
            method_parts.append(f"expanded_limit_{requested_limit}")

    background_counts = {}
    identity_methods = {}
    background_statuses = {}
    completed = 0
    executor = ThreadPoolExecutor(max_workers=min(max_workers, len(jobs)))
    try:
        futures = {
            executor.submit(fetch, target_type, relation_types, record): (
                target_type,
                relation_types,
                record,
            )
            for target_type, relation_types, record in jobs
        }
        for future in as_completed(futures):
            key, count, method, status = future.result()
            if key in background_counts:
                raise RuntimeError(f"Duplicate KG background result for {key}")
            background_counts[key] = count
            identity_methods[key] = method
            background_statuses[key] = status
            completed += 1
            if completed == 1 or completed % 250 == 0 or completed == len(jobs):
                print(f"KG-wide enrichment targets: {completed}/{len(jobs)}")
    finally:
        # If one exact request hard-fails, cancel jobs that have not started so
        # the caller receives the real error instead of silently draining the
        # entire queue. Already-running requests are allowed to finish.
        executor.shutdown(wait=True, cancel_futures=True)

    unresolved = []
    for record in targets.to_dict("records"):
        key = (record["path.2.type"], record["path.2.id"])
        background = int(background_counts.get(key, 0))
        match_count = int(record["path.2.match_count"])
        if background < match_count:
            unresolved.append(
                f"{key}: background={background}, match={match_count}"
            )
    if unresolved:
        raise RuntimeError(
            "KG relation background is incomplete; refusing approximate enrichment. "
            + "; ".join(unresolved[:10])
        )
    return background_counts, identity_methods, background_statuses


def enrich(df, target_group_df=None, prune_parent=False, **kwargs):
    """
    Groups the given DataFrame by specified columns, processes the retained columns,
    and returns a new DataFrame with statistical results.

    This function enhances the input DataFrame by grouping it based on specified columns
    and processing the retained columns. If a pre-grouped DataFrame is provided,
    new computed indicators are added to it.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be processed.
        target_group_df (pandas.DataFrame, optional): A pre-grouped DataFrame with 'target' values.
                                                      Used when additional sorting indicators need
                                                      to be added to the existing results.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        group_col (str, optional): The column name to group by. The default is 'path.2.name',
                                   and it can be passed via **kwargs.
        retained_col (list of str, optional): A list of column names to retain and process.
                                              This can be passed via **kwargs.

    Returns:
        pandas.DataFrame: A new DataFrame containing the processed statistical information
                          and computed results.
    """
    target_index = kwargs.pop("target_index", "auto")
    background_method = kwargs.pop("background_method", "auto")
    pvalue_cutoff = kwargs.pop("pvalue_cutoff", None)
    relation_batch_size = kwargs.pop("relation_batch_size", None)
    relation_workers = kwargs.pop("relation_workers", None)
    relation_retries = kwargs.pop("relation_retries", None)
    relation_server_row_cap = kwargs.pop("relation_server_row_cap", None)
    if background_method not in {
        "auto",
        "cypher",
        "relation",
        "query_subset_degraded",
    }:
        raise ValueError(
            "background_method must be auto, cypher, relation, or "
            "query_subset_degraded"
        )
    if pvalue_cutoff is not None:
        pvalue_cutoff = float(pvalue_cutoff)
        if not 0.0 < pvalue_cutoff <= 1.0:
            raise ValueError("pvalue_cutoff must be within (0, 1]")
    df, target_prefix = _apply_target_index(df, target_index)
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        # Stable IDs are the lossless default identity for enrichment. Labels
        # are not guaranteed to be unique across KG namespaces or entities.
        if "group_col" not in kwargs and "path.2.id" in df.columns:
            kwargs["group_col"] = "path.2.id"
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df
    relation_col = _relation_column(df)
    source_type = ':' + '|'.join(list(df['path.0.type'].dropna().unique()))
    relation_type = ':' + '|'.join(list(df[relation_col].dropna().unique()))
    target_type = ':' + '|'.join(list(df['path.2.type'].dropna().unique()))

    source_set_length = df['path.0.id'].nunique()
    df_count = _match_count(df, target_group_df=df_count, **kwargs)
    df_count["path.2.enrich_pvalue_lower_bound"] = [
        fisher_exact(
            [[int(match_count), 0], [source_set_length, 20000]],
            alternative="greater",
        ).pvalue
        for match_count in df_count["path.2.match_count"]
    ]

    target_set = list(df_count['path.2.id'])
    targetid2bg_count = None
    degraded_counts = None
    if background_method == "query_subset_degraded":
        deduplicated = _drop_duplicate_paths(df, _path_key_columns(df))
        degraded_counts = (
            deduplicated.assign(
                **{
                    "path.2.id": deduplicated["path.2.id"].astype(str),
                    "path.2.type": deduplicated["path.2.type"].astype(str),
                }
            )
            .groupby(["path.2.type", "path.2.id"], dropna=False)
            .size()
            .astype(int)
            .to_dict()
        )
    elif background_method in {"auto", "cypher"}:
        try:
            cypher_counts = query_cypher(
                f"MATCH (n{target_type})-[r{relation_type}]-(m{source_type}) "
                "WHERE n.id IN $target_set return n.id, count(m)",
                parameters={"target_set": target_set},
            )
            targetid2bg_count = {
                str(target_id): int(count)
                for target_id, count in zip(
                    cypher_counts["n.id"], cypher_counts["count(m)"]
                )
            }
        except Exception as exc:
            if background_method == "cypher" or not _missing_cypher_endpoint(exc):
                raise
            print(
                "Cypher endpoint is unavailable; using exact reverse relation "
                "queries for KG background counts."
            )

    background_values = pd.Series(pd.NA, index=df_count.index, dtype="Int64")
    background_sources = pd.Series("", index=df_count.index, dtype="object")
    enrichment_status = pd.Series("exact", index=df_count.index, dtype="object")
    evaluated_mask = pd.Series(True, index=df_count.index, dtype=bool)

    if degraded_counts is not None:
        enrichment_status[:] = "degraded_query_subset_not_kg_wide"
        for index, row in df_count.iterrows():
            key = (str(row["path.2.type"]), str(row["path.2.id"]))
            background_values.loc[index] = int(degraded_counts.get(key, 0))
            background_sources.loc[index] = "query_subset_degraded_not_kg_wide"
    elif targetid2bg_count is not None:
        for index, row in df_count.iterrows():
            target_id = str(row["path.2.id"])
            background_values.loc[index] = int(targetid2bg_count.get(target_id, 0))
            background_sources.loc[index] = "cypher"
    else:
        if pvalue_cutoff is not None:
            evaluated_mask = (
                df_count["path.2.enrich_pvalue_lower_bound"] < pvalue_cutoff
            )
            enrichment_status.loc[~evaluated_mask] = "screened_exact_bound"
            background_sources.loc[~evaluated_mask] = "not_queried_exact_bound"

        evaluated_targets = df_count.loc[evaluated_mask]
        relation_counts, identity_methods, relation_statuses = _relation_background_counts(
            df,
            evaluated_targets,
            relation_col,
            batch_size=relation_batch_size,
            max_workers=relation_workers,
            max_retries=relation_retries,
            pvalue_cutoff=pvalue_cutoff,
            source_set_length=source_set_length,
            server_row_cap=relation_server_row_cap,
        )
        for index, row in evaluated_targets.iterrows():
            key = (str(row["path.2.type"]), str(row["path.2.id"]))
            background_values.loc[index] = int(relation_counts.get(key, 0))
            background_sources.loc[index] = identity_methods.get(key, "relation_id")
            enrichment_status.loc[index] = relation_statuses.get(key, "exact")

    df_count["path.2.background_count"] = background_values
    df_count["path.2.background_method"] = background_sources
    df_count["path.2.enrichment_status"] = enrichment_status

    incomplete = df_count.loc[
        evaluated_mask
        & (
            df_count["path.2.background_count"].isna()
            | (
                df_count["path.2.background_count"].astype("Int64")
                < df_count["path.2.match_count"].astype("Int64")
            )
        )
    ]
    if not incomplete.empty:
        raise RuntimeError(
            "KG background counts are incomplete or smaller than observed matches; "
            "refusing approximate enrichment"
        )

    enrich_statistics = pd.Series(float("nan"), index=df_count.index, dtype=float)
    enrich_pvalues = pd.Series(float("nan"), index=df_count.index, dtype=float)
    statistical_mask = evaluated_mask & enrichment_status.isin(
        ["exact", "degraded_query_subset_not_kg_wide"]
    )
    for index, row in df_count.loc[statistical_mask].iterrows():
        result = fisher_exact(
            [
                [int(row["path.2.match_count"]), int(row["path.2.background_count"])],
                [source_set_length, 20000],
            ],
            alternative="greater",
        )
        enrich_statistics.loc[index] = result.statistic
        enrich_pvalues.loc[index] = result.pvalue

    df_count["path.2.enrich_statistic"] = enrich_statistics
    df_count["path.2.enrich_pvalue"] = enrich_pvalues
    if statistical_mask.all():
        df_count["path.2.enrich_fdr"] = _benjamini_hochberg(enrich_pvalues)
        df_count["path.2.enrich_fdr_status"] = (
            "degraded_query_subset_all_targets"
            if background_method == "query_subset_degraded"
            else "exact_all_targets"
        )
    else:
        df_count["path.2.enrich_fdr"] = float("nan")
        df_count["path.2.enrich_fdr_status"] = (
            "not_computed_after_kg_wide_bound_or_hub_screen"
        )

    return df_count.sort_values(
        by='path.2.enrich_pvalue', ascending=True, na_position="last"
    )



WEAK_RELATION_PATTERNS = ("co_mentioned", "co-mentioned", "co occurrence", "co-occurrence", "text_mined", "text-mined", "literature_mined", "abstract")
RELATION_QUALITY_WEIGHTS = {
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


def _bounded_score(value, default=0.5):
    try:
        if pd.isna(value):
            return default
        value = float(value)
        if 1.0 < value <= 100.0:
            value /= 100.0
        return max(0.0, min(1.0, value))
    except Exception:
        return default


def _is_weak_row(row):
    text = " ".join(str(row.get(col, "")) for col in [
        "path.1.relation", "path.1.original_relation", "path.1.condition", "path.1.info_source",
    ]).casefold()
    return any(pattern in text for pattern in WEAK_RELATION_PATTERNS)


def _row_evidence_quality(row):
    relation = str(row.get("path.1.relation", "")).strip()
    relation_weight = RELATION_QUALITY_WEIGHTS.get(relation, 0.5)
    if _is_weak_row(row):
        relation_weight = min(relation_weight, 0.35)
    confidence = _bounded_score(row.get("path.1.relation_confidence"), 0.5)
    source_weight = source_quality_weight(row.get("path.1.info_source"))
    return relation_weight * (0.5 + 0.5 * confidence) * source_weight


def _evidence_quality_count(df, target_group_df=None, prune_parent=False, **kwargs):
    """Aggregate source/relation-aware evidence quality, with weak co-mentioned cap."""
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df
    scores = df.apply(_row_evidence_quality, axis=1)
    target_col = _stable_identity_column(df, "path.2")
    output_target_col = target_col if target_col in df_count else "path.2.name"
    target2score = scores.groupby(df[target_col]).sum().to_dict()
    df_count["path.2.evidence_quality_count"] = [target2score.get(x, 0.0) for x in df_count[output_target_col]]
    return df_count.sort_values(by="path.2.evidence_quality_count", ascending=False)


def rank(df, target_group_df=None,
                metrics=['match_count', 'match_probability', 'info_source_count', 'info_source_probability', 'evidence_quality_count'],
                ascending=False,
                prune_parent=False,
                **kwargs):
    """
    Ranks data using multiple metrics and applies a voting mechanism to determine the final ranking score.

    This function sorts the given DataFrame based on multiple specified metrics, assigns ranking scores,
    and aggregates them to produce a final ranking score.

    Args:
        df (pandas.DataFrame): The input DataFrame containing data to be ranked.
        target_group_df (pandas.DataFrame, optional): A pre-processed grouped DataFrame.
                                                      If None (default), the ranking metrics are computed.
        metrics (list of str, optional): A list of metric names used for ranking.
                                         The default includes 'match_count', 'match_probability',
                                         'info_source_count', and 'info_source_probability'.
        ascending (bool or list of bool, optional): Sorting order. If a single boolean is provided,
                                                    all metrics follow the same sorting order.
                                                    If a list is given, each metric follows its respective order.
        prune_parent (bool, optional): Whether to prune the parent nodes. The default value is False.
        kwargs: Additional parameters passed to the `_query_df2count_df` function.

    Returns:
        pandas.DataFrame: A DataFrame with ranking scores, including a new column 'path.2.rank_voting'
                          representing the final ranking score.
    """
    target_index = kwargs.pop("target_index", "auto")
    df, target_prefix = _apply_target_index(df, target_index)
    if prune_parent:
        df = _prune_parent(df)
    if target_group_df is None:
        df_count = _query_df2count_df(df, **kwargs)
    else:
        df_count = target_group_df
    metric_functions = {
        'match_count': _match_count,
        'match_probability': _match_probability,
        'info_source_count': _info_source_count,
        'info_source_probability': _info_source_probability,
        'evidence_quality_count': _evidence_quality_count,
    }
    unsupported = [metric for metric in metrics if metric not in metric_functions]
    if unsupported:
        raise ValueError(
            f"Unsupported ranking metric(s): {unsupported}. "
            f"Supported metrics: {list(metric_functions)}"
        )

    for metric in metrics:
        df_count = metric_functions[metric](df, df_count)

    if isinstance(ascending, bool):
        ascending_list = [ascending] * len(metrics)
    elif isinstance(ascending, list):
        ascending_list = [bool(x) for x in ascending]
    else:
        raise ValueError('ascending must be bool type or a list of bool type')
    rank_list = []
    for met, asc in zip(metrics, ascending_list):
        if met == 'enrich':
            met = 'enrich_pvalue'
            asc = True
        rank_list.append(df_count['path.2.' + met].rank(ascending=asc))

    rank_list = pd.DataFrame(rank_list).T
    df_count['path.2.rank_voting'] = list(rank_list.mean(axis=1))
    return df_count.sort_values(by='path.2.rank_voting', ascending=True)


def _cli_main():
    """CLI entry point for brick rk command."""
    parser = argparse.ArgumentParser(
        description="brick rk - Rank or enrich DataFrame from brick qr query results",
        usage='brick rk -i <input.csv> [options]'
    )

    parser.add_argument("-i", "--input", required=True, help="Input CSV file (from brick qr)")
    parser.add_argument("--rank", action="store_true", help="Apply ranking analysis (default)")
    parser.add_argument("--enrich", action="store_true", help="Apply enrichment analysis")
    parser.add_argument("-o", "--output", help="Output CSV file")

    # rank() function parameters
    parser.add_argument("--metrics", help="Metrics for ranking, comma-separated (default: match_count,match_probability,info_source_count,info_source_probability,evidence_quality_count)")
    parser.add_argument("--ascending", help="Sort order: 'true', 'false', or comma-separated list for each metric")
    parser.add_argument("--prune-parent", action="store_true", help="Prune parent nodes before analysis")
    parser.add_argument("--target-index", default="auto", help="Target node prefix for ranking/enrichment: auto, path.2, path.4, or path.6")

    # enrich() function parameters
    parser.add_argument("--target-group-df", help="Path to pre-grouped DataFrame CSV (optional)")
    parser.add_argument(
        "--background-method",
        choices=["auto", "cypher", "relation", "query_subset_degraded"],
        default="auto",
        help=(
            "KG background provider. auto falls back to exact relation queries only "
            "when Cypher is unavailable. query_subset_degraded uses only deduplicated "
            "rows from the primary query and is explicitly not KG-wide."
        ),
    )
    parser.add_argument(
        "--pvalue-cutoff",
        type=float,
        help="Optional exact lower-bound screen for a predeclared raw p-value cutoff.",
    )
    parser.add_argument(
        "--relation-batch-size",
        type=int,
        help="Targets per exact reverse-relation background request (default: 10).",
    )
    parser.add_argument(
        "--relation-workers",
        type=int,
        help="Concurrent exact background requests (default: 4).",
    )
    parser.add_argument(
        "--relation-retries",
        type=int,
        help="Transient retries per exact background request (default: 2).",
    )

    args = parser.parse_args()

    if not args.rank and not args.enrich:
        args.rank = True  # default to ranking

    if args.rank and args.enrich:
        print("Error: Cannot use --rank and --enrich together", file=sys.stderr)
        sys.exit(1)

    try:
        # Load input DataFrame
        print(f"Loading: {args.input}")
        df = pd.read_csv(args.input)

        if df.empty:
            print("Error: Input file is empty", file=sys.stderr)
            sys.exit(1)

        print(f"Loaded {len(df)} rows")

        # Parse rank parameters
        rank_kwargs = {"prune_parent": args.prune_parent, "target_index": args.target_index}
        if args.metrics:
            rank_kwargs["metrics"] = [m.strip() for m in args.metrics.split(",")]
        if args.ascending:
            asc_parts = [a.strip().lower() for a in args.ascending.split(",")]
            if len(asc_parts) == 1:
                rank_kwargs["ascending"] = asc_parts[0] == "true"
            else:
                rank_kwargs["ascending"] = [a == "true" for a in asc_parts]

        # Load target_group_df if provided
        target_group_df = None
        if args.target_group_df:
            print(f"Loading target group df: {args.target_group_df}")
            target_group_df = pd.read_csv(args.target_group_df)

        # Apply analysis
        if args.enrich:
            print("Applying enrichment analysis...")
            rank_kwargs.update(
                {
                    "background_method": args.background_method,
                    "pvalue_cutoff": args.pvalue_cutoff,
                    "relation_batch_size": args.relation_batch_size,
                    "relation_workers": args.relation_workers,
                    "relation_retries": args.relation_retries,
                }
            )
            result = enrich(df, target_group_df=target_group_df, **rank_kwargs)
            key_col = 'path.2.enrich_pvalue'
            print("\nTop enriched (by p-value):")
        else:
            print("Applying ranking analysis...")
            result = rank(df, target_group_df=target_group_df, **rank_kwargs)
            key_col = 'path.2.rank_voting'
            print("\nTop ranked (by voting score):")

        # Display results
        display_cols = ['path.2.name', 'path.2.type']
        if key_col in result.columns:
            display_cols.append(key_col)
        available = [c for c in display_cols if c in result.columns]
        print(result[available].head(10).to_string())

        # Save output
        if args.output:
            result.to_csv(args.output, index=False)
            print(f"\nSaved to {args.output}")

        # Also print summary
        print(f"\nTotal results: {len(result)}")

    except FileNotFoundError as e:
        print(f"Error: File not found: {e.filename}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _cli_main()
