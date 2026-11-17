import argparse
import ast
import sys

import pandas as pd
import networkx as nx
from scipy import sparse
import scanpy as sc
import numpy as np


def _coerce_gene_list(value):
    """Normalize a scalar or stringified list from CSV into gene names."""
    if isinstance(value, list):
        return [str(x) for x in value]
    if isinstance(value, tuple) or isinstance(value, set):
        return [str(x) for x in value]
    if pd.isna(value):
        return []
    if isinstance(value, str):
        text = value.strip()
        if text.startswith("[") and text.endswith("]"):
            try:
                parsed = ast.literal_eval(text)
                if isinstance(parsed, (list, tuple, set)):
                    return [str(x) for x in parsed]
            except (SyntaxError, ValueError):
                pass
        return [text]
    return [str(value)]

def query_df2nxg(df, directed=False, reverse_direct=False):
    """
    Converts a query result DataFrame (`df`) into a NetworkX graph (`nx.Graph` or `nx.DiGraph`).

    This function transforms a DataFrame containing knowledge graph query results into
    a network structure where nodes represent entities and edges represent relationships
    between entities.

    Args:
        df (pd.DataFrame): A DataFrame containing query path results, typically retrieved
                           from Neo4j or other knowledge graph databases.
                           The DataFrame must include the following columns:
                           - `path.0.name`: Name of the starting node in the path
                           - `path.0.type`: Type of the starting node
                           - `path.2.name`: Name of the ending node in the path
                           - `path.2.type`: Type of the ending node
        directed (bool, optional): Whether to create a directed graph. Defaults to `False` (undirected graph).
        reverse_direct (bool, optional): If `directed=True`, this parameter determines whether to
                                         reverse the edge direction (making the endpoint point to the start point).
                                         Defaults to `False`.

    Returns:
        networkx.Graph or networkx.DiGraph: A NetworkX graph representing the query result.
    """
    if directed:
        nxg = nx.DiGraph()
    else:
        nxg = nx.Graph()

    for x, y in zip(df['path.0.name'], df['path.0.type']):
        if y == 'Gene':
            nxg.add_node(x, type=y)
        else:
            nxg.add_node(x, type="KG_" + y)

    for x, y in zip(df['path.2.name'], df['path.2.type']):
        if y == 'Gene':
            nxg.add_node(x, type=y)
        else:
            nxg.add_node(x, type="KG_" + y)

    for x, y in zip(df['path.0.name'], df['path.2.name']):
        if reverse_direct:
            nxg.add_edge(y, x)
        else:
            nxg.add_edge(x, y)
    return nxg


def count_df2nxg(df, directed=False, reverse_direct=False):
    """
    Converts a query result DataFrame (`df`) into a NetworkX graph (`nx.Graph` or `nx.DiGraph`).

    This function processes associations between genes and entities in a knowledge graph,
    constructing a network representation:
    - `path.0.name` may contain multiple genes (in list format).
    - `path.2.name` represents target entities in the knowledge graph (e.g., diseases, drugs).
    - The resulting graph structure helps explore relationships between genes and
      knowledge graph entities.

    Args:
        df (pd.DataFrame): A DataFrame containing query path results, typically retrieved
                           from a knowledge graph database. The DataFrame must include:
                           - `path.0.name`: A list of genes, where each cell contains multiple genes.
                           - `path.2.name`: Entity names in the knowledge graph (e.g., diseases, drugs).
                           - `path.2.type`: The type of knowledge graph entity (e.g., `'Disease'`, `'Drug'`).
        directed (bool, optional): Whether to create a directed graph. Defaults to `False` (undirected graph).
        reverse_direct (bool, optional): If `directed=True`, this parameter determines whether to
                                         reverse the edge direction (making the endpoint point to the start point).
                                         Defaults to `False`.

    Returns:
        networkx.Graph or networkx.DiGraph: A NetworkX graph representing the query result.
    """
    if directed:
        nxg = nx.DiGraph()
    else:
        nxg = nx.Graph()

    df = df.copy()
    df['path.0.name'] = df['path.0.name'].apply(_coerce_gene_list)

    tmp = set()
    for x in df['path.0.name']:
        for y in x:
            tmp.add(y)
    for x in tmp:
        nxg.add_node(x, type='Gene')
    for x, y in zip(df['path.2.name'], df['path.2.type']):
        nxg.add_node(x, type="KG_"+y)

    for x, y in zip(df['path.0.name'], df['path.2.name']):
        for z in x:
            if reverse_direct:
                nxg.add_edge(y, z)
            else:
                nxg.add_edge(z, y)
    return nxg


def df2nxg(df, directed=False, reverse_direct=False):
    """
    Converts a query result DataFrame (`df`) into a NetworkX graph (`nx.Graph` or `nx.DiGraph`).

    This function determines whether the input DataFrame represents gene-entity associations
    or general entity relationships and selects the appropriate conversion method:
    - If `path.0.name` contains a list (indicating multiple genes per entry),
      it calls `count_df2nxg()` to construct a gene-to-entity network.
    - Otherwise, it calls `query_df2nxg()` to create a general knowledge graph.

    Args:
        df (pd.DataFrame): A DataFrame containing query path results from a knowledge graph database.
                           It must contain the column `path.0.name`.
        directed (bool, optional): Whether to create a directed graph. Defaults to `False` (undirected graph).
        reverse_direct (bool, optional): If `directed=True`, this parameter determines whether to
                                         reverse the edge direction (making the endpoint point to the start point).
                                         Defaults to `False`.

    Returns:
        networkx.Graph or networkx.DiGraph: A NetworkX graph representing the query result.
    """
    tmp = list(df['path.0.name'])[0]
    if isinstance(tmp, list) or (isinstance(tmp, str) and tmp.strip().startswith("[") and tmp.strip().endswith("]")):
        return count_df2nxg(df, directed=directed, reverse_direct=reverse_direct)
    else:
        return query_df2nxg(df, directed=directed, reverse_direct=reverse_direct)


def adata2nxg(adata, directed=False, reverse_direct=False, cell_neighbor=False):
    """
    Converts an AnnData (`adata`) object into a NetworkX graph (`nx.Graph` or `nx.DiGraph`).

    This function constructs a network representation of single-cell transcriptomics data, including:
    - Cell-gene expression relationships (extracted from `adata.X`).
    - Optional cell-cell adjacency relationships (extracted from `adata.obsp['connectivities']` if `cell_neighbor=True`).

    Args:
        adata (AnnData): A single-cell data object that must contain:
            - `adata.X`: Cell-gene expression matrix (in sparse matrix format).
            - `adata.obs_names`: Cell names.
            - `adata.var_names`: Gene names.
            - `adata.obsp['connectivities']` (optional): Cell adjacency matrix (used when `cell_neighbor=True`).
        directed (bool, optional): Whether to create a directed graph. Defaults to `False` (undirected graph).
        reverse_direct (bool, optional): If `directed=True`, this parameter determines whether to reverse
                                         the edge direction (making the endpoint point to the start point).
                                         Defaults to `False`.
        cell_neighbor (bool, optional): Whether to add cell-cell adjacency relationships based on
                                        `adata.obsp['connectivities']`. Defaults to `False`.

    Returns:
        networkx.Graph or networkx.DiGraph: A NetworkX graph representing the single-cell data.
    """
    if sparse.issparse(adata.X):
        expr = adata.X.tocoo()
    else:
        expr = sparse.coo_array(adata.X)
        print('adata.X seems not a sparse matrix, please check whether the adata.X is raw count or log normalized count')
    Omics_edges = []
    Omics_nodes = {}

    celllist = adata.obs_names
    genelist = adata.var_names

    for x in genelist:
        Omics_nodes[x] = {'type': 'Gene'}

    for x in celllist:
        Omics_nodes[x] = {'type': 'Omics_Cell'}

    for x, y in zip(expr.row, expr.col):
        Omics_edges.append([celllist[x], genelist[y]])

    if cell_neighbor:
        if 'connectivities' in adata.obsp:
            neighbor = adata.obsp['connectivities'].tocoo()
            for x, y in zip(neighbor.row, neighbor.col):
                Omics_edges.append([celllist[x], celllist[y]])
        else:
            print('the adata must run sc.pp.neighbor first when cell_neighbor=True, skipping cell neighbor edges')

    Omics_nodes = [[x, y] for x, y in Omics_nodes.items()]
    if directed:
        nxg = nx.DiGraph()
    else:
        nxg = nx.Graph()
    #nxg = nx.Graph()
    nxg.add_nodes_from(Omics_nodes)
    if reverse_direct:
        nxg.add_edges_from(Omics_edges)
    else:
        nxg.add_edges_from([x[::-1] for x in Omics_edges])
    return nxg


def combine_graph(adata, df, directed=False, omics_reverse_direct=False, kg_reverse_direct=False, cell_neighbor=False):
    """
    Constructs an integrated network graph by combining single-cell omics data (`adata`) and a knowledge graph (`df`).

    This function merges single-cell transcriptomics data with knowledge graph data to create a unified `networkx` graph
    for further analysis and visualization.

    Args:
        adata (AnnData): An `AnnData` object containing single-cell transcriptomics data.
        df (pandas.DataFrame): A DataFrame containing edge data from a knowledge graph.
        directed (bool, optional): Whether to create a directed graph. Defaults to `False`.
        omics_reverse_direct (bool, optional): Whether to reverse the edge direction in the omics network. Defaults to `False`.
        kg_reverse_direct (bool, optional): Whether to reverse the edge direction in the knowledge graph network. Defaults to `False`.
        cell_neighbor (bool, optional): Whether to include cell-cell adjacency information (adds edges between neighboring cells). Defaults to `False`.

    Returns:
        networkx.Graph or networkx.DiGraph: An integrated network graph combining single-cell data and the knowledge graph.
    """
    Omics_Graph = adata2nxg(adata, directed=directed, reverse_direct=omics_reverse_direct, cell_neighbor=cell_neighbor)
    KG_Graph = df2nxg(df, directed=directed, reverse_direct=kg_reverse_direct)
    return nx.compose(Omics_Graph, KG_Graph)


def _prepare_integrated_graph(adata_sub_path, query_df_path, cell_neighbor=False):
    """Load inputs, retain query genes present in AnnData, and compose the graph."""
    adata = sc.read(adata_sub_path)
    df = pd.read_csv(query_df_path)
    required_columns = {"path.0.name", "path.0.type", "path.2.name", "path.2.type"}
    missing_columns = sorted(required_columns - set(df.columns))
    if missing_columns:
        raise ValueError(
            "The query CSV is missing required columns: " + ", ".join(missing_columns)
        )
    if df.empty:
        raise ValueError("The query CSV contains no rows.")

    # gene subset progress
    h5ad_gene_list = list(adata.var_names)
    h5ad_gene_set = set(h5ad_gene_list)

    delete_row = set()
    if "path.0.type" in df:
        source_gene_missing = (
            (df["path.0.type"] == "Gene")
            & (~df["path.0.name"].apply(lambda value: any(g in h5ad_gene_set for g in _coerce_gene_list(value))))
        )
        for x in df.loc[source_gene_missing].index:
            delete_row.add(x)

    if "path.2.type" in df:
        target_gene_missing = (
            (df["path.2.type"] == "Gene")
            & (~df["path.2.name"].apply(lambda value: any(g in h5ad_gene_set for g in _coerce_gene_list(value))))
        )
        for x in df.loc[target_gene_missing].index:
            delete_row.add(x)

    df_sub = df.loc[[x for x in df.index if x not in delete_row]]
    if df_sub.empty:
        raise ValueError("No query rows contain genes present in the AnnData var_names.")

    validate_gene_set = set()

    for x in df_sub["path.0.name"]:
        for y in _coerce_gene_list(x):
            validate_gene_set.add(y)

    for x in df_sub.loc[df_sub["path.2.type"] == "Gene", "path.2.name"]:
        for y in _coerce_gene_list(x):
            validate_gene_set.add(y)

    # Preserve AnnData variable order for deterministic graph construction.
    validate_gene_list = [x for x in h5ad_gene_list if x in validate_gene_set]
    if not validate_gene_list:
        raise ValueError("No query genes were found in the AnnData var_names.")
    adata_sub = adata[:, validate_gene_list].copy()
    nxg = combine_graph(adata_sub, df_sub, cell_neighbor=cell_neighbor)
    return adata_sub, nxg


def _graph_obs(nxg):
    node_names = list(nxg.nodes)
    obs = pd.DataFrame(
        {
            "entity_type": [str(nxg.nodes[node].get("type", "Unknown")) for node in node_names],
            "entity_name": [str(node) for node in node_names],
        },
        index=pd.Index(node_names, dtype=str),
    )
    return node_names, obs


def _pairwise_distance(left, right, metric="euclidean"):
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    if metric == "euclidean":
        left_squared = np.sum(left * left, axis=1, keepdims=True)
        right_squared = np.sum(right * right, axis=1, keepdims=True).T
        return np.sqrt(
            np.maximum(left_squared + right_squared - 2.0 * (left @ right.T), 0.0)
        )
    if metric == "cosine":
        left_norm = np.linalg.norm(left, axis=1, keepdims=True)
        right_norm = np.linalg.norm(right, axis=1, keepdims=True)
        left_unit = np.divide(left, left_norm, out=np.zeros_like(left), where=left_norm > 0)
        right_unit = np.divide(right, right_norm, out=np.zeros_like(right), where=right_norm > 0)
        return np.clip(1.0 - left_unit @ right_unit.T, 0.0, 2.0)
    raise ValueError("distance metric must be 'euclidean' or 'cosine'.")


def _add_kg_omics_distance(adata_graph, embedding_key, metric):
    entity_types = adata_graph.obs["entity_type"].astype(str)
    kg_mask = entity_types.str.startswith("KG_")
    omics_mask = entity_types.str.startswith("Omics_")
    if not kg_mask.any():
        raise ValueError("The integrated graph contains no KG entity nodes.")
    if not omics_mask.any():
        raise ValueError("The integrated graph contains no omics cell nodes.")

    representation = np.asarray(adata_graph.obsm[embedding_key], dtype=np.float64)
    kg_representation = representation[np.asarray(kg_mask)]
    omics_representation = representation[np.asarray(omics_mask)]
    distances = _pairwise_distance(kg_representation, omics_representation, metric=metric)
    adata_graph.uns["KG_Omics_distance"] = pd.DataFrame(
        distances.T,
        index=adata_graph.obs_names[omics_mask],
        columns=adata_graph.obs_names[kg_mask],
    )


def _copy_omics_obs(adata_graph, adata_sub):
    for column in adata_sub.obs.columns:
        values = adata_sub.obs[column].reindex(adata_graph.obs_names)
        adata_graph.obs[column] = values.to_numpy()


def UMAP_graph(adata_sub_path, query_df_path, cell_neighbor=False, distance_metric="euclidean"):
    """Build the legacy adjacency/PCA/UMAP integrated graph embedding."""
    adata_sub, nxg = _prepare_integrated_graph(
        adata_sub_path,
        query_df_path,
        cell_neighbor=cell_neighbor,
    )

    # Build adjacency matrix from the merged graph.
    node_names, obs = _graph_obs(nxg)
    adj_matrix = sparse.csr_matrix(nx.adjacency_matrix(nxg, nodelist=node_names))

    adata_graph = sc.AnnData(adj_matrix, obs=obs, var=obs)
    sc.pp.pca(adata_graph)
    adata_graph.obsp['connectivities'] = adj_matrix
    adata_graph.obsp["distances"] = adj_matrix
    adata_graph.uns["neighbors"] = {'connectivities_key': 'connectivities',
                                    'distances_key': 'distances',
                                    'params': {'n_neighbors': int(len(adj_matrix.data) / adj_matrix.shape[0]),
                                    'method': 'umap',
                                    'random_state': 0,
                                    'metric': 'euclidean'}}
    sc.tl.umap(adata_graph)
    _add_kg_omics_distance(adata_graph, "X_pca", metric=distance_metric)
    _copy_omics_obs(adata_graph, adata_sub)
    adata_graph.uns["brick_embedding"] = {
        "method": "umap",
        "embedding_key": "X_pca",
        "distance_metric": distance_metric,
        "node_count": int(nxg.number_of_nodes()),
        "edge_count": int(nxg.number_of_edges()),
        "cell_neighbor": bool(cell_neighbor),
    }
    return adata_graph


def GRACE_graph(
    adata_sub_path,
    query_df_path,
    *,
    cell_neighbor=False,
    distance_metric="cosine",
    hidden_dim=32,
    num_layers=2,
    learning_rate=0.0001,
    epochs=300,
    edge_drop_rate=0.3,
    feature_mask_rate=0.3,
    temperature=0.2,
    feature_mode="projected_identity",
    feature_dim=64,
    identity_max_nodes=4096,
    max_nodes=5000,
    device="auto",
    seed=20241204,
    verbose=True,
):
    """Build an integrated graph and learn a GRACE contrastive embedding."""
    from grace import train_grace

    adata_sub, nxg = _prepare_integrated_graph(
        adata_sub_path,
        query_df_path,
        cell_neighbor=cell_neighbor,
    )
    result = train_grace(
        nxg,
        hidden_dim=hidden_dim,
        num_layers=num_layers,
        learning_rate=learning_rate,
        epochs=epochs,
        edge_drop_rate=edge_drop_rate,
        feature_mask_rate=feature_mask_rate,
        temperature=temperature,
        feature_mode=feature_mode,
        feature_dim=feature_dim,
        identity_max_nodes=identity_max_nodes,
        max_nodes=max_nodes,
        device=device,
        seed=seed,
        verbose=verbose,
    )

    node_names, obs = _graph_obs(nxg)
    adata_graph = sc.AnnData(X=result.embedding.copy(), obs=obs)
    adata_graph.obsm["X_grace"] = result.embedding
    adjacency = sparse.csr_matrix(nx.adjacency_matrix(nxg, nodelist=node_names))
    adata_graph.obsp["integrated_adjacency"] = adjacency

    if adata_graph.n_obs < 3:
        raise ValueError("GRACE UMAP requires an integrated graph with at least 3 nodes.")
    neighbor_count = min(15, adata_graph.n_obs - 1)
    sc.pp.neighbors(
        adata_graph,
        n_neighbors=neighbor_count,
        use_rep="X_grace",
        random_state=seed,
    )
    sc.tl.umap(adata_graph, random_state=seed)
    _add_kg_omics_distance(adata_graph, "X_grace", metric=distance_metric)
    _copy_omics_obs(adata_graph, adata_sub)

    grace_metadata = dict(result.metadata)
    grace_metadata["losses"] = np.asarray(result.losses, dtype=np.float32)
    grace_metadata["distance_metric"] = distance_metric
    grace_metadata["cell_neighbor"] = bool(cell_neighbor)
    adata_graph.uns["grace"] = grace_metadata
    adata_graph.uns["brick_embedding"] = {
        "method": "grace",
        "embedding_key": "X_grace",
        "distance_metric": distance_metric,
        "node_count": int(nxg.number_of_nodes()),
        "edge_count": int(nxg.number_of_edges()),
        "cell_neighbor": bool(cell_neighbor),
        "seed": int(seed),
    }
    return adata_graph


def embed_graph(adata_path, query_path, method="umap", **kwargs):
    """Dispatch to the requested integrated graph embedding method."""
    if method == "umap":
        allowed = {"cell_neighbor", "distance_metric"}
        return UMAP_graph(
            adata_path,
            query_path,
            **{key: value for key, value in kwargs.items() if key in allowed},
        )
    if method == "grace":
        return GRACE_graph(adata_path, query_path, **kwargs)
    raise ValueError("method must be 'umap' or 'grace'.")

def _cli_main():
    """CLI entry point for brick emb command."""
    parser = argparse.ArgumentParser(
        description=(
            "brick emb - Build an integrated omics/KG graph and embed it with "
            "the legacy PCA/UMAP path or GRACE graph contrastive learning"
        ),
        usage=(
            "brick emb -a <adata.h5ad> -i <query.csv/rank.csv> "
            "[-m umap|grace] -o <adata_graph.h5ad>"
        ),
    )

    parser.add_argument(
        "-a",
        "--adata",
        required=True,
        help="Input AnnData file (.h5ad)",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input CSV file from brick qr or brick rk (query.csv or rank.csv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output integrated AnnData file",
    )
    parser.add_argument(
        "-m",
        "--method",
        choices=("umap", "grace"),
        default="umap",
        help="Embedding method (default: umap; use grace for graph contrastive learning)",
    )
    parser.add_argument(
        "--cell-neighbor",
        action="store_true",
        help="Include adata.obsp['connectivities'] cell-cell edges when available",
    )
    parser.add_argument(
        "--distance-metric",
        choices=("euclidean", "cosine"),
        default=None,
        help="KG-to-omics distance (default: euclidean for umap, cosine for grace)",
    )
    grace = parser.add_argument_group("GRACE options")
    grace.add_argument("--grace-epochs", type=int, default=300)
    grace.add_argument("--grace-hidden-dim", type=int, default=32)
    grace.add_argument("--grace-layers", type=int, default=2)
    grace.add_argument("--grace-learning-rate", type=float, default=0.0001)
    grace.add_argument("--grace-edge-drop", type=float, default=0.3)
    grace.add_argument("--grace-feature-mask", type=float, default=0.3)
    grace.add_argument("--grace-temperature", type=float, default=0.2)
    grace.add_argument(
        "--grace-feature-mode",
        choices=("projected_identity", "identity"),
        default="projected_identity",
        help="Scalable projected identity (default) or upstream dense one-hot identity",
    )
    grace.add_argument("--grace-feature-dim", type=int, default=64)
    grace.add_argument("--grace-identity-max-nodes", type=int, default=4096)
    grace.add_argument(
        "--grace-max-nodes",
        type=int,
        default=5000,
        help="Safety cap for the all-node contrastive objective; 0 disables the cap",
    )
    grace.add_argument(
        "--grace-device",
        choices=("auto", "cpu", "cuda", "mps"),
        default="auto",
    )
    grace.add_argument("--seed", type=int, default=20241204)
    grace.add_argument("--quiet", action="store_true", help="Hide GRACE epoch progress")

    args = parser.parse_args()

    try:
        print(f"Loading AnnData: {args.adata}")
        print(f"Loading query CSV: {args.input}")
        distance_metric = args.distance_metric or (
            "cosine" if args.method == "grace" else "euclidean"
        )
        print(f"Embedding method: {args.method}")
        adata_graph = embed_graph(
            args.adata,
            args.input,
            method=args.method,
            cell_neighbor=args.cell_neighbor,
            distance_metric=distance_metric,
            hidden_dim=args.grace_hidden_dim,
            num_layers=args.grace_layers,
            learning_rate=args.grace_learning_rate,
            epochs=args.grace_epochs,
            edge_drop_rate=args.grace_edge_drop,
            feature_mask_rate=args.grace_feature_mask,
            temperature=args.grace_temperature,
            feature_mode=args.grace_feature_mode,
            feature_dim=args.grace_feature_dim,
            identity_max_nodes=args.grace_identity_max_nodes,
            max_nodes=args.grace_max_nodes,
            device=args.grace_device,
            seed=args.seed,
            verbose=not args.quiet,
        )

        adata_graph.write(args.output)
        print(f"Saved graph AnnData to {args.output}")

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
