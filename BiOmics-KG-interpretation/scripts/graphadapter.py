"""Agent_brick graph adapter module - Unified graph conversion and I/O.

Consolidates brick preprocess.py graph functions into 4 unified APIs:
- read_graph:  Read graph files (.svg, .pkl, .graphml)
- save_graph:  Save graphs to files (.pkl, .svg, .pdf, .graphml)
- data2nxg:    Convert data (DataFrame/AnnData/PAGA) to NetworkX graph
- combine_graph: Combine omics and knowledge graphs
"""

import os
import argparse
import sys
import pickle
import ast
import xml.etree.ElementTree as ET

import networkx as nx
import pandas as pd
from networkx.drawing.nx_pydot import to_pydot
from scipy import sparse


# ─── Internal helpers ────────────────────────────────────────────────

def _remove_namespace(tree):
    """Remove namespaces from XML tags for easier SVG parsing."""
    for elem in tree.iter():
        if '}' in elem.tag:
            elem.tag = elem.tag.split('}', 1)[1]
    return tree


def _coerce_list(value):
    """Normalize a scalar or stringified list from CSV into a Python list."""
    if isinstance(value, list):
        return value
    if isinstance(value, tuple) or isinstance(value, set):
        return list(value)
    if pd.isna(value):
        return []
    if isinstance(value, str):
        text = value.strip()
        if text.startswith('[') and text.endswith(']'):
            try:
                parsed = ast.literal_eval(text)
                if isinstance(parsed, (list, tuple, set)):
                    return list(parsed)
            except (SyntaxError, ValueError):
                pass
        return [text]
    return [value]


def _is_list_like_cell(value):
    """Return True for actual lists or CSV cells that look like list literals."""
    if isinstance(value, (list, tuple, set)):
        return True
    return isinstance(value, str) and value.strip().startswith('[') and value.strip().endswith(']')


def _query_df2nxg(df, directed=False, reverse_direct=False):
    """Convert query result DataFrame (entity-entity pairs) to NetworkX graph."""
    nxg = nx.DiGraph() if directed else nx.Graph()

    for x, y in zip(df['path.0.name'], df['path.0.type']):
        nxg.add_node(x, type=y if y == 'Gene' else "KG_" + y)

    for x, y in zip(df['path.2.name'], df['path.2.type']):
        nxg.add_node(x, type=y if y == 'Gene' else "KG_" + y)

    for x, y in zip(df['path.0.name'], df['path.2.name']):
        if reverse_direct:
            nxg.add_edge(y, x)
        else:
            nxg.add_edge(x, y)
    return nxg


def _count_df2nxg(df, directed=False, reverse_direct=False):
    """Convert count DataFrame (gene-list ->entity) to NetworkX graph."""
    nxg = nx.DiGraph() if directed else nx.Graph()
    df = df.copy()
    df['path.0.name'] = df['path.0.name'].apply(_coerce_list)

    tmp = set()
    for x in df['path.0.name']:
        for y in x:
            tmp.add(y)
    for x in tmp:
        nxg.add_node(x, type='Gene')

    for x, y in zip(df['path.2.name'], df['path.2.type']):
        nxg.add_node(x, type="KG_" + y)

    for x, y in zip(df['path.0.name'], df['path.2.name']):
        for z in x:
            if reverse_direct:
                nxg.add_edge(y, z)
            else:
                nxg.add_edge(z, y)
    return nxg


def _read_svg_graph(svg_file, node_class='node', edge_class='edge',
                    default_node_color='black', default_edge_color='black'):
    """Parse SVG file into NetworkX MultiDiGraph."""
    tree = ET.parse(svg_file)
    tree = _remove_namespace(tree)
    root = tree.getroot()

    G = nx.MultiDiGraph()

    for node in root.findall(f".//g[@class='{node_class}']"):
        title = node.find('title')
        title_text = title.text.strip() if title is not None else None
        shape, color = None, default_node_color
        for shape_tag in ['ellipse', 'polygon', 'circle', 'rect']:
            element = node.find(shape_tag)
            if element is not None:
                shape = shape_tag
                color = element.attrib.get('stroke', default_node_color)
                break
        if title_text:
            G.add_node(title_text, shape=shape, color=color)

    for edge in root.findall(f".//g[@class='{edge_class}']"):
        title = edge.find('title')
        title_text = title.text.strip() if title is not None else None
        path = edge.find('path')
        color = default_edge_color
        if path is not None:
            color = path.attrib.get('stroke', default_edge_color)
        if title_text and '->' in title_text:
            source, target = map(str.strip, title_text.split('->'))
            G.add_edge(source, target, color=color)

    return G


# ─── Public API ──────────────────────────────────────────────────────

def read_graph(file_path, **kwargs):
    """
    Read a graph file and return a NetworkX graph object.

    Supports .svg, .pkl, .gpickle, .graphml formats.
    For .svg files, optional kwargs: node_class, edge_class,
    default_node_color, default_edge_color.

    Args:
        file_path: Path to the graph file.
        **kwargs: Additional parameters (for SVG: node_class, edge_class,
                  default_node_color, default_edge_color).

    Returns:
        networkx.Graph or networkx.DiGraph or networkx.MultiDiGraph:
            The parsed graph object.

    Raises:
        ValueError: If file format is unsupported.
    """
    ext = os.path.splitext(file_path)[-1].lower()
    if ext == '.svg':
        return _read_svg_graph(file_path, **kwargs)
    elif ext in ('.pkl', '.gpickle'):
        with open(file_path, 'rb') as f:
            return pickle.load(f)
    elif ext == '.graphml':
        return nx.read_graphml(file_path)
    else:
        raise ValueError(f"Unsupported file format: {ext}. Use .svg, .pkl, .gpickle, or .graphml")


def save_graph(graph, file_name):
    """
    Save a NetworkX graph to file.

    Supports .pkl, .svg, .pdf, .graphml formats.
    Automatically handles edge styling for visual formats (svg/pdf).

    Args:
        graph: NetworkX graph object.
        file_name: Output file name with extension (e.g., 'network.pkl').

    Raises:
        ValueError: If file format is unsupported.
    """
    file_format = file_name.split(".")[-1]

    # Apply edge styling for visual formats
    for u, v, data in graph.edges(data=True):
        if "color" in data:
            if data["color"] == "green":
                data["arrowhead"] = "normal"
            else:
                data["arrowhead"] = "none"

    if file_format == "pkl":
        with open(file_name, "wb") as f:
            pickle.dump(graph, f)
    elif file_format == "svg":
        to_pydot(graph).write_svg(file_name)
    elif file_format == "pdf":
        to_pydot(graph).write_pdf(file_name)
    elif file_format == "graphml":
        nx.write_graphml(graph, file_name)
    else:
        raise ValueError(f"Unsupported format: {file_format}. Use 'pkl', 'svg', 'pdf', or 'graphml'.")

    print(f"Graph successfully saved to {file_name}")


def data2nxg(data, data_type='df', directed=False, reverse_direct=False,
             cell_neighbor=False, **kwargs):
    """
    Convert data to a NetworkX graph.

    Unified interface covering DataFrame, AnnData, and PAGA conversions.

    Args:
        data: Input data object.
        data_type: Type of input data:
            - 'df': DataFrame from query results (covers df2nxg).
                    Auto-detects entity-entity vs gene-list format.
            - 'adata': AnnData single-cell object (covers adata2nxg).
            - 'paga': AnnData with PAGA results (covers paga2nx).
        directed: Whether to create a directed graph (for 'df' and 'adata').
        reverse_direct: Whether to reverse edge direction (for 'df' and 'adata').
        cell_neighbor: Whether to add cell-cell adjacency edges
                       from adata.obsp['connectivities'] (for 'adata' only).
        **kwargs: Reserved for future extensions.

    Returns:
        networkx.Graph or networkx.DiGraph or networkx.MultiDiGraph:
            The converted graph.

    Raises:
        ValueError: If data_type is unsupported.
    """
    if data_type == 'df':
        tmp = list(data['path.0.name'])[0]
        if _is_list_like_cell(tmp):
            return _count_df2nxg(data, directed=directed, reverse_direct=reverse_direct)
        else:
            return _query_df2nxg(data, directed=directed, reverse_direct=reverse_direct)

    elif data_type == 'adata':
        nxg = nx.DiGraph() if directed else nx.Graph()

        if sparse.issparse(data.X):
            expr = data.X.tocoo()
        else:
            expr = sparse.coo_array(data.X)
            print('adata.X seems not a sparse matrix, please check whether it is raw count or log normalized count')

        celllist = data.obs_names
        genelist = data.var_names

        for x in genelist:
            nxg.add_node(x, type='Gene')
        for x in celllist:
            nxg.add_node(x, type='Omics_Cell')

        for x, y in zip(expr.row, expr.col):
            if reverse_direct:
                nxg.add_edge(genelist[y], celllist[x])
            else:
                nxg.add_edge(celllist[x], genelist[y])

        if cell_neighbor:
            if 'connectivities' in data.obsp:
                neighbor = data.obsp['connectivities'].tocoo()
                for x, y in zip(neighbor.row, neighbor.col):
                    nxg.add_edge(celllist[x], celllist[y])
            else:
                print('adata must run sc.pp.neighbors first when cell_neighbor=True, skipping cell neighbor edges')

        return nxg

    elif data_type == 'paga':
        baseline_size_max = 10
        baseline_size = 20
        baseline_weight_max = 5
        baseline_weight = 10

        adj = pd.DataFrame(data.uns['paga']['connectivities_tree'].todense())
        group = data.uns['paga']['groups']
        categories = data.obs[group].cat.categories

        adj.columns = categories
        adj.index = categories
        adj = adj.loc[categories, categories]

        G = nx.from_pandas_adjacency(adj, create_using=nx.MultiDiGraph)

        name2color = {x: y for x, y in zip(categories, data.uns[f'{group}_colors'])}

        tmp = data.obs.copy()
        tmp['0'] = 0
        name2size = tmp[[group, '0']].groupby(group).count().to_dict()['0']
        max_size = max(name2size.values())
        name2size = {x: y * baseline_size_max / max_size + baseline_size
                     for x, y in name2size.items()}

        for x in G.nodes:
            G.nodes[x]['original_name'] = x
            G.nodes[x]['color'] = name2color[x]
            G.nodes[x]['size'] = name2size[x]

        for edge in G.edges:
            e = G.edges[edge]
            e['weight'] = e['weight'] * baseline_weight_max + baseline_weight
            e['color'] = 'black'
            e['arrows'] = False

        return G

    else:
        raise ValueError(f"Unsupported data_type: {data_type}. Use 'df', 'adata', or 'paga'")


def combine_graph(adata, df, directed=False, omics_reverse_direct=False,
                  kg_reverse_direct=False, cell_neighbor=False):
    """
    Construct an integrated network graph by combining single-cell omics
    data (adata) and a knowledge graph (df).

    Args:
        adata: AnnData object containing single-cell transcriptomics data.
        df: DataFrame containing edge data from a knowledge graph.
        directed: Whether to create a directed graph.
        omics_reverse_direct: Whether to reverse edge direction in the
                              omics network.
        kg_reverse_direct: Whether to reverse edge direction in the
                           knowledge graph network.
        cell_neighbor: Whether to include cell-cell adjacency information.

    Returns:
        networkx.Graph or networkx.DiGraph: Combined graph merging
            single-cell data and knowledge graph.
    """
    Omics_Graph = data2nxg(adata, data_type='adata', directed=directed,
                           reverse_direct=omics_reverse_direct,
                           cell_neighbor=cell_neighbor)
    KG_Graph = data2nxg(df, data_type='df', directed=directed,
                        reverse_direct=kg_reverse_direct)
    return nx.compose(Omics_Graph, KG_Graph)


# ─── Internal helpers for complete_result ────────────────────────────

def _filter_graph_paths(qg, node_type_filter=None, reverse_relations=None):
    """Extract normalized paths from a query-result NetworkX graph.

    Returns deduplicated tuples ``(source_name, relation, target_name, metadata)``.
    Metadata keeps edge provenance fields such as info_source, original_relation,
    condition, and relation_confidence so downstream graph analysis does not
    silently drop evidence strength/review context.
    """
    reverse_relations = set(reverse_relations or [])

    if node_type_filter:
        id_to_name = {
            n: qg.nodes[n].get('name', n)
            for n in qg.nodes
            if 'name' in qg.nodes[n]
            and qg.nodes[n].get('type') == node_type_filter
        }
    else:
        id_to_name = {
            n: qg.nodes[n].get('name', n)
            for n in qg.nodes
            if 'name' in qg.nodes[n]
        }

    paths = []
    seen = set()
    for source_id, target_id, edge_data in qg.edges(data=True):
        if source_id not in id_to_name or target_id not in id_to_name:
            continue

        source_name = id_to_name[source_id]
        target_name = id_to_name[target_id]
        relation = edge_data.get('relation', 'unknown')
        metadata = {
            key: value
            for key, value in edge_data.items()
            if key in {
                'info_source', 'info_source_count', 'info_source_length', 'relation_confidence',
                'confidence', 'original_relation', 'condition', 'evidence_text', 'edge_support_class',
                'evidence_strength', 'source_quality', 'manual_review_required', 'path_id', 'edge_position',
            }
        }

        if relation in reverse_relations:
            item = (target_name, relation, source_name, metadata)
        else:
            item = (source_name, relation, target_name, metadata)

        key = (item[0], item[1], item[2], tuple(sorted((str(k), str(v)) for k, v in metadata.items())))
        if key not in seen:
            seen.add(key)
            paths.append(item)

    return paths

def _path_edge_indices(df):
    return sorted(
        {int(match.group(1)) for col in df.columns for match in [__import__('re').match(r"^path\.(\d+)\.relation$", col)] if match}
    )


def _edge_attrs_from_row(row, edge_idx=None):
    """Preserve provenance metadata on graph edges."""
    if edge_idx is None:
        candidates = {
            'info_source': ['path.1.info_source', 'info_source'],
            'info_source_length': ['path.1.info_source_length', 'info_source_count', 'info_source_length'],
            'relation_confidence': ['path.1.relation_confidence', 'confidence', 'relation_confidence'],
            'original_relation': ['path.1.original_relation', 'original_relation', 'evidence_text'],
            'condition': ['path.1.condition', 'condition'],
            'edge_support_class': ['edge_support_class'],
            'evidence_strength': ['evidence_strength'],
            'review_required': ['review_required'],
        }
    else:
        candidates = {
            'info_source': [f'path.{edge_idx}.info_source', 'info_source'],
            'info_source_length': [f'path.{edge_idx}.info_source_length', 'info_source_count', 'info_source_length'],
            'relation_confidence': [f'path.{edge_idx}.relation_confidence', 'confidence', 'relation_confidence'],
            'original_relation': [f'path.{edge_idx}.original_relation', 'original_relation', 'evidence_text'],
            'condition': [f'path.{edge_idx}.condition', 'condition'],
            'edge_support_class': ['edge_support_class'],
            'evidence_strength': ['evidence_strength'],
            'review_required': ['review_required'],
        }
    attrs = {}
    for out_key, cols in candidates.items():
        for col in cols:
            if col in row.index and pd.notna(row.get(col)):
                attrs[out_key] = row.get(col)
                break
    return attrs


def _df_to_paths(df, relation_col=None):
    """Extract normalized paths from a brick query-result DataFrame.

    One-hop and multi-hop brick qr outputs are converted into edge tuples
    ``(source_name, relation, target_name, metadata)``.  Metadata includes
    original_relation/info_source/condition/confidence whenever present.
    """
    if relation_col is None:
        for candidate in ('path.1.relation', 'path.1.name', 'path.1'):
            if candidate in df.columns:
                relation_col = candidate
                break

    edge_indices = []
    for col in df.columns:
        import re
        match = re.match(r'^path\.(\d+)\.(?:relation|name)$', str(col))
        if match:
            idx = int(match.group(1))
            if idx % 2 == 1 and (f'path.{idx - 1}.name' in df.columns) and (f'path.{idx + 1}.name' in df.columns):
                edge_indices.append(idx)
    edge_indices = sorted(set(edge_indices))

    if not edge_indices:
        required = ['path.0.name', 'path.2.name']
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(
                f"DataFrame missing required columns: {missing}. "
                f"Available columns: {list(df.columns)}"
            )
        edge_indices = [1]

    paths = []
    seen = set()
    for _, row in df.iterrows():
        for edge_idx in edge_indices:
            source_col = f'path.{edge_idx - 1}.name'
            target_col = f'path.{edge_idx + 1}.name'
            rel_col = f'path.{edge_idx}.relation' if f'path.{edge_idx}.relation' in df.columns else relation_col
            if source_col not in df.columns or target_col not in df.columns:
                continue
            source = row.get(source_col)
            target = row.get(target_col)
            relation = row.get(rel_col, 'unknown') if rel_col else 'unknown'
            if pd.isna(source) or pd.isna(target) or pd.isna(relation):
                continue
            if not str(source).strip() or not str(target).strip() or not str(relation).strip():
                continue
            metadata = {}
            for out_key, candidates in {
                'info_source': [f'path.{edge_idx}.info_source', 'info_source'],
                'info_source_count': [f'path.{edge_idx}.info_source_length', 'info_source_count', 'info_source_length'],
                'relation_confidence': [f'path.{edge_idx}.relation_confidence', 'relation_confidence', 'confidence'],
                'condition': [f'path.{edge_idx}.condition', 'condition'],
                'original_relation': [f'path.{edge_idx}.original_relation', 'original_relation', 'evidence_text'],
                'evidence_strength': ['evidence_strength'],
                'edge_support_class': ['edge_support_class'],
                'source_quality': ['source_quality'],
                'manual_review_required': ['manual_review_required'],
                'path_id': ['path_id'],
                'edge_position': ['edge_position'],
            }.items():
                for col in candidates:
                    if col in df.columns and pd.notna(row.get(col)):
                        metadata[out_key] = row.get(col)
                        break
            # Restore the true KG edge direction when the query API provides
            # it (path.N.reversed): path column order is traversal order and
            # may oppose the actual relationship direction (e.g. is_a points
            # child -> parent regardless of traversal direction).
            rev_col = f'path.{edge_idx}.reversed'
            if rev_col in df.columns:
                rev_val = row.get(rev_col)
                if pd.notna(rev_val) and str(rev_val).strip().lower() in ('true', '1', 'yes'):
                    source, target = target, source
            item = (source, relation, target, metadata)
            key = (str(source), str(relation), str(target), tuple(sorted((str(k), str(v)) for k, v in metadata.items())))
            if key not in seen:
                seen.add(key)
                paths.append(item)
    return paths

def complete_result(original_graph, query_results=None,
                    node_attrs=None, edge_attrs=None,
                    node_type_filter=None, reverse_relations=None,
                    auto_filter=True, highlight_disconnected=True,
                    disconnected_edge_attrs=None):
    """Merge an original graph with query results into an annotated graph.

    This function replicates the workflow in ``Trajectory.ipynb``
    (``G`` ->``nx_new_paga``) but is generalized for any original graph
    and any type of query result.

    ``query_results`` is auto-detected and converted to a list of
    ``(source, relation, target)`` tuples:

    * **NetworkX Graph** (e.g. output of ``query_shortest_path`` with
      ``return_type='graph'``): automatically passed through
      :func:`_filter_graph_paths` when ``auto_filter=True``.
    * **pandas DataFrame** (e.g. output of ``brick qr`` CSV):
      converted by :func:`_df_to_paths`.
    * **list of tuples**: used directly.
    * **None**: only the original graph is processed (useful for
      highlighting disconnected edges).

    The returned graph is a :class:`networkx.MultiDiGraph` that contains:

    1. All original nodes and their attributes.
    2. New nodes discovered through query results (assigned ``node_attrs``).
    3. Query-result edges (assigned ``edge_attrs``).
    4. *(Optional)* Red edges for original edges that remain disconnected
       after merging.

    Args:
        original_graph (nx.Graph): Base graph (e.g. PAGA graph).
        query_results: Query results.  Acceptable types:
            - ``nx.Graph`` / ``nx.MultiDiGraph``
            - ``pd.DataFrame``
            - ``list`` of ``(source, relation, target)`` tuples
            - ``None``
        node_attrs (dict, optional): Default attributes for nodes that are
            introduced by ``query_results`` but absent from
            ``original_graph``.  Example:
            ``{'color': 'gray', 'shape': 'hexagon', 'size': 12}``.
        edge_attrs (dict, optional): Attributes applied to every edge that
            comes from ``query_results``.  Example:
            ``{'color': 'green', 'weight': 8, 'arrows': 'to'}``.
        node_type_filter (str, optional): When ``query_results`` is a graph,
            only keep nodes whose ``type`` attribute matches this value.
            Typical value: ``'Cell'``.
        reverse_relations (list, optional): When ``query_results`` is a graph,
            relations in this list have their direction reversed.
            Typical value: ``['is_a']``.
        auto_filter (bool): If ``True`` and ``query_results`` is a graph,
        highlight_disconnected (bool): If ``True``, original edges that do
            **not** participate in any path within the merged graph are
            added back with a distinct red style (as in
            ``Trajectory.ipynb``).
        disconnected_edge_attrs (dict, optional): Attributes for the
            highlighted disconnected edges.  Defaults to
            ``{'color': 'red', 'weight': 2, 'arrows': 'none'}``.

    Returns:
        nx.MultiDiGraph: Consolidated annotated graph.

    Examples:
        Trajectory use-case (equivalent to the notebook workflow)::

            all_paths = []
            for s, t in itertools.combinations(paga_G.nodes, 2):
                qg = brick.qr.query_shortest_path(
                    source_entity_name=s, target_entity_name=t,
                    source_entity_type='Cell', target_entity_type='Cell',
                    relation='is_a|develope_into',
                    multi_hop=(1, 10), return_type='graph')
                all_paths.extend(qg)   # list of graphs

            nx_new = complete_result(
                original_graph=paga_G,
                query_results=all_paths,          # auto-filtered
                node_attrs={'color': 'gray', 'shape': 'hexagon', 'size': 12},
                edge_attrs={'color': 'green', 'weight': 8, 'arrows': 'to'},
                node_type_filter='Cell',
                reverse_relations=['is_a'])

        CLI-equivalent use-case (CSV from ``brick qr``)::

            df = pd.read_csv('query_result.csv')
            nx_new = complete_result(
                original_graph=paga_G,
                query_results=df,
                node_attrs={'color': 'gray'},
                edge_attrs={'color': 'green', 'arrows': 'to'})
    """
    node_attrs = node_attrs or {}
    edge_attrs = edge_attrs or {}
    disconnected_edge_attrs = disconnected_edge_attrs or {
        'color': 'red', 'weight': 2, 'arrows': 'none'
    }

    # ── 1. Normalise query_results ->list of (source, relation, target) ──
    filtered_paths = []

    if query_results is None:
        pass

    elif isinstance(query_results, list):
        # Distinguish: list of graphs vs list of tuples
        if len(query_results) > 0 and isinstance(query_results[0], nx.Graph):
            for qg in query_results:
                if auto_filter:
                    filtered_paths.extend(
                        _filter_graph_paths(
                            qg,
                            node_type_filter=node_type_filter,
                            reverse_relations=reverse_relations
                        )
                    )
                else:
                    # Use raw edge names when auto_filter=False
                    for s, t, d in qg.edges(data=True):
                        filtered_paths.append(
                            (str(s), d.get('relation', 'unknown'), str(t), dict(d))
                        )
        else:
            # Assume already a list of tuples
            filtered_paths = list(query_results)

    elif isinstance(query_results, nx.Graph):
        if auto_filter:
            filtered_paths = _filter_graph_paths(
                query_results,
                node_type_filter=node_type_filter,
                reverse_relations=reverse_relations
            )
        else:
            for s, t, d in query_results.edges(data=True):
                filtered_paths.append(
                    (str(s), d.get('relation', 'unknown'), str(t), dict(d))
                )

    elif isinstance(query_results, pd.DataFrame):
        filtered_paths = _df_to_paths(query_results)

    else:
        raise TypeError(
            f"Unsupported query_results type: {type(query_results).__name__}. "
            "Expected nx.Graph, pd.DataFrame, list, or None."
        )

    # ── 2. Build merged MultiDiGraph ──
    G = nx.MultiDiGraph()

    # Preserve original nodes and their attributes
    for n, data in original_graph.nodes(data=True):
        G.add_node(n, **data)

    # Add query-result nodes and edges. Keep KG evidence metadata on edges so
    # downstream visualization/reasoning can show provenance and weak-evidence notes.
    for item in filtered_paths:
        if len(item) >= 4:
            s, relation, t, evidence_attrs = item[0], item[1], item[2], dict(item[3] or {})
        else:
            s, relation, t = item[:3]
            evidence_attrs = {}
        if s not in G:
            G.add_node(s, **node_attrs)
        if t not in G:
            G.add_node(t, **node_attrs)
        attrs = dict(edge_attrs)
        attrs.update({k: v for k, v in evidence_attrs.items() if v is not None})
        attrs.setdefault('relation', relation)
        G.add_edge(s, t, **attrs)

    # Highlight disconnected original edges
    if highlight_disconnected:
        for s, t in original_graph.edges():
            if not (nx.has_path(G, s, t) or nx.has_path(G, t, s)):
                G.add_edge(s, t, **disconnected_edge_attrs)

    return G


# ─── CLI ─────────────────────────────────────────────────────────────

def _cli_main():
    """CLI entry point for brick ga command."""
    parser = argparse.ArgumentParser(
        description="brick ga - Graph adapter: read, convert, save, and combine graphs",
        usage='brick ga <command> [options]'
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # --- read command ---
    read_parser = subparsers.add_parser('read', help='Read a graph file and display info')
    read_parser.add_argument("-i", "--input", required=True,
                             help="Input graph file (.svg/.pkl/.graphml)")

    # --- convert command ---
    convert_parser = subparsers.add_parser('convert',
                                           help='Convert graph file to another format')
    convert_parser.add_argument("-i", "--input", required=True,
                                help="Input graph file")
    convert_parser.add_argument("-o", "--output", required=True,
                                help="Output graph file (.pkl/.svg/.pdf/.graphml)")

    # --- data2nxg command ---
    data_parser = subparsers.add_parser('data2nxg',
                                        help='Convert data to NetworkX graph file')
    data_parser.add_argument("-i", "--input", required=True,
                             help="Input file (CSV for df, .h5ad for adata/paga)")
    data_parser.add_argument("-o", "--output", required=True,
                             help="Output graph file (.pkl/.svg/.graphml)")
    data_parser.add_argument("-t", "--type", required=True,
                             choices=['df', 'adata', 'paga'],
                             help="Input data type: df (DataFrame), adata (AnnData), or paga (PAGA)")
    data_parser.add_argument("--directed", action="store_true",
                             help="Create directed graph")
    data_parser.add_argument("--reverse", action="store_true",
                             help="Reverse edge direction")
    data_parser.add_argument("--cell-neighbor", action="store_true",
                             help="Add cell-cell adjacency edges (adata only)")

    # --- combine command ---
    combine_parser = subparsers.add_parser('combine',
                                           help='Combine omics and knowledge graphs')
    combine_parser.add_argument("-a", "--adata", required=True,
                                help="Input AnnData file (.h5ad)")
    combine_parser.add_argument("-d", "--df", required=True,
                                help="Input DataFrame CSV file (query result)")
    combine_parser.add_argument("-o", "--output", required=True,
                                help="Output graph file (.pkl/.svg/.graphml)")
    combine_parser.add_argument("--directed", action="store_true",
                                help="Create directed graph")
    combine_parser.add_argument("--omics-reverse", action="store_true",
                                help="Reverse edge direction in omics network")
    combine_parser.add_argument("--kg-reverse", action="store_true",
                                help="Reverse edge direction in knowledge graph")
    combine_parser.add_argument("--cell-neighbor", action="store_true",
                                help="Add cell-cell adjacency edges")

    # --- complete_result command ---
    cr_parser = subparsers.add_parser(
        'complete_result',
        help='Merge original graph with query results into an annotated graph'
    )
    cr_parser.add_argument("-g", "--graph", required=True,
                           help="Original graph file (.graphml/.pkl)")
    cr_parser.add_argument("-q", "--query", default=None,
                           help="Query result file (.csv/.graphml/.pkl). "
                                "If omitted, only processes the original graph.")
    cr_parser.add_argument("-o", "--output", required=True,
                           help="Output graph file (.graphml/.pkl)")
    cr_parser.add_argument("--node-attrs", type=str, default=None,
                           help='JSON dict for new-node attributes, e.g. '
                                '\'{"color":"gray","shape":"hexagon","size":12}\'')
    cr_parser.add_argument("--edge-attrs", type=str, default=None,
                           help='JSON dict for query-result edge attributes, e.g. '
                                '\'{"color":"green","weight":8,"arrows":"to"}\'')
    cr_parser.add_argument("--disconnected-attrs", type=str, default=None,
                           help='JSON dict for disconnected-edge attributes, e.g. '
                                '\'{"color":"red","weight":2}\'')
    cr_parser.add_argument("--node-type-filter", type=str, default=None,
                           help="Filter graph query results by node type (e.g. Cell)")
    cr_parser.add_argument("--reverse-relations", type=str, default=None,
                           help="Comma-separated relations to reverse (e.g. is_a)")
    cr_parser.add_argument("--no-auto-filter", action="store_true",
                           help="Disable auto-filter for graph query results")
    cr_parser.add_argument("--no-highlight-disconnected", action="store_true",
                           help="Do not highlight disconnected original edges")

    # --- info command ---
    info_parser = subparsers.add_parser('info', help='Show graph file details')
    info_parser.add_argument("-i", "--input", required=True,
                             help="Input graph file")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(0)

    try:
        if args.command == 'read':
            G = read_graph(args.input)
            print(f"Graph loaded: {args.input}")
            print(f"  Nodes: {G.number_of_nodes()}")
            print(f"  Edges: {G.number_of_edges()}")

        elif args.command == 'convert':
            G = read_graph(args.input)
            print(f"Read graph from {args.input} "
                  f"({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")
            save_graph(G, args.output)

        elif args.command == 'data2nxg':
            print(f"Loading data from {args.input} (type: {args.type})")
            if args.type == 'df':
                data = pd.read_csv(args.input)
                if data.empty:
                    print("Error: Input DataFrame is empty", file=sys.stderr)
                    sys.exit(1)
                print(f"Loaded {len(data)} rows")
            else:  # adata or paga
                import scanpy as sc
                data = sc.read(args.input)
                print(f"Loaded AnnData: {data.n_obs} cells x {data.n_vars} genes")
            G = data2nxg(data, data_type=args.type,
                         directed=args.directed, reverse_direct=args.reverse,
                         cell_neighbor=args.cell_neighbor)
            print(f"Converted to graph: {G.number_of_nodes()} nodes, "
                  f"{G.number_of_edges()} edges")
            save_graph(G, args.output)

        elif args.command == 'combine':
            import scanpy as sc
            print(f"Loading AnnData: {args.adata}")
            adata = sc.read(args.adata)
            print(f"Loading DataFrame: {args.df}")
            df = pd.read_csv(args.df)
            print(f"Loaded {len(df)} rows")
            G = combine_graph(adata, df, directed=args.directed,
                              omics_reverse_direct=args.omics_reverse,
                              kg_reverse_direct=args.kg_reverse,
                              cell_neighbor=args.cell_neighbor)
            print(f"Combined graph: {G.number_of_nodes()} nodes, "
                  f"{G.number_of_edges()} edges")
            save_graph(G, args.output)

        elif args.command == 'complete_result':
            import json

            # Load original graph
            print(f"Loading original graph: {args.graph}")
            original_graph = read_graph(args.graph)
            print(f"  Original: {original_graph.number_of_nodes()} nodes, "
                  f"{original_graph.number_of_edges()} edges")

            # Load and auto-detect query results
            query_results = None
            if args.query:
                print(f"Loading query results: {args.query}")
                ext = os.path.splitext(args.query)[-1].lower()
                if ext == '.csv':
                    query_results = pd.read_csv(args.query)
                    print(f"  Detected DataFrame: {len(query_results)} rows")
                elif ext in ('.graphml', '.pkl', '.svg'):
                    query_results = read_graph(args.query)
                    print(f"  Detected Graph: {query_results.number_of_nodes()} nodes, "
                          f"{query_results.number_of_edges()} edges")
                else:
                    print(f"Warning: Unknown extension {ext}, attempting graph read")
                    query_results = read_graph(args.query)

            # Parse JSON options
            node_attrs = json.loads(args.node_attrs) if args.node_attrs else {}
            edge_attrs = json.loads(args.edge_attrs) if args.edge_attrs else {}
            disconnected_attrs = (
                json.loads(args.disconnected_attrs) if args.disconnected_attrs else None
            )
            reverse_relations = (
                [r.strip() for r in args.reverse_relations.split(',')]
                if args.reverse_relations else None
            )

            # Build merged graph
            G = complete_result(
                original_graph=original_graph,
                query_results=query_results,
                node_attrs=node_attrs,
                edge_attrs=edge_attrs,
                node_type_filter=args.node_type_filter,
                reverse_relations=reverse_relations,
                auto_filter=not args.no_auto_filter,
                highlight_disconnected=not args.no_highlight_disconnected,
                disconnected_edge_attrs=disconnected_attrs
            )

            print(f"Merged graph: {G.number_of_nodes()} nodes, "
                  f"{G.number_of_edges()} edges")
            save_graph(G, args.output)

        elif args.command == 'info':
            G = read_graph(args.input)
            print(f"File: {args.input}")
            print(f"  Type: {type(G).__name__}")
            print(f"  Nodes: {G.number_of_nodes()}")
            print(f"  Edges: {G.number_of_edges()}")
            if G.number_of_nodes() > 0:
                sample_nodes = list(G.nodes(data=True))[:5]
                print(f"  Sample nodes: {sample_nodes}")
            if G.number_of_edges() > 0:
                sample_edges = list(G.edges(data=True))[:5]
                print(f"  Sample edges: {sample_edges}")

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _cli_main()
