#!/usr/bin/env python3
"""brick pl - plotting CLI for brick knowledge graph visualization."""

import argparse
import json
import os
import pickle
import sys

import matplotlib
import networkx as nx
import pandas as pd

matplotlib.use("Agg")


def _load_graph(path):
    """Load a graph file supported by brick pl."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    try:
        if path.endswith(".graphml"):
            return nx.read_graphml(path)
        if path.endswith(".gml"):
            return nx.read_gml(path)
        if path.endswith(".gpickle"):
            with open(path, "rb") as f:
                return pickle.load(f)
        if path.endswith((".json", ".json.gz")):
            from networkx.readwrite import json_graph

            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            return json_graph.node_link_graph(data)
    except Exception as e:
        raise ValueError(f"Error loading graph file: {e}")

    raise ValueError(
        f"Unsupported file format: {path}\n"
        "Supported formats: .graphml, .gml, .gpickle, .json"
    )


def _ensure_node_names(nxg):
    """Populate missing display names with node IDs."""
    if nxg.number_of_nodes() == 0:
        return

    first_node = list(nxg.nodes())[0]
    if "name" not in nxg.nodes[first_node]:
        print("Warning: Nodes missing 'name' attribute, using node ID as label")
        for node in nxg.nodes():
            nxg.nodes[node]["name"] = str(node)


def main():
    parser = argparse.ArgumentParser(
        description="brick pl - visualization tools for brick knowledge graph",
        usage="brick pl <command> [options]",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available plot types")

    interact_parser = subparsers.add_parser(
        "interact", help="Interactive network visualization (HTML)"
    )
    interact_parser.add_argument(
        "-i", "--input", required=True, help="Input graph file (graphml, gml, etc.)"
    )
    interact_parser.add_argument(
        "-o", "--output", default="network.html", help="Output HTML file"
    )
    interact_parser.add_argument("--width", type=int, default=800, help="Canvas width")
    interact_parser.add_argument("--height", type=int, default=600, help="Canvas height")

    static_parser = subparsers.add_parser(
        "static", help="Static network visualization (PNG/PDF)"
    )
    static_parser.add_argument("-i", "--input", required=True, help="Input graph file")
    static_parser.add_argument(
        "-o", "--output", default="network.png", help="Output image file"
    )
    static_parser.add_argument(
        "--layout",
        default="spring",
        choices=["spring", "tree", "circular", "kamada_kawai", "shell", "spectral"],
        help="Layout algorithm",
    )
    static_parser.add_argument(
        "--root", help="Root node for hierarchical layout (requires --layout tree)"
    )

    radial_parser = subparsers.add_parser("radial", help="Radial enrichment plot")
    radial_parser.add_argument(
        "-i", "--input", required=True, help="Input CSV file with enrichment results"
    )
    radial_parser.add_argument(
        "-o", "--output", default="radial.png", help="Output image file"
    )
    radial_parser.add_argument(
        "--top-k", type=int, default=7, help="Top K entities to show"
    )

    tf_parser = subparsers.add_parser("tf", help="TF-regulon network visualization")
    tf_parser.add_argument(
        "-i", "--input", required=True, help="Input CSV file with TF regulon data"
    )
    tf_parser.add_argument("--tf", required=True, help="TF names, comma-separated")
    tf_parser.add_argument("--verified", help="Verified genes file, one gene per line")
    tf_parser.add_argument(
        "-o", "--output", default="tf_regulon.png", help="Output image file"
    )

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        if args.command == "interact":
            cmd_interact(args)
        elif args.command == "static":
            cmd_static(args)
        elif args.command == "radial":
            cmd_radial(args)
        elif args.command == "tf":
            cmd_tf(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def cmd_interact(args):
    """Interactive network visualization."""
    from .interact_graph import interact_visualize_network

    nxg = _load_graph(args.input)
    _ensure_node_names(nxg)

    print(
        f"Loaded graph with {nxg.number_of_nodes()} nodes and "
        f"{nxg.number_of_edges()} edges"
    )
    interact_visualize_network(
        nxg,
        width=args.width,
        height=args.height,
        notebook=False,
        save=args.output,
    )
    print(f"Saved interactive visualization to {args.output}")


def cmd_static(args):
    """Static network visualization."""
    from .static_graph import static_visualize_network

    nxg = _load_graph(args.input)
    _ensure_node_names(nxg)

    print(
        f"Loaded graph with {nxg.number_of_nodes()} nodes and "
        f"{nxg.number_of_edges()} edges"
    )

    if not isinstance(nxg, (nx.MultiGraph, nx.MultiDiGraph)):
        nxg = nx.MultiDiGraph(nxg)

    static_visualize_network(
        nxg,
        layout=args.layout,
        root=args.root,
        save_path=args.output,
    )
    print(f"Saved static visualization to {args.output}")


def cmd_radial(args):
    """Radial enrichment plot."""
    from .radial_enrich import radial_enrich_plot

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        raise ValueError(f"Error reading CSV file: {e}")

    if len(df) == 0:
        raise ValueError("Input CSV file is empty")

    print(f"Loaded {len(df)} records")

    required_cols = [
        "path.0.name",
        "path.2.name",
        "path.2.type",
        "path.2.enrich_pvalue",
    ]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Missing required columns: {missing_cols}\n"
            f"Available columns: {list(df.columns)}\n"
            "Please ensure the input is from brick rk --enrich output"
        )

    try:
        fig = radial_enrich_plot(df, top_k=args.top_k)
        fig.savefig(args.output, dpi=300, bbox_inches="tight")
        print(f"Saved radial plot to {args.output}")
    except Exception as e:
        raise RuntimeError(f"Error generating radial plot: {e}")


def cmd_tf(args):
    """TF-regulon visualization."""
    from .visulize_tf_regulon import visualize_tf_regulon

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    try:
        df = pd.read_csv(args.input)
        if len(df.columns) == 0:
            raise ValueError("Input CSV file has no columns")
        df.set_index(df.columns[0], inplace=True)
    except Exception as e:
        raise ValueError(f"Error reading CSV file: {e}")

    if len(df) == 0:
        raise ValueError("Input CSV file is empty")

    if "Regulon" not in df.columns:
        raise ValueError(
            "Input file must contain 'Regulon' column\n"
            f"Available columns: {list(df.columns)}\n"
            "Please ensure the input is a valid TF regulon file"
        )

    tf_list = [tf.strip() for tf in args.tf.split(",") if tf.strip()]
    if not tf_list:
        raise ValueError("No transcription factors specified. Use --tf to specify TF names")

    missing_tfs = [tf for tf in tf_list if tf not in df.index]
    if missing_tfs:
        raise ValueError(
            f"TFs not found in input file: {missing_tfs}\n"
            f"Available TFs: {list(df.index)}"
        )

    verified_genes = []
    if args.verified:
        if not os.path.exists(args.verified):
            raise FileNotFoundError(f"Verified genes file not found: {args.verified}")
        try:
            with open(args.verified, "r", encoding="utf-8") as f:
                verified_genes = [line.strip() for line in f if line.strip()]
        except Exception as e:
            raise ValueError(f"Error reading verified genes file: {e}")

    print(f"Visualizing TF: {tf_list}")
    print(f"Found {len(verified_genes)} verified genes")

    try:
        visualize_tf_regulon(tf_list, df, verified_genes, save_path=args.output)
        print(f"Saved TF regulon plot to {args.output}")
    except Exception as e:
        raise RuntimeError(f"Error generating TF regulon plot: {e}")


if __name__ == "__main__":
    main()
