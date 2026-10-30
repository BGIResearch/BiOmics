"""Agent_BRICK - A lightweight HTTP client for BRICK Knowledge Graph API."""

__version__ = "0.2.0"

# Lazy imports to avoid circular dependencies
__all__ = ["query_relation", "rank", "enrich", "read_graph", "save_graph", "data2nxg", "combine_graph", "embedgraph", "grace_embedding"]

def __getattr__(name):
    if name == "query_relation":
        from querygraph import query_relation
        return query_relation
    elif name == "rank":
        from rankgraph import rank
        return rank
    elif name == "enrich":
        from rankgraph import enrich
        return enrich
    elif name == "embedgraph":
        from embedgraph import UMAP_graph
        return UMAP_graph
    elif name == "grace_embedding":
        from embedgraph import GRACE_graph
        return GRACE_graph
    elif name == "read_graph":
        from graphadapter import read_graph
        return read_graph
    elif name == "save_graph":
        from graphadapter import save_graph
        return save_graph
    elif name == "data2nxg":
        from graphadapter import data2nxg
        return data2nxg
    elif name == "combine_graph":
        from graphadapter import combine_graph
        return combine_graph
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
