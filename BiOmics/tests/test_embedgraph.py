import importlib.util
from pathlib import Path
import sys
import types

import numpy as np


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"


def _load_embedgraph_without_scanpy_runtime():
    module_name = "brick_embedgraph_test_module"
    sys.modules.pop(module_name, None)
    previous_scanpy = sys.modules.get("scanpy")
    if previous_scanpy is None:
        sys.modules["scanpy"] = types.ModuleType("scanpy")
    try:
        spec = importlib.util.spec_from_file_location(
            module_name,
            SCRIPTS_DIR / "embedgraph.py",
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        if previous_scanpy is None:
            sys.modules.pop("scanpy", None)
        else:
            sys.modules["scanpy"] = previous_scanpy


def test_pairwise_distance_supports_euclidean_and_cosine():
    embedgraph = _load_embedgraph_without_scanpy_runtime()
    left = np.asarray([[1.0, 0.0], [0.0, 1.0]])
    right = np.asarray([[1.0, 0.0]])

    np.testing.assert_allclose(
        embedgraph._pairwise_distance(left, right, metric="euclidean"),
        [[0.0], [np.sqrt(2.0)]],
    )
    np.testing.assert_allclose(
        embedgraph._pairwise_distance(left, right, metric="cosine"),
        [[0.0], [1.0]],
    )


def test_embed_dispatch_keeps_legacy_method_isolated_from_grace_options(monkeypatch):
    embedgraph = _load_embedgraph_without_scanpy_runtime()
    captured = {}

    def fake_umap(adata_path, query_path, **kwargs):
        captured.update(kwargs)
        return adata_path, query_path

    monkeypatch.setattr(embedgraph, "UMAP_graph", fake_umap)
    result = embedgraph.embed_graph(
        "input.h5ad",
        "query.csv",
        method="umap",
        cell_neighbor=True,
        distance_metric="euclidean",
        epochs=999,
    )

    assert result == ("input.h5ad", "query.csv")
    assert captured == {"cell_neighbor": True, "distance_metric": "euclidean"}
