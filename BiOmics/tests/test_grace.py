from pathlib import Path
import sys

import networkx as nx
import numpy as np
import pytest


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from grace import build_node_features, train_grace


def _small_integrated_graph():
    graph = nx.Graph()
    graph.add_node("cell-1", type="Omics_Cell")
    graph.add_node("G1", type="Gene")
    graph.add_node("G2", type="Gene")
    graph.add_node("candidate", type="KG_Cell")
    graph.add_edges_from(
        [("cell-1", "G1"), ("cell-1", "G2"), ("G1", "candidate")]
    )
    return graph


def test_projected_identity_features_are_deterministic_and_typed():
    graph = _small_integrated_graph()
    first, nodes, types = build_node_features(graph, feature_dim=8, seed=17)
    second, second_nodes, second_types = build_node_features(graph, feature_dim=8, seed=17)

    assert nodes == list(graph.nodes)
    assert second_nodes == nodes
    assert types == ["Gene", "KG_Cell", "Omics_Cell"]
    assert second_types == types
    assert first.shape == (4, 8 + 1 + len(types))
    np.testing.assert_allclose(first, second)
    np.testing.assert_allclose(first[:, -len(types):].sum(axis=1), 1.0)


def test_identity_features_match_upstream_small_graph_behavior():
    features, nodes, _ = build_node_features(
        _small_integrated_graph(),
        feature_mode="identity",
        identity_max_nodes=10,
    )

    assert nodes == list(_small_integrated_graph().nodes)
    np.testing.assert_array_equal(features, np.eye(4, dtype=np.float32))


def test_identity_features_have_quadratic_memory_guard():
    with pytest.raises(ValueError, match="Dense identity features"):
        build_node_features(
            _small_integrated_graph(),
            feature_mode="identity",
            identity_max_nodes=3,
        )

def test_grace_node_safety_limit_fails_before_optional_imports():
    with pytest.raises(ValueError, match="above the GRACE safety limit"):
        train_grace(_small_integrated_graph(), max_nodes=3, epochs=1, verbose=False)


def test_grace_hyperparameters_are_validated_before_optional_imports():
    with pytest.raises(ValueError, match="edge_drop_rate"):
        train_grace(
            _small_integrated_graph(),
            edge_drop_rate=1.0,
            epochs=1,
            verbose=False,
        )
