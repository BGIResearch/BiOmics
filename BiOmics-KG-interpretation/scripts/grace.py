"""GRACE graph-contrastive embedding backend for ``brick emb``.

The training design follows the BRICK GRACE implementation (two independently
augmented graph views, a GCN encoder, and a local-to-local InfoNCE objective),
while keeping the heavyweight PyTorch stack optional.  The default node
features use a seeded random projection of node identity plus node-type and
degree features so integrated omics graphs do not require an O(N^2) identity
matrix.  Use ``feature_mode="identity"`` to reproduce the upstream one-hot
feature construction on small graphs.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import networkx as nx
import numpy as np


UPSTREAM_GRACE_URL = "https://github.com/caolei2-BGI/BRICK/blob/main/embedgraph/GRACE.py"


class GraceDependencyError(ImportError):
    """Raised when the optional GRACE runtime is not installed."""


@dataclass
class GraceResult:
    """Learned node representations and reproducibility metadata."""

    embedding: np.ndarray
    losses: List[float]
    metadata: Dict[str, Any]


def build_node_features(
    nxg: nx.Graph,
    feature_mode: str = "projected_identity",
    feature_dim: int = 64,
    seed: int = 20241204,
    identity_max_nodes: int = 4096,
) -> Tuple[np.ndarray, List[Any], List[str]]:
    """Build deterministic node features in NetworkX iteration order.

    ``projected_identity`` is a scalable random projection of the upstream
    one-hot node identity, augmented with standardized degree and node-type
    indicators.  ``identity`` creates the exact dense one-hot matrix and is
    guarded because its memory use is quadratic in the number of nodes.
    """

    nodes = list(nxg.nodes)
    node_count = len(nodes)
    if node_count == 0:
        raise ValueError("GRACE requires a graph with at least one node.")

    node_types = [str(nxg.nodes[node].get("type", "Unknown")) for node in nodes]
    type_names = sorted(set(node_types))

    if feature_mode == "identity":
        if identity_max_nodes > 0 and node_count > identity_max_nodes:
            raise ValueError(
                "Dense identity features would allocate an "
                f"{node_count} x {node_count} matrix. Use "
                "--grace-feature-mode projected_identity, or explicitly raise "
                "--grace-identity-max-nodes for a small, well-provisioned run."
            )
        features = np.eye(node_count, dtype=np.float32)
        return features, nodes, type_names

    if feature_mode != "projected_identity":
        raise ValueError(
            "feature_mode must be 'projected_identity' or 'identity', got "
            f"{feature_mode!r}."
        )
    if feature_dim < 1:
        raise ValueError("feature_dim must be at least 1.")

    rng = np.random.default_rng(seed)
    projected_identity = rng.standard_normal((node_count, feature_dim)).astype(np.float32)
    projected_identity /= np.float32(math.sqrt(feature_dim))

    degrees = np.asarray([nxg.degree(node) for node in nodes], dtype=np.float32)
    log_degrees = np.log1p(degrees)
    degree_scale = float(log_degrees.std())
    if degree_scale > 0:
        log_degrees = (log_degrees - float(log_degrees.mean())) / degree_scale
    else:
        log_degrees = np.zeros_like(log_degrees)

    type_index = {name: index for index, name in enumerate(type_names)}
    type_features = np.zeros((node_count, len(type_names)), dtype=np.float32)
    for row, node_type in enumerate(node_types):
        type_features[row, type_index[node_type]] = 1.0

    features = np.concatenate(
        [projected_identity, log_degrees[:, None], type_features], axis=1
    ).astype(np.float32, copy=False)
    return features, nodes, type_names


def _load_grace_dependencies():
    try:
        import torch
        import torch.nn.functional as functional
        from torch.optim import Adam
        from torch_geometric.nn import GCNConv
        from torch_geometric.utils import from_networkx
    except ImportError as exc:
        raise GraceDependencyError(
            "GRACE requires optional dependencies. Install them with "
            "`pip install 'BiOmics-BRICK[grace]'` (packages: torch, "
            "and torch-geometric)."
        ) from exc

    return {
        "torch": torch,
        "functional": functional,
        "Adam": Adam,
        "GCNConv": GCNConv,
        "from_networkx": from_networkx,
    }


def _resolve_device(torch, requested: str) -> str:
    requested = str(requested).lower()
    if requested == "auto":
        if torch.cuda.is_available():
            return "cuda"
        mps = getattr(getattr(torch, "backends", None), "mps", None)
        if mps is not None and mps.is_available():
            return "mps"
        return "cpu"
    if requested == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA was requested for GRACE, but torch.cuda.is_available() is false.")
    if requested == "mps":
        mps = getattr(getattr(torch, "backends", None), "mps", None)
        if mps is None or not mps.is_available():
            raise ValueError("MPS was requested for GRACE, but it is not available.")
    if requested not in {"cpu", "cuda", "mps"}:
        raise ValueError("device must be one of: auto, cpu, cuda, mps.")
    return requested


def _seed_runtime(torch, seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    cudnn = getattr(getattr(torch, "backends", None), "cudnn", None)
    if cudnn is not None:
        cudnn.deterministic = True
        cudnn.benchmark = False


def _make_model_classes(dependencies):
    torch = dependencies["torch"]
    functional = dependencies["functional"]
    GCNConv = dependencies["GCNConv"]

    class GraphEncoder(torch.nn.Module):
        def __init__(self, input_dim: int, hidden_dim: int, num_layers: int):
            super().__init__()
            if num_layers < 1:
                raise ValueError("num_layers must be at least 1.")
            dimensions = [input_dim] + [hidden_dim] * num_layers
            self.layers = torch.nn.ModuleList(
                [
                    GCNConv(dimensions[index], dimensions[index + 1], cached=False)
                    for index in range(num_layers)
                ]
            )

        def forward(self, x, edge_index, edge_weight=None):
            representation = x
            for convolution in self.layers:
                representation = convolution(representation, edge_index, edge_weight)
                representation = functional.relu(representation)
            return representation

    class GraceEncoder(torch.nn.Module):
        def __init__(self, encoder, hidden_dim: int):
            super().__init__()
            self.encoder = encoder
            self.projection_in = torch.nn.Linear(hidden_dim, hidden_dim)
            self.projection_out = torch.nn.Linear(hidden_dim, hidden_dim)

        @staticmethod
        def augment(x, edge_index, edge_weight, edge_drop_rate, feature_mask_rate):
            feature_keep = torch.rand(x.size(1), device=x.device) >= feature_mask_rate
            augmented_x = x * feature_keep.to(dtype=x.dtype).unsqueeze(0)
            edge_keep = torch.rand(edge_index.size(1), device=edge_index.device) >= edge_drop_rate
            augmented_edge_index = edge_index[:, edge_keep]
            augmented_edge_weight = None
            if edge_weight is not None:
                augmented_edge_weight = edge_weight[edge_keep]
            return augmented_x, augmented_edge_index, augmented_edge_weight

        def augmented_views(
            self,
            x,
            edge_index,
            edge_weight,
            edge_drop_rate,
            feature_mask_rate,
        ):
            first = self.augment(
                x,
                edge_index,
                edge_weight,
                edge_drop_rate,
                feature_mask_rate,
            )
            second = self.augment(
                x,
                edge_index,
                edge_weight,
                edge_drop_rate,
                feature_mask_rate,
            )
            return self.encoder(*first), self.encoder(*second)

        def project(self, representation):
            return self.projection_out(functional.elu(self.projection_in(representation)))

    return GraphEncoder, GraceEncoder


def _dual_view_info_nce(first, second, temperature, functional, torch):
    """Symmetric GRACE local-to-local InfoNCE with intra-view negatives."""

    first = functional.normalize(first, p=2, dim=1)
    second = functional.normalize(second, p=2, dim=1)
    node_count = first.size(0)
    labels = torch.arange(node_count, device=first.device)

    cross_first = first @ second.T / temperature
    within_first = first @ first.T / temperature
    within_first.fill_diagonal_(float("-inf"))
    first_logits = torch.cat([cross_first, within_first], dim=1)

    cross_second = second @ first.T / temperature
    within_second = second @ second.T / temperature
    within_second.fill_diagonal_(float("-inf"))
    second_logits = torch.cat([cross_second, within_second], dim=1)

    return 0.5 * (
        functional.cross_entropy(first_logits, labels)
        + functional.cross_entropy(second_logits, labels)
    )


def train_grace(
    nxg: nx.Graph,
    *,
    hidden_dim: int = 32,
    num_layers: int = 2,
    learning_rate: float = 0.0001,
    epochs: int = 300,
    edge_drop_rate: float = 0.3,
    feature_mask_rate: float = 0.3,
    temperature: float = 0.2,
    feature_mode: str = "projected_identity",
    feature_dim: int = 64,
    identity_max_nodes: int = 4096,
    max_nodes: int = 5000,
    device: str = "auto",
    seed: int = 20241204,
    verbose: bool = True,
) -> GraceResult:
    """Train GRACE and return a node embedding aligned to ``nxg.nodes``.

    The PyGCL local-to-local InfoNCE objective forms node-pair similarity
    matrices, so ``max_nodes`` guards against unexpectedly large memory use.
    Set it to 0 only when the available memory has been checked explicitly.
    """

    node_count = nxg.number_of_nodes()
    if max_nodes > 0 and node_count > max_nodes:
        raise ValueError(
            f"The integrated graph has {node_count} nodes, above the GRACE safety "
            f"limit of {max_nodes}. Subset the analysis, increase --grace-max-nodes, "
            "or use --method umap."
        )
    if epochs < 1:
        raise ValueError("epochs must be at least 1.")
    if hidden_dim < 1:
        raise ValueError("hidden_dim must be at least 1.")
    if learning_rate <= 0:
        raise ValueError("learning_rate must be positive.")
    if temperature <= 0:
        raise ValueError("temperature must be positive.")
    for name, value in {
        "edge_drop_rate": edge_drop_rate,
        "feature_mask_rate": feature_mask_rate,
    }.items():
        if not 0 <= value < 1:
            raise ValueError(f"{name} must be in [0, 1).")

    dependencies = _load_grace_dependencies()
    torch = dependencies["torch"]
    runtime_device = _resolve_device(torch, device)
    _seed_runtime(torch, seed)

    features, ordered_nodes, type_names = build_node_features(
        nxg,
        feature_mode=feature_mode,
        feature_dim=feature_dim,
        seed=seed,
        identity_max_nodes=identity_max_nodes,
    )
    if ordered_nodes != list(nxg.nodes):
        raise RuntimeError("NetworkX node order changed while constructing GRACE features.")

    data = dependencies["from_networkx"](nxg)
    data.x = torch.from_numpy(features)
    data = data.to(runtime_device)

    GraphEncoder, GraceEncoder = _make_model_classes(dependencies)
    graph_encoder = GraphEncoder(features.shape[1], hidden_dim, num_layers).to(runtime_device)
    encoder = GraceEncoder(graph_encoder, hidden_dim).to(runtime_device)
    optimizer = dependencies["Adam"](encoder.parameters(), lr=learning_rate)

    losses: List[float] = []
    progress_interval = max(1, epochs // 10)
    for epoch in range(1, epochs + 1):
        encoder.train()
        optimizer.zero_grad()
        first_view, second_view = encoder.augmented_views(
            data.x,
            data.edge_index,
            None,
            edge_drop_rate,
            feature_mask_rate,
        )
        loss = _dual_view_info_nce(
            encoder.project(first_view),
            encoder.project(second_view),
            temperature,
            dependencies["functional"],
            torch,
        )
        loss.backward()
        optimizer.step()
        loss_value = float(loss.detach().cpu().item())
        if not math.isfinite(loss_value):
            raise RuntimeError(f"GRACE loss became non-finite at epoch {epoch}.")
        losses.append(loss_value)
        if verbose and (epoch == 1 or epoch == epochs or epoch % progress_interval == 0):
            print(f"GRACE epoch {epoch}/{epochs}: loss={loss_value:.6f}")

    encoder.eval()
    with torch.no_grad():
        embedding = encoder.encoder(data.x, data.edge_index, None)
    embedding_array = embedding.detach().cpu().numpy().astype(np.float32, copy=False)

    metadata = {
        "method": "grace",
        "backend": "pytorch_geometric",
        "upstream_source": UPSTREAM_GRACE_URL,
        "node_count": int(node_count),
        "edge_count": int(nxg.number_of_edges()),
        "input_feature_count": int(features.shape[1]),
        "node_type_count": int(len(type_names)),
        "hidden_dim": int(hidden_dim),
        "num_layers": int(num_layers),
        "learning_rate": float(learning_rate),
        "epochs": int(epochs),
        "edge_drop_rate": float(edge_drop_rate),
        "feature_mask_rate": float(feature_mask_rate),
        "temperature": float(temperature),
        "feature_mode": feature_mode,
        "feature_dim": int(feature_dim),
        "device": runtime_device,
        "seed": int(seed),
        "final_loss": float(losses[-1]),
    }
    return GraceResult(embedding=embedding_array, losses=losses, metadata=metadata)
