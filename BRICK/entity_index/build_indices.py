from __future__ import annotations

import argparse
import logging
from typing import Iterable, Optional

from .config import NODE_CONFIGS
from .embedding import EmbeddingClient, EmbeddingConfig, build_batch_embedder
from .es_writer import ESConfig, index_artifacts
from .pipeline import build_index_artifacts


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Elasticsearch indices for entity catalog")
    parser.add_argument(
        "--nodes",
        nargs="*",
        default=None,
        help="Specific node types to index (defaults to all configured nodes)",
    )
    parser.add_argument(
        "--sample-limit",
        type=int,
        default=1000,
        help="Number of documents to index per channel (default: 1000; use -1 for all)",
    )
    parser.add_argument(
        "--no-recreate",
        action="store_true",
        help="Do not drop existing indices before indexing",
    )
    parser.add_argument(
        "--skip-vector",
        action="store_true",
        help="Skip building vector indices even if the node is configured for them",
    )
    parser.add_argument(
        "--embedding-dim",
        type=int,
        default=1024,
        help="Embedding dimension (default: 1024 for text-embedding-v3)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10,
        help="Embedding batch size (default: 32)",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Maximum embedding retries per batch (default: 3)",
    )
    parser.add_argument(
        "--index-prefix",
        type=str,
        default=None,
        help="Custom index prefix to avoid clashing with shared clusters (overrides ES_INDEX_PREFIX).",
    )
    parser.add_argument(
        "--catalog-sample-limit",
        type=int,
        default=None,
        help="TEMP: limit catalog records per node before embedding for quick testing; remove after validation.",
    )
    return parser.parse_args()


def _need_vectors(nodes: Optional[Iterable[str]]) -> bool:
    selected = nodes or NODE_CONFIGS.keys()
    return any(NODE_CONFIGS[node].enable_vector for node in selected)


def main() -> None:
    import os
    print("当前 ES_CONFIG", os.environ.get("ES_CONFIG"))
    args = _parse_args()
    node_names = args.nodes
    sample_limit = None if args.sample_limit == -1 else args.sample_limit

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    embedding_fn = None
    include_vectors = not args.skip_vector

    if include_vectors and _need_vectors(node_names):
        embed_config = EmbeddingConfig.from_env()
        embed_client = EmbeddingClient(embed_config)
        embedding_fn = build_batch_embedder(
            embed_client,
            batch_size=args.batch_size,
            max_retries=args.max_retries,
        )

    artifacts = build_index_artifacts(
        node_names=node_names,
        embedding_fn=embedding_fn,
        embedding_dimension=args.embedding_dim,
        include_vectors=include_vectors,
        catalog_sample_limit=args.catalog_sample_limit,
    )

    es_cfg = ESConfig.from_env()
    es_client = es_cfg.create_client()
    index_artifacts(
        es_client,
        artifacts,
        recreate=not args.no_recreate,
        sample_limit=sample_limit,
        index_prefix=args.index_prefix,
    )


if __name__ == "__main__":
    main()
