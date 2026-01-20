"""
Scripts module for Query Agent.

This module provides various query functions for interacting with
the BRICK knowledge graph system.
"""

from .Most_related_query import most_related_query as mrq
from .Relation_query import relation_query as rq
from .Terminology_query import terminology_query as tq

__all__ = [
    "mrq",
    "rq", 
    "tq",
]

__version__ = "0.1.0"
