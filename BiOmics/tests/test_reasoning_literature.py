import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from reasoning import _parse_list_field, collect_co_mentioned_literature


def test_parse_list_field_handles_lists_strings_and_empty():
    assert _parse_list_field("['PMID:1', 'PMID:2']") == ["PMID:1", "PMID:2"]
    assert _parse_list_field("[]") == []
    assert _parse_list_field("nan") == []
    assert _parse_list_field(None) == []
    assert _parse_list_field(["PMID:3"]) == ["PMID:3"]
    assert _parse_list_field("plain text") == ["plain text"]


def test_collect_co_mentioned_literature_aligns_pmids_and_sentences():
    df = pd.DataFrame([
        {  # co-mentioned row with 2 aligned PMIDs/sentences
            "path.0.name": "TCF7L2",
            "path.1.relation": "co_mentioned",
            "path.1.info_source": "['PMID:111', 'PMID:222']",
            "path.1.original_relation": "['sentence one', 'sentence two']",
            "path.2.name": "type 2 diabetes mellitus",
        },
        {  # curated row must be ignored
            "path.0.name": "TCF7L2",
            "path.1.relation": "causally_related_to",
            "path.1.info_source": "['GWAS_Catalog']",
            "path.1.original_relation": "['curated edge']",
            "path.2.name": "type 2 diabetes mellitus",
        },
        {  # duplicate co-mentioned row must be deduped
            "path.0.name": "TCF7L2",
            "path.1.relation": "co_mentioned",
            "path.1.info_source": "['PMID:111']",
            "path.1.original_relation": "['sentence one']",
            "path.2.name": "type 2 diabetes mellitus",
        },
    ])
    entries = collect_co_mentioned_literature(df)
    assert len(entries) == 2
    assert entries[0]["pmid"] == "PMID:111"
    assert entries[0]["sentence"] == "sentence one"
    assert entries[1]["pmid"] == "PMID:222"
    assert entries[1]["sentence"] == "sentence two"
    assert entries[0]["source"] == "TCF7L2"
    assert entries[0]["target"] == "type 2 diabetes mellitus"


def test_collect_co_mentioned_literature_handles_shorter_sentence_list():
    df = pd.DataFrame([
        {
            "path.0.name": "A",
            "path.1.relation": "co_mentioned",
            "path.1.info_source": "['PMID:1', 'PMID:2']",
            "path.1.original_relation": "['only one sentence']",
            "path.2.name": "B",
        }
    ])
    entries = collect_co_mentioned_literature(df)
    assert len(entries) == 2
    assert entries[1]["sentence"] == ""
