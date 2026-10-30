import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from reasoning import normalize_evidence_df, score_evidence_df
from rankgraph import infer_terminal_node_prefix, rank


def test_multihop_normalization_preserves_edge_metadata_and_weak_co_mentions():
    df = pd.DataFrame([
        {
            "path.0.name": "GeneA",
            "path.0.type": "Gene",
            "path.1.relation": "co_mentioned",
            "path.1.info_source": "PubMed text-mined",
            "path.1.info_source_length": 1,
            "path.1.relation_confidence": 0.9,
            "path.1.original_relation": "GeneA and ProteinB appeared in the same abstract",
            "path.1.condition": "disease context",
            "path.2.name": "ProteinB",
            "path.2.type": "Protein",
            "path.3.relation": "participate_in",
            "path.3.info_source": "Reactome",
            "path.3.relation_confidence": 0.95,
            "path.3.original_relation": "ProteinB participates in PathwayC",
            "path.4.name": "PathwayC",
            "path.4.type": "Pathway",
        },
        {
            "path.0.name": "GeneA",
            "path.0.type": "Gene",
            "path.1.relation": "pharmaco_effect_on",
            "path.1.info_source": "DrugBank",
            "path.1.relation_confidence": 0.92,
            "path.1.original_relation": "GeneA targeted by DrugD",
            "path.2.name": "DrugD",
            "path.2.type": "Chemical",
        },
    ])
    evidence = normalize_evidence_df(df, provenance_file="unit.csv")
    assert len(evidence) == 3
    assert {"path_hops", "edge_position", "original_relation", "condition", "evidence_group", "answer_role"}.issubset(evidence.columns)
    weak = evidence[evidence["relation"].eq("co_mentioned")].iloc[0]
    assert weak["manual_review_required"] is True or str(weak["manual_review_required"]).lower() == "true"
    assert str(weak["evidence_strength"]).startswith("weak")
    assert weak["evidence_group"] == "weak_trace"
    assert weak["answer_role"] == "weak_trace"

    scored = score_evidence_df(evidence)
    weak_score = scored.loc[scored["relation"].astype(str).str.contains("co_mentioned"), "evidence_score"].max()
    assert weak_score <= 0.40


def test_mixed_onehop_multihop_rows_have_row_specific_path_hops():
    df = pd.DataFrame([
        {
            "path.0.name": "GeneA",
            "path.0.type": "Gene",
            "path.1.relation": "pharmaco_effect_on",
            "path.1.info_source": "DrugBank",
            "path.1.original_relation": "GeneA targeted by DrugD",
            "path.2.name": "DrugD",
            "path.2.type": "Chemical",
        },
        {
            "path.0.name": "GeneA",
            "path.0.type": "Gene",
            "path.1.relation": "interact_with",
            "path.1.info_source": "BioGRID_hs",
            "path.1.original_relation": "GeneA interacts with ProteinB",
            "path.2.name": "ProteinB",
            "path.2.type": "Protein",
            "path.3.relation": "participate_in",
            "path.3.info_source": "Reactome",
            "path.3.original_relation": "ProteinB participates in PathwayC",
            "path.4.name": "PathwayC",
            "path.4.type": "Pathway",
        },
    ])
    evidence = normalize_evidence_df(df, provenance_file="mixed.csv")
    onehop = evidence[evidence["target"].eq("DrugD")].iloc[0]
    assert onehop["path_hops"] == 1
    assert onehop["evidence_kind"] == "kg_edge"
    path_edges = evidence[evidence["path_hops"].eq(2)]
    assert len(path_edges) == 2


def test_inhibition_is_effect_direction_not_negative_evidence():
    evidence = pd.DataFrame([
        {
            "source": "DrugA",
            "relation": "pharmaco_effect_on",
            "target": "EGFR",
            "confidence": 0.9,
            "info_source": "DrugBank",
            "info_source_count": 1,
            "evidence_text": "DrugA inhibits EGFR kinase activity.",
            "evidence_strength": "direct",
        }
    ])
    scored = score_evidence_df(evidence)
    row = scored.iloc[0]
    assert row["evidence_polarity"] == "positive"
    assert row["effect_direction"] == "inhibitory"


def test_rankgraph_auto_uses_terminal_multihop_node():
    df = pd.DataFrame([
        {
            "path.0.name": "GeneA",
            "path.0.id": "G1",
            "path.0.type": "Gene",
            "path.1.relation": "interact_with",
            "path.1.info_source": "BioGRID_hs",
            "path.1.info_source_length": 1,
            "path.1.relation_confidence": 0.8,
            "path.2.name": "ProteinB",
            "path.2.id": "P1",
            "path.2.type": "Protein",
            "path.3.relation": "participate_in",
            "path.3.info_source": "Reactome",
            "path.3.info_source_length": 1,
            "path.3.relation_confidence": 0.9,
            "path.4.name": "PathwayC",
            "path.4.id": "PW1",
            "path.4.type": "Pathway",
        }
    ])
    assert infer_terminal_node_prefix(df) == "path.4"
    ranked = rank(df, metrics=["match_count"], target_index="auto")
    assert ranked.iloc[0]["path.2.name"] == "PathwayC"


def test_explain_reports_global_strength_summary_and_co_mentioned_only_warning(tmp_path):
    evidence = pd.DataFrame([
        {
            "evidence_id": "e1",
            "source": "GeneA",
            "relation": "co_mentioned",
            "target": "DiseaseX",
            "confidence": 0.6,
            "info_source": "PubMed text-mined",
            "info_source_count": 1,
            "original_relation": "GeneA and DiseaseX appeared in the same sentence.",
            "evidence_text": "GeneA and DiseaseX appeared in the same sentence.",
        }
    ])
    path = tmp_path / "evidence.csv"
    evidence.to_csv(path, index=False)
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "reasoning.py"),
            "explain",
            "-e",
            str(path),
            "--target",
            "DiseaseX",
        ],
        cwd=str(ROOT),
        text=True,
        check=True,
        capture_output=True,
    )
    assert "## Evidence strength summary" in result.stdout
    assert "weak_co_mentioned_review_required" in result.stdout
    assert "Review required: this claim is supported only by co-mentioned/text-mined evidence." in result.stdout
