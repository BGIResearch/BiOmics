import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from evidence_packet import build_evidence_packet, build_optionwise_evidence_packets, classify_evidence_row


def test_evidence_packet_separates_strong_weak_and_background():
    rows = [
        {
            "source": "DrugA",
            "source_type": "Chemical",
            "relation": "pharmaco_effect_on",
            "target": "DiseaseX",
            "target_type": "Disease",
            "info_source": "DrugBank",
            "evidence_text": "DrugA has a pharmacological effect on DiseaseX.",
        },
        {
            "source": "DrugB",
            "source_type": "Chemical",
            "relation": "co_mentioned",
            "target": "DiseaseX",
            "target_type": "Disease",
            "info_source": "PubMed text-mined",
            "evidence_text": "DrugB and DiseaseX appeared in one abstract.",
        },
        {
            "source": "DrugC",
            "source_type": "Chemical",
            "relation": "pharmaco_effect_on",
            "target": "ProteinY",
            "target_type": "Protein",
            "info_source": "DrugBank",
            "evidence_text": "DrugC affects ProteinY.",
        },
        {
            "source": "DrugD",
            "source_type": "Chemical",
            "relation": "no_relation",
            "target": "DiseaseX",
            "target_type": "Disease",
        },
    ]
    packet = build_evidence_packet(rows, expected_answer_type="Disease", question_scope_terms=["DiseaseX"])

    assert [row["source"] for row in packet["strong_curated"]] == ["DrugA"]
    assert {row["source"] for row in packet["weak_trace"]} == {"DrugB"}
    assert {row["source"] for row in packet["drop"]} == {"DrugD"}
    assert [row["source"] for row in packet["background"]] == ["DrugC"]
    assert packet["strong_curated"][0]["answer_role"] == "candidate_answer"
    assert all(not row["candidate_eligible"] for row in packet["weak_trace"])
    assert packet["diagnostics"]["has_strong_evidence"] is True


def test_no_relation_and_co_mentioned_are_never_strong_candidates():
    for relation in ["co_mentioned", "mentioned"]:
        row = {
            "source": "OptionA",
            "source_type": "Gene",
            "relation": relation,
            "target": "DiseaseX",
            "target_type": "Disease",
        }
        classification = classify_evidence_row(
            row,
            expected_answer_type="Gene",
            candidate_terms=["OptionA"],
            candidate_side="source",
        )
        assert classification["evidence_group"] == "weak_trace"
        assert classification["candidate_eligible"] is False
    no_relation = classify_evidence_row(
        {
            "source": "OptionA",
            "source_type": "Gene",
            "relation": "no_relation",
            "target": "DiseaseX",
            "target_type": "Disease",
        },
        expected_answer_type="Gene",
        candidate_terms=["OptionA"],
        candidate_side="source",
    )
    assert no_relation["evidence_group"] == "drop"
    assert no_relation["candidate_eligible"] is False


def test_optionwise_packets_support_mcq_direct_only_fallback():
    option_rows = {
        "TERT": pd.DataFrame([
            {
                "path.0.name": "TERT",
                "path.0.type": "Gene",
                "path.1.relation": "causally_related_to",
                "path.1.info_source": "ClinVar",
                "path.1.original_relation": "TERT is causally related to DiseaseX.",
                "path.2.name": "DiseaseX",
                "path.2.type": "Disease",
            }
        ]),
        "HLA-B": pd.DataFrame([
            {
                "path.0.name": "HLA-B",
                "path.0.type": "Gene",
                "path.1.relation": "co_mentioned",
                "path.1.info_source": "PubMed text-mined",
                "path.2.name": "DiseaseX",
                "path.2.type": "Disease",
            }
        ]),
    }
    result = build_optionwise_evidence_packets(
        option_rows,
        expected_answer_type="Gene",
        question_scope_terms=["DiseaseX"],
    )
    assert result["fallback_native"] is False
    assert result["option_strong_counts"]["TERT"] == 1
    assert result["option_strong_counts"]["HLA-B"] == 0
    assert result["option_weak_counts"]["HLA-B"] == 1
    assert result["prompt_groups"] == ["strong_curated"]

    weak_only = build_optionwise_evidence_packets(
        {"HLA-B": option_rows["HLA-B"]},
        expected_answer_type="Gene",
        question_scope_terms=["DiseaseX"],
    )
    assert weak_only["fallback_native"] is True
    assert weak_only["prompt_groups"] == []


def test_snippet_supported_co_mentioned_is_capped_to_factoid_or_list_only():
    rows = [
        {
            "source": "DrugA",
            "source_type": "Chemical",
            "relation": "co_mentioned",
            "target": "DiseaseX",
            "target_type": "Disease",
            "info_source": "PubMed",
            "evidence_text": "DrugA treatment was associated with improved DiseaseX response.",
        },
        {
            "source": "DrugB",
            "source_type": "Chemical",
            "relation": "co_mentioned",
            "target": "DiseaseX",
            "target_type": "Disease",
            "info_source": "PubMed",
            "evidence_text": "DrugB and DiseaseX appeared in one abstract.",
        },
    ]
    packet = build_evidence_packet(
        rows,
        expected_answer_type="Disease",
        question_scope_terms=["DrugA"],
        candidate_side="target",
        task_type="ShortEntity",
        weak_evidence_policy="snippet_supported",
        allow_weak_final_types="factoid,list",
    )
    assert packet["diagnostics"]["snippet_supported_weak_rows"] == 1
    assert packet["diagnostics"]["weak_trace_rows"] == 1
    assert packet["snippet_supported_weak"][0]["candidate_eligible"] is True
    assert packet["snippet_supported_weak"][0]["evidence_tier"] == "E2"

    mcq_packet = build_evidence_packet(
        rows[:1],
        expected_answer_type="Chemical",
        question_scope_terms=["DiseaseX"],
        candidate_terms=["DrugA"],
        candidate_side="source",
        task_type="ClosedChoice",
        weak_evidence_policy="snippet_supported",
        allow_weak_final_types="factoid,list",
    )
    assert mcq_packet["diagnostics"]["snippet_supported_weak_rows"] == 1
    assert mcq_packet["snippet_supported_weak"][0]["candidate_eligible"] is False


def test_optionwise_scope_match_all_requires_all_question_targets():
    option_rows = {
        "PTPN22": pd.DataFrame([
            {
                "path.0.name": "PTPN22",
                "path.0.type": "Gene",
                "path.1.relation": "functionally_related_to",
                "path.1.info_source": "AGRKB",
                "path.1.original_relation": "PTPN22 is functionally related to head and neck cancer.",
                "path.2.name": "head and neck cancer",
                "path.2.type": "Disease",
            }
        ]),
        "TERT": pd.DataFrame([
            {
                "path.0.name": "TERT",
                "path.0.type": "Gene",
                "path.1.relation": "causally_related_to",
                "path.1.info_source": "ClinVar",
                "path.1.original_relation": "TERT is causally related to head and neck cancer.",
                "path.2.name": "head and neck cancer",
                "path.2.type": "Disease",
            },
            {
                "path.0.name": "TERT",
                "path.0.type": "Gene",
                "path.1.relation": "causally_related_to",
                "path.1.info_source": "ClinVar",
                "path.1.original_relation": "TERT is causally related to uveal melanoma.",
                "path.2.name": "uveal melanoma",
                "path.2.type": "Disease",
            },
        ]),
        "HLA-B": pd.DataFrame([
            {
                "path.0.name": "HLA-B",
                "path.0.type": "Gene",
                "path.1.relation": "functionally_related_to",
                "path.1.info_source": "AGRKB",
                "path.1.original_relation": "HLA-B is related to head and neck cancer and uveal melanoma.",
                "path.2.name": "uveal melanoma",
                "path.2.type": "Disease",
            }
        ]),
    }
    result = build_optionwise_evidence_packets(
        option_rows,
        expected_answer_type="Gene",
        question_scope_terms=["head and neck cancer", "uveal melanoma"],
        question_scope_match="all",
    )
    assert result["option_strong_counts"]["PTPN22"] == 0
    assert result["option_partial_strong_counts"]["PTPN22"] == 1
    assert result["option_scope_pass"]["PTPN22"] is False
    assert result["option_strong_counts"]["TERT"] == 2
    assert result["option_scope_pass"]["TERT"] is True
    assert result["option_covered_targets"]["TERT"] == ["head and neck cancer", "uveal melanoma"]
    assert result["option_strong_counts"]["HLA-B"] == 1
    assert result["option_scope_pass"]["HLA-B"] is True


def test_optionwise_stage_policy_demotes_drop_target_candidates():
    option_rows = {
        "TERT": pd.DataFrame([
            {
                "path.0.name": "TERT",
                "path.0.type": "Gene",
                "path.1.relation": "causally_related_to",
                "path.1.info_source": "ClinVar",
                "path.1.original_relation": "TERT is causally related to Crohn's disease.",
                "path.2.name": "Crohn's disease",
                "path.2.type": "Disease",
                "query.stage": "fallback_drop_target_no_relation",
            }
        ]),
        "IL23R": pd.DataFrame([
            {
                "path.0.name": "IL23R",
                "path.0.type": "Gene",
                "path.1.relation": "disease_relation",
                "path.1.info_source": "OMIM",
                "path.1.original_relation": "IL23R is related to Crohn's disease.",
                "path.2.name": "Crohn's disease",
                "path.2.type": "Disease",
                "query.stage": "initial",
            }
        ]),
    }
    result = build_optionwise_evidence_packets(
        option_rows,
        expected_answer_type="Gene",
        question_scope_terms=["Crohn's disease"],
        candidate_stage_policy="target_retained_only",
    )
    assert result["option_strong_counts"]["TERT"] == 0
    assert result["option_background_counts"]["TERT"] == 1
    assert "stage is not candidate-eligible" in result["packets"]["TERT"]["background"][0]["boundary_reason"]
    assert result["option_strong_counts"]["IL23R"] == 1


def test_optionwise_ranking_margin_controls_prompt_options():
    option_rows = {
        "TERT": pd.DataFrame([
            {
                "path.0.name": "TERT",
                "path.0.type": "Gene",
                "path.1.relation": "causally_related_to",
                "path.1.info_source": "ClinVar",
                "path.1.original_relation": "TERT is causally related to lung carcinoma and lung adenocarcinoma.",
                "path.2.name": "lung carcinoma",
                "path.2.type": "Disease",
                "query.stage": "initial",
            }
        ]),
        "TP53": pd.DataFrame([
            {
                "path.0.name": "TP53",
                "path.0.type": "Gene",
                "path.1.relation": "functionally_related_to",
                "path.1.info_source": "AGRKB",
                "path.1.original_relation": "TP53 is functionally related to lung carcinoma and lung adenocarcinoma.",
                "path.2.name": "lung carcinoma",
                "path.2.type": "Disease",
                "query.stage": "fallback_relax_target_type",
            }
        ]),
    }
    result = build_optionwise_evidence_packets(
        option_rows,
        expected_answer_type="Gene",
        question_scope_terms=["lung carcinoma", "lung adenocarcinoma"],
        question_scope_match="all",
        scope_aggregation="option",
        evidence_ranking=True,
        evidence_margin=0.15,
    )
    assert result["fallback_native"] is False
    assert result["prompt_options"] == ["TERT"]
    assert result["best_option"] == "TERT"

    blocked = build_optionwise_evidence_packets(
        option_rows,
        expected_answer_type="Gene",
        question_scope_terms=["lung carcinoma", "lung adenocarcinoma"],
        question_scope_match="all",
        scope_aggregation="option",
        evidence_ranking=True,
        evidence_margin=0.9,
    )
    assert blocked["fallback_native"] is True
    assert blocked["prompt_options"] == []


def test_reason_packet_cli_outputs_boundary_groups(tmp_path):
    path = tmp_path / "kg.csv"
    pd.DataFrame([
        {
            "path.0.name": "DrugA",
            "path.0.type": "Chemical",
            "path.1.relation": "pharmaco_effect_on",
            "path.1.info_source": "DrugBank",
            "path.1.original_relation": "DrugA treats DiseaseX.",
            "path.2.name": "DiseaseX",
            "path.2.type": "Disease",
        },
        {
            "path.0.name": "DrugB",
            "path.0.type": "Chemical",
            "path.1.relation": "co_mentioned",
            "path.1.info_source": "PubMed text-mined",
            "path.2.name": "DiseaseX",
            "path.2.type": "Disease",
        },
    ]).to_csv(path, index=False)

    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "reasoning.py"),
            "packet",
            "-i",
            str(path),
            "--expected-answer-type",
            "Disease",
            "--scope-term",
            "DiseaseX",
        ],
        cwd=str(ROOT),
        text=True,
        check=True,
        capture_output=True,
    )
    packet = json.loads(result.stdout)
    assert packet["diagnostics"]["strong_curated_rows"] == 1
    assert packet["diagnostics"]["weak_trace_rows"] == 1
