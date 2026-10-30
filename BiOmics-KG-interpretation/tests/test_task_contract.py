import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import biohopr_paths
from biohopr_eval import build_runner_command
from evaluation_adapters import exact_match, selected_strong_overlap
from reasoning import reason_with_contract
from renderers import render_answer
from task_contract import infer_bioasq_contract, infer_biomix_mcq_contract, infer_task_contract


def test_biomix_mcq_maps_to_closed_choice_contract():
    contract = infer_biomix_mcq_contract({
        "id": 1,
        "question": (
            "Out of the given list, which Gene is associated with head and neck cancer "
            "and uveal melanoma. Given list is: ABO, CACNA2D1, PSCA, TERT, SULT1B1"
        ),
    })
    assert contract.task_type == "ClosedChoice"
    assert contract.expected_answer_type == "Gene"
    assert contract.target_mentions == ["head and neck cancer", "uveal melanoma"]
    assert contract.options == ["ABO", "CACNA2D1", "PSCA", "TERT", "SULT1B1"]
    assert contract.evidence_policy["option_anchor"] == "hard"


def test_bioasq_types_map_to_general_contracts():
    assert infer_bioasq_contract({"type": "yesno", "body": "Is X true?"}).task_type == "BinaryClaim"
    assert infer_bioasq_contract({"type": "factoid", "body": "Which gene?"}).task_type == "ShortEntity"
    assert infer_bioasq_contract({"type": "list", "body": "List genes."}).task_type == "EntityList"
    assert infer_bioasq_contract({"type": "summary", "body": "Summarize."}).task_type == "MechanisticSummary"


def test_reason_closed_choice_uses_ranked_prompt_option():
    contract = infer_task_contract(
        {
            "dataset": "biomix",
            "task": "mcq",
            "question": "Out of the given list, which Gene is associated with disease X. Given list is: A, B",
        }
    )
    answer = reason_with_contract(
        contract,
        {
            "fallback_native": False,
            "best_option": "B",
            "prompt_options": ["B"],
            "evidence_score_margin": 0.3,
        },
        native_answer="A",
    )
    assert answer["answer"] == "B"
    assert answer["reasoning_policy"] == "closed_choice_evidence_rank"


def test_renderers_and_exact_match_follow_contract_shape():
    contract = infer_bioasq_contract({"type": "list", "body": "List genes."})
    rendered = render_answer(["TP53", "TERT"], contract)
    assert rendered == [["TP53"], ["TERT"]]
    result = exact_match(rendered, [["TP53"], ["TERT"]], contract)
    assert result["correct"] is True


def test_selected_strong_overlap_counts_override_outcomes():
    rows = [
        {"id": 1, "answer": "A", "kg": {"selected_option_has_strong_evidence": True}},
        {"id": 2, "answer": "B", "kg": {"selected_option_has_strong_evidence": True}},
        {"id": 3, "answer": "C", "kg": {"selected_option_has_strong_evidence": False}},
    ]
    overlap = selected_strong_overlap(rows, native_answers={1: "A", 2: "A"}, gold_answers={1: "A", 2: "B"})
    assert overlap["native_correct_evidence_correct"] == 1
    assert overlap["native_wrong_evidence_correct"] == 1
    assert overlap["selected_strong_questions"] == 2


def test_biohopr_question_parser_extracts_path_contract():
    contract = infer_task_contract(
        {
            "dataset": "biohopr",
            "task": "2hop_single",
            "question": "Name a drug that can treat a disease that is associated with gene/protein TNXB.",
        }
    )
    assert contract.task_type == "PathAnswerSet"
    assert contract.source_mentions == ["TNXB"]
    assert contract.expected_answer_type == "Chemical"
    assert contract.evidence_policy["source_type"] == "Gene"
    assert contract.evidence_policy["target_type"] == "Chemical"
    assert contract.evidence_policy["intermediate_type"] == "Disease"
    assert contract.metadata["bridge_entity"] is None


def test_biohopr_balanced_policy_is_extracted_question_only():
    balanced = biohopr_paths.biohopr_benchmark_policy("balanced")
    strict = biohopr_paths.biohopr_benchmark_policy("strict_e2e_question_only")
    assert balanced["query_source"] == "extracted"
    assert strict["query_source"] == "extracted"
    assert biohopr_paths.biohopr_benchmark_policy("extracted")["query_source"] == "extracted"


def test_biohopr_surface_entity_candidates_are_question_only_aliases():
    candidates = biohopr_paths.biohopr_entity_candidates(
        "IL-6",
        "Gene",
        mode="surface",
        max_candidates=8,
    )
    assert "IL-6" in candidates
    assert "IL6" in candidates
    assert "interleukin 6" in candidates


def test_biohopr_type_fallback_filter_keeps_equivalent_types():
    contract = infer_task_contract(
        {
            "dataset": "biohopr",
            "task": "2hop_single",
            "question": "Name a disease that is related to a gene/protein that is associated with drug Aspirin.",
        }
    )
    df = pd.DataFrame(
        [
            {
                "path.0.name": "Aspirin",
                "path.0.type": "Chemical",
                "path.1.relation": "pharmaco_effect_on",
                "path.2.name": "PTGS1",
                "path.2.type": "Protein",
                "path.3.relation": "causally_related_to",
                "path.4.name": "pain",
                "path.4.type": "Phenotype",
            }
        ]
    )
    assert biohopr_paths.filter_path_rows(df, contract).empty
    assert len(biohopr_paths.filter_path_rows(df, contract, allow_type_fallback=True)) == 1


def test_biohopr_scorecutoff_policy_remains_explicit_ablation():
    scorecutoff = biohopr_paths.biohopr_benchmark_policy("strict_e2e_entitylink_fallback_scorecutoff")
    assert scorecutoff["query_source"] == "extracted"
    assert scorecutoff["entity_linking"] == "surface"
    assert scorecutoff["type_fallback"] is True
    assert scorecutoff["discover_bridge_fallback"] is True
    assert scorecutoff["relax_source_type_fallback"] is False
    assert scorecutoff["relax_target_type_fallback"] is False
    assert scorecutoff["answer_selection_mode"] == "score_cutoff"
    assert scorecutoff["kg_add_max_2hop"] <= 8


def test_biohopr_optimized_policy_uses_soft_rerank_not_hard_cutoff():
    optimized = biohopr_paths.biohopr_benchmark_policy("optimized")
    assert optimized["query_source"] == "extracted"
    assert optimized["hop_aware"] is True
    assert optimized["entity_linking_1hop"] == "off"
    assert optimized["entity_linking_2hop"] == "surface"
    assert optimized["type_fallback_1hop"] is False
    assert optimized["type_fallback_2hop"] is True
    assert optimized["discover_bridge_fallback_1hop"] is False
    assert optimized["discover_bridge_fallback_2hop"] is True
    assert optimized["answer_selection_mode"] == "top_k"
    assert optimized["max_answers"] == 40
    assert optimized["kg_add_max_2hop"] == 30
    assert optimized["kg_add_max_weak_2hop"] == 6


def test_biohopr_eval_balanced_command_does_not_use_gold_hops(tmp_path):
    runner = tmp_path / "run_agent_pipeline" / "brickkg" / "run_biohopr_graph_answers.py"
    runner.parent.mkdir(parents=True)
    runner.write_text("", encoding="utf-8")
    args = argparse.Namespace(
        preset="balanced",
        benchmark_root=tmp_path,
        method=None,
        split_group="single",
        splits=None,
        workers=1,
        limit=None,
        kg_request_timeout=1.0,
        native_answer_dir=None,
        gold_dir=None,
        overwrite=False,
        no_evaluate=True,
        dry_run=True,
    )
    cmd, manifest = build_runner_command(args)
    assert manifest["method"] == "strict_e2e_question_only"
    assert cmd[cmd.index("--query-source") + 1] == "extracted"
    assert "gold_hops" not in cmd


def test_biohopr_eval_optimized_command_passes_repair_flags(tmp_path):
    runner = tmp_path / "run_agent_pipeline" / "brickkg" / "run_biohopr_graph_answers.py"
    runner.parent.mkdir(parents=True)
    runner.write_text("", encoding="utf-8")
    args = argparse.Namespace(
        preset="strict_e2e_entitylink_fallback_scorecutoff",
        benchmark_root=tmp_path,
        method=None,
        split_group="single",
        splits=None,
        workers=1,
        limit=None,
        kg_request_timeout=1.0,
        native_answer_dir=None,
        gold_dir=None,
        overwrite=False,
        no_evaluate=True,
        dry_run=True,
    )
    cmd, manifest = build_runner_command(args)
    assert manifest["policy"]["query_source"] == "extracted"
    assert cmd[cmd.index("--entity-linking") + 1] == "surface"
    assert cmd[cmd.index("--kg-add-max-weak-2hop") + 1] == "0"
    assert "--type-fallback" in cmd
    assert "--discover-bridge-fallback" in cmd
    assert "--relax-source-type-fallback" not in cmd
    assert "--relax-target-type-fallback" not in cmd
    assert "gold_hops" not in cmd


def test_biohopr_eval_hop_aware_command_passes_split_specific_repair_flags(tmp_path):
    runner = tmp_path / "run_agent_pipeline" / "brickkg" / "run_biohopr_graph_answers.py"
    runner.parent.mkdir(parents=True)
    runner.write_text("", encoding="utf-8")
    args = argparse.Namespace(
        preset="strict_e2e_hop_aware",
        benchmark_root=tmp_path,
        method=None,
        split_group="single",
        splits=None,
        workers=1,
        limit=None,
        kg_request_timeout=1.0,
        native_answer_dir=None,
        gold_dir=None,
        overwrite=False,
        no_evaluate=True,
        dry_run=True,
    )
    cmd, manifest = build_runner_command(args)
    assert manifest["policy"]["query_source"] == "extracted"
    assert cmd[cmd.index("--query-source") + 1] == "extracted"
    assert cmd[cmd.index("--entity-linking") + 1] == "off"
    assert cmd[cmd.index("--entity-linking-1hop") + 1] == "off"
    assert cmd[cmd.index("--entity-linking-2hop") + 1] == "surface"
    assert "--type-fallback" not in cmd
    assert "--discover-bridge-fallback" not in cmd
    assert "--type-fallback-2hop" in cmd
    assert "--discover-bridge-fallback-2hop" in cmd
    assert "--type-fallback-1hop" not in cmd
    assert "--discover-bridge-fallback-1hop" not in cmd
    assert cmd[cmd.index("--kg-add-max-weak-2hop") + 1] == "6"
    assert "gold_hops" not in cmd
