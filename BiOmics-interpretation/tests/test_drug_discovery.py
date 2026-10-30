import ast
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import networkx as nx
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

try:
    from drugdiscovery import chemical_scope_table, cmd_build_graph, freeze_ranking
except ImportError:  # Local pre-deployment validation copy.
    from remote_drugdiscovery import chemical_scope_table, cmd_build_graph, freeze_ranking
from agent_checker import validate_agent_plan
from rankgraph import rank


def _query_rows():
    return pd.DataFrame(
        [
            {
                "path.0.id": "G1",
                "path.0.name": "GENE1",
                "path.0.type": "Gene",
                "path.1.relation": "pharmaco_effect_on",
                "path.1.info_source": "DrugBank",
                "path.1.info_source_length": 1,
                "path.2.id": "CHEM:1",
                "path.2.name": "Shared label",
                "path.2.type": "Chemical",
            },
            {
                "path.0.id": "G2",
                "path.0.name": "GENE2",
                "path.0.type": "Gene",
                "path.1.relation": "interact_with",
                "path.1.info_source": "CHEBI",
                "path.1.info_source_length": 1,
                "path.2.id": "CHEM:2",
                "path.2.name": "Shared label",
                "path.2.type": "Chemical",
            },
        ]
    )


def _run(*args, check=True):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SCRIPTS) + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, *map(str, args)],
        cwd=str(ROOT),
        env=env,
        text=True,
        capture_output=True,
        check=check,
    )


def test_rankgraph_keeps_same_label_different_stable_ids_separate():
    result = rank(_query_rows(), metrics=["match_count"])
    assert set(result["path.2.id"]) == {"CHEM:1", "CHEM:2"}
    assert len(result) == 2


@pytest.mark.skip(reason='requires batch-2 merge: wangbiao rankgraph stable-ID / brickembed')
def test_build_graph_qualifies_colliding_labels_by_stable_id(monkeypatch, tmp_path):
    class FakeAnnData:
        var_names = pd.Index(["GENE1"])
        obs_names = pd.Index(["CELL1"])

        def __getitem__(self, _item):
            return self

        def copy(self):
            return self

    enrichment = pd.DataFrame(
        [
            {
                "path.0.id": "['GENE:1']",
                "path.0.name": "['GENE1']",
                "path.0.type": "Gene",
                "path.2.id": "DRUGBANK:DB01593",
                "path.2.name": "Zinc",
                "path.2.type": "Chemical",
                "path.2.enrich_pvalue": 0.001,
            },
            {
                "path.0.id": "['GENE:1']",
                "path.0.name": "['GENE1']",
                "path.0.type": "Gene",
                "path.2.id": "MESH:D015032",
                "path.2.name": "Zinc",
                "path.2.type": "Chemical",
                "path.2.enrich_pvalue": 0.002,
            },
            {
                "path.0.id": "['GENE:1']",
                "path.0.name": "['GENE1']",
                "path.0.type": "Gene",
                "path.2.id": "MONDO:0100096",
                "path.2.name": "COVID-19",
                "path.2.type": "Disease",
                "path.2.enrich_pvalue": 0.003,
            },
        ]
    )

    captured = {}

    def fake_combine_graph(_adata, frame):
        captured["frame"] = frame.copy()
        graph = nx.Graph()
        for _, row in frame.iterrows():
            for source in ast.literal_eval(row["path.0.name"]):
                graph.add_node(source, type="Gene")
                graph.add_edge(source, str(row["path.2.name"]))
            graph.add_node(str(row["path.2.name"]), type=f"KG_{row['path.2.type']}")
        return graph

    monkeypatch.setitem(sys.modules, "scanpy", SimpleNamespace(read_h5ad=lambda _path: FakeAnnData()))
    monkeypatch.setitem(sys.modules, "brickembed", SimpleNamespace(combine_graph=fake_combine_graph))

    input_path = tmp_path / "strict.csv"
    genes_path = tmp_path / "genes.txt"
    graph_path = tmp_path / "graph.graphml"
    identity_path = tmp_path / "nodes.csv"
    audit_path = tmp_path / "audit.json"
    enrichment.to_csv(input_path, index=False)
    genes_path.write_text("GENE1\n", encoding="utf-8")
    args = SimpleNamespace(
        input=input_path,
        genes=genes_path,
        adata=tmp_path / "input.h5ad",
        output=graph_path,
        identity_output=identity_path,
        audit_output=audit_path,
        pvalue_cutoff=0.01,
        disease="COVID-19",
    )

    cmd_build_graph(args)

    graph_labels = set(captured["frame"]["path.2.name"])
    assert "Zinc [DRUGBANK:DB01593]" in graph_labels
    assert "Zinc [MESH:D015032]" in graph_labels
    assert "Zinc" not in graph_labels
    identity = pd.read_csv(identity_path).set_index("stable_id")
    assert identity.loc["DRUGBANK:DB01593", "display_label"] == "Zinc"
    assert identity.loc["MESH:D015032", "display_label"] == "Zinc"
    assert identity.loc["DRUGBANK:DB01593", "stable_id_count"] == 1
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    qualification = audit["stable_id_label_qualification"]
    assert qualification["qualified_entity_count"] == 2
    assert qualification["qualified_input_row_count"] == 2
    assert qualification["selection_or_filtering_used"] is False


def test_chemical_scope_is_metadata_based_and_non_filtering():
    scope = chemical_scope_table(_query_rows())
    by_id = scope.set_index("path.2.id")
    assert by_id.loc["CHEM:1", "chemical_scope"] == "curated_pharmacologic_entity"
    assert by_id.loc["CHEM:2", "chemical_scope"] == "chemical_entity_unspecified"
    assert scope["formal_candidate_eligible"].all()
    assert scope["scope_is_diagnostic_only"].all()


def test_extract_preserves_stable_ids_and_scope(tmp_path):
    ranked = rank(_query_rows(), metrics=["match_count"])
    ranked_path = tmp_path / "ranked.csv"
    scope_path = tmp_path / "scope.csv"
    output = tmp_path / "primary.csv"
    ranked.to_csv(ranked_path, index=False)
    chemical_scope_table(_query_rows()).to_csv(scope_path, index=False)
    _run(
        SCRIPTS / "drugdiscovery.py",
        "extract",
        "-i",
        ranked_path,
        "--scope",
        scope_path,
        "--top-n",
        "0",
        "-o",
        output,
    )
    result = pd.read_csv(output)
    assert list(result["candidate_rank"]) == [1, 2]
    assert set(result["path.2.id"]) == {"CHEM:1", "CHEM:2"}
    assert "chemical_scope" in result


def test_validation_uses_stable_ids_not_matching_labels(tmp_path):
    ranking = pd.DataFrame(
        [
            {"path.2.id": "CHEM:1", "path.2.name": "same", "path.2.type": "Chemical", "consensus_rank": 1},
            {"path.2.id": "CHEM:2", "path.2.name": "real hit", "path.2.type": "Chemical", "consensus_rank": 2},
        ]
    )
    reference = pd.DataFrame(
        [
            {"path.2.id": "CHEM:X", "path.2.name": "same", "path.2.type": "Chemical"},
            {"path.2.id": "CHEM:2", "path.2.name": "different label", "path.2.type": "Chemical"},
        ]
    )
    ranking_path = tmp_path / "ranking.csv"
    reference_path = tmp_path / "reference.csv"
    freeze_path = tmp_path / "ranking_frozen.sha256"
    output = tmp_path / "validation.csv"
    audit = tmp_path / "identity.csv"
    ranking.to_csv(ranking_path, index=False)
    freeze_ranking(ranking_path, freeze_path)
    time.sleep(0.01)
    reference.to_csv(reference_path, index=False)
    _run(
        SCRIPTS / "drugdiscovery.py",
        "validate",
        "--ranking",
        ranking_path,
        "--freeze",
        freeze_path,
        "--reference",
        reference_path,
        "--cutoffs",
        "1,2",
        "--identity-audit",
        audit,
        "-o",
        output,
    )
    metrics = pd.read_csv(output).set_index("cutoff")
    assert metrics.loc[1, "hit_count"] == 0
    assert metrics.loc[2, "hit_count"] == 1
    identity = pd.read_csv(audit)
    assert bool(identity.loc[0, "exact_label_hit_diagnostic"])
    assert not bool(identity.loc[0, "stable_id_hit"])
    assert bool(identity.loc[1, "stable_id_hit"])


def test_workflow_rank_is_separate_complete_and_post_freeze(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("TP53\nEGFR\n", encoding="utf-8")
    output = tmp_path / "workflow"
    _run(
        SCRIPTS / "workflow.py",
        "drug",
        "--genes",
        genes,
        "--method",
        "rank",
        "--disease",
        "COVID-19",
        "--evaluation-only",
        "--agent-plan",
        "-o",
        output,
    )
    commands = (output / "commands.sh").read_text(encoding="utf-8")
    assert f'brick rk -i {output / "query_chemical.csv"}' in commands
    assert f'brick rk -i {output / "query_disease.csv"}' in commands
    assert f'brick rk -i {output / "combined_query.csv"}' not in commands
    assert "brick drug scope" in commands
    assert "brick drug extract" in commands
    assert "brick drug freeze" in commands
    assert "brick drug validate" in commands
    assert commands.index("brick drug freeze") < commands.index("brick qr COVID-19")
    assert commands.index("brick qr COVID-19") < commands.index("brick drug validate")
    agent_plan = json.loads((output / "agent_plan.json").read_text(encoding="utf-8"))
    merge_step = next(step for step in agent_plan["steps"] if step["command"].startswith("brick plan merge-csv"))
    assert merge_step["type"] == "utility"
    assert validate_agent_plan(agent_plan, base_dir=output, workdir=output) == []


def test_workflow_blocks_reference_without_evaluation_only(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("TP53\n", encoding="utf-8")
    result = _run(
        SCRIPTS / "workflow.py",
        "drug",
        "--genes",
        genes,
        "--disease",
        "COVID-19",
        "-o",
        tmp_path / "workflow",
        check=False,
    )
    assert result.returncode != 0
    assert "--evaluation-only" in result.stderr


def test_workflow_grace_has_exact_background_cuda_and_freeze_before_validation(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("TP53\n", encoding="utf-8")
    output = tmp_path / "grace"
    _run(
        SCRIPTS / "workflow.py",
        "drug",
        "--genes",
        genes,
        "--adata",
        tmp_path / "input.h5ad",
        "--method",
        "grace",
        "--disease",
        "COVID-19",
        "--evaluation-only",
        "--agent-plan",
        "-o",
        output,
    )
    commands = (output / "commands.sh").read_text(encoding="utf-8")
    assert "--background-method auto" in commands
    assert "--pvalue-cutoff 0.0001" in commands
    assert "query_subset_degraded" not in commands
    assert "brick-gpu drug grace-rank" in commands
    assert "--device cuda" in commands
    assert commands.index("brick-gpu drug grace-rank") < commands.index("brick qr COVID-19")
    agent_plan = json.loads((output / "agent_plan.json").read_text(encoding="utf-8"))
    assert validate_agent_plan(agent_plan, base_dir=output, workdir=output) == []


def test_brickplan_drug_keeps_validation_after_freeze(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("TP53\n", encoding="utf-8")
    output = tmp_path / "plan.json"
    _run(
        SCRIPTS / "brickplan.py",
        "template",
        "--task",
        "drug",
        "--genes",
        genes,
        "--disease",
        "COVID-19",
        "--evaluation-only",
        "-o",
        output,
    )
    text = output.read_text(encoding="utf-8")
    assert text.index('"id": "freeze_primary_ranking"') < text.index('"id": "query_known_drugs_post_freeze"')
    assert text.index('"id": "query_known_drugs_post_freeze"') < text.index('"id": "validate_frozen_ranking"')
