import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from kgpolicy import (
    DEFAULT_CATALOG,
    KNOWN_KG_RELATIONS,
    STRONG_CURATED_RELATIONS,
    WEAK_TRACE_RELATIONS,
    resolve_relation_aliases,
    select_kg_policy,
)


def _subprocess_env():
    env = os.environ.copy()
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = str(SCRIPTS) if not existing else str(SCRIPTS) + os.pathsep + existing
    env.setdefault("PYTHONIOENCODING", "utf-8")
    return env


def test_drug_policy_routes_to_curated_sources():
    policy = select_kg_policy(
        question="Which drugs target upregulated macrophage genes?",
        source_type="Gene",
        target_type="Chemical",
        kg_policy="auto",
    )
    assert policy["use_kg"] is True
    assert policy["intent"] == "drug_target_or_repositioning"
    assert "DrugBank" in policy["primary_info_sources"]
    assert "pharmaco_effect_on" in policy["relations"]
    assert policy["hops"] == 1


def test_auto_policy_does_not_guess_for_vague_gene_query():
    policy = select_kg_policy(
        question="TP53 neighbors",
        source_type="Gene",
        kg_policy="auto",
    )
    assert policy["intent"] == "ambiguous"
    assert policy["kg_policy"] == "manual_required"
    assert policy["info_sources"] == []
    assert policy["broad_query"] is True


def test_relation_aliases_resolve_to_known_kg_relations():
    resolved, unknown, aliases = resolve_relation_aliases("marker_of,regulate,associated_with")
    assert "produce" in resolved
    assert "regulated" in resolved
    assert "associated_with" not in resolved
    assert "marker_of" not in resolved
    assert "regulate" not in resolved
    assert aliases["marker_of"] == ["produce"]
    assert aliases["associated_with"] == []
    assert unknown == []


def test_associated_with_is_dropped():
    resolved, unknown, aliases = resolve_relation_aliases("associated_with")
    assert resolved == []
    assert aliases["associated_with"] == []
    assert unknown == []


def test_no_relation_is_blocked():
    resolved, unknown, aliases = resolve_relation_aliases("no_relation")
    assert resolved == []
    assert aliases == {}
    assert unknown == []


def test_known_relations_include_supplementary_table_schema_terms():
    assert "chemical_relation" in KNOWN_KG_RELATIONS
    resolved, unknown, aliases = resolve_relation_aliases("chemical_relation")
    assert resolved == ["chemical_relation"]
    assert aliases == {}
    assert unknown == []


def test_schema_declares_evidence_boundary_relation_sets():
    assert "pharmaco_effect_on" in STRONG_CURATED_RELATIONS
    assert "causally_related_to" in STRONG_CURATED_RELATIONS
    assert "co_mentioned" in WEAK_TRACE_RELATIONS
    assert "mentioned" in WEAK_TRACE_RELATIONS
    assert "no_relation" in WEAK_TRACE_RELATIONS
    assert "no_relation" not in STRONG_CURATED_RELATIONS


def test_default_policy_catalog_matches_reference_json():
    reference = json.loads((ROOT / "references" / "kg_source_policy.json").read_text(encoding="utf-8"))
    assert DEFAULT_CATALOG == reference


def test_policy_export_default_matches_reference_json(tmp_path):
    output = tmp_path / "kg_source_policy.json"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "kgpolicy.py"),
            "export-default",
            "-o",
            str(output),
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    exported = json.loads(output.read_text(encoding="utf-8"))
    reference = json.loads((ROOT / "references" / "kg_source_policy.json").read_text(encoding="utf-8"))
    assert exported == reference


def test_querygraph_dry_run_uses_policy_without_api(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("GeneA\nGeneB\n", encoding="utf-8")
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "querygraph.py"),
            "-i",
            str(genes),
            "--st",
            "Gene",
            "--tt",
            "Chemical",
            "--kg-policy",
            "auto",
            "--dry-run-plan",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    plan = json.loads(result.stdout)
    assert "--info-source" in plan["command"]
    assert "DrugBank" in plan["info_sources"]
    assert "pharmaco_effect_on" in plan["relations"]


def test_querygraph_question_defaults_to_auto_policy():
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "querygraph.py"),
            "TP53",
            "--st",
            "Gene",
            "--tt",
            "Chemical",
            "--question",
            "which drugs target TP53",
            "--dry-run-plan",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    plan = json.loads(result.stdout)
    assert "DrugBank" in plan["info_sources"]
    assert plan["policy"]["kg_policy"] == "auto"


def test_querygraph_manual_relation_alias_is_resolved():
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "querygraph.py"),
            "GeneA",
            "--st",
            "Gene",
            "--tt",
            "Cell",
            "-r",
            "marker_of",
            "--dry-run-plan",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    plan = json.loads(result.stdout)
    assert plan["relations"] == ["produce"]


def test_querygraph_dry_run_drops_associated_with_and_softens_natural_language_target():
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "querygraph.py"),
            "colorectal cancer",
            "--st",
            "Disease",
            "--target",
            "consensus molecular subtype",
            "--tt",
            "Disease",
            "-r",
            "associated_with",
            "--dry-run-plan",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    plan = json.loads(result.stdout)
    assert plan["target_query_mode"] == "auto_soft"
    assert plan["target_entities_for_query_preview"] == []
    assert plan["relations"] in (None, [])
    assert plan["policy"]["relation_aliases"]["associated_with"] == []
    assert plan["recommended_action"] == "query_kg"
    assert plan["fallback_enabled"] is True
    assert [stage["stage"] for stage in plan["fallback_stages"]] == ["initial"]


def test_querygraph_dry_run_shows_fallback_stage_plan_for_hard_target():
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "querygraph.py"),
            "DrugA",
            "--st",
            "Chemical",
            "--target",
            "DiseaseX",
            "--target-mode",
            "strict",
            "--tt",
            "Disease",
            "-r",
            "pharmaco_effect_on",
            "--dry-run-plan",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    plan = json.loads(result.stdout)
    assert [stage["stage"] for stage in plan["fallback_stages"]] == [
        "initial",
        "fallback_no_relation",
        "fallback_relax_target_type",
        "fallback_drop_target",
        "fallback_drop_target_no_relation",
    ]


def test_querygraph_api_payload_uses_live_info_source_key(monkeypatch):
    import querygraph

    captured = {}

    class DummyResponse:
        def raise_for_status(self):
            return None

        def json(self):
            return {"success": True, "data": {"rows": [], "columns": ["path.1.info_source"]}}

    def fake_post_api(url, payload):
        captured.update(payload)
        return DummyResponse()

    monkeypatch.setattr(querygraph, "_post_api", fake_post_api)
    querygraph._query_relation_single_batch(
        ["TP53"],
        relation=["pharmaco_effect_on"],
        source_entity_type="Gene",
        target_entity_type="Chemical",
        relation_info_source=["DrugBank"],
        limit_number=1,
    )
    assert captured["info_source"] == ["DrugBank"]
    assert "relation_info_source" not in captured


def test_query_relation_drops_associated_with_before_api(monkeypatch):
    import querygraph

    captured = {}

    def fake_single_batch(**kwargs):
        captured.update(kwargs)
        return querygraph.pd.DataFrame()

    monkeypatch.setattr(querygraph, "_query_relation_single_batch", fake_single_batch)
    querygraph.query_relation(
        source_entity_set=["colorectal cancer"],
        target_entity_set=["consensus molecular subtype"],
        relation="associated_with",
        source_entity_type="Disease",
        target_entity_type="Disease",
        limit_number=1,
    )
    assert captured["relation"] is None


def test_query_relation_with_fallback_selects_first_nonempty_stage(monkeypatch):
    import querygraph

    calls = []

    def fake_query_relation(**kwargs):
        calls.append(kwargs)
        if kwargs["target_entity_set"] is None and kwargs["relation"] == ["pharmaco_effect_on"]:
            return querygraph.pd.DataFrame([
                {
                    "path.0.name": "DrugA",
                    "path.0.type": "Chemical",
                    "path.1.relation": "pharmaco_effect_on",
                    "path.2.name": "DiseaseX",
                    "path.2.type": "Disease",
                }
            ])
        return querygraph.pd.DataFrame()

    monkeypatch.setattr(querygraph, "query_relation", fake_query_relation)
    df, diagnostics = querygraph.query_relation_with_fallback(
        source_entity_set=["DrugA"],
        target_entity_set=["DiseaseX"],
        relation="pharmaco_effect_on",
        source_entity_type="Chemical",
        target_entity_type="Disease",
        target_mode="strict",
        return_diagnostics=True,
    )
    assert not df.empty
    assert diagnostics["selected_stage"] == "fallback_drop_target"
    assert diagnostics["fallback_drop_target_nonempty"] is True
    assert [call["target_entity_set"] for call in calls] == [
        ["DiseaseX"],
        ["DiseaseX"],
        ["DiseaseX"],
        None,
    ]


def test_query_relation_with_fallback_can_wait_for_candidate_eligible(monkeypatch):
    import querygraph

    calls = []

    def fake_query_relation(**kwargs):
        calls.append(kwargs)
        if (
            kwargs["target_entity_set"] == ["DiseaseX"]
            and kwargs["relation"] == ["pharmaco_effect_on"]
            and kwargs["target_entity_type"] == "Disease"
        ):
            return querygraph.pd.DataFrame([
                {
                    "path.0.name": "DrugA",
                    "path.0.type": "Chemical",
                    "path.1.relation": "co_mentioned",
                    "path.2.name": "DiseaseX",
                    "path.2.type": "Disease",
                }
            ])
        if kwargs["target_entity_set"] is None and kwargs["relation"] == ["pharmaco_effect_on"]:
            return querygraph.pd.DataFrame([
                {
                    "path.0.name": "DrugA",
                    "path.0.type": "Chemical",
                    "path.1.relation": "pharmaco_effect_on",
                    "path.2.name": "DiseaseX",
                    "path.2.type": "Disease",
                }
            ])
        return querygraph.pd.DataFrame()

    def candidate_evaluator(df):
        return "pharmaco_effect_on" in set(df["path.1.relation"].astype(str))

    monkeypatch.setattr(querygraph, "query_relation", fake_query_relation)
    df, diagnostics = querygraph.query_relation_with_fallback(
        source_entity_set=["DrugA"],
        target_entity_set=["DiseaseX"],
        relation="pharmaco_effect_on",
        source_entity_type="Chemical",
        target_entity_type="Disease",
        target_mode="strict",
        stop_condition="candidate_eligible",
        candidate_evaluator=candidate_evaluator,
        return_diagnostics=True,
    )
    assert len(df) == 2
    assert diagnostics["stopped_by"] == "candidate_eligible"
    assert diagnostics["selected_stage"] == "fallback_drop_target"
    assert diagnostics["selected_nonempty_stages"] == ["initial", "fallback_drop_target"]
    assert [call["target_entity_set"] for call in calls] == [
        ["DiseaseX"],
        ["DiseaseX"],
        ["DiseaseX"],
        None,
    ]


def test_querygraph_blocks_unscoped_broad_query(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("GeneA\n", encoding="utf-8")
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "querygraph.py"),
            "-i",
            str(genes),
            "--st",
            "Gene",
            "--dry-run-plan",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        capture_output=True,
    )
    assert result.returncode != 0
    assert "Broad KG neighbor query" in result.stderr


def test_querygraph_dry_run_strips_utf8_bom(tmp_path):
    genes = tmp_path / "genes_bom.txt"
    genes.write_text("GeneA\nGeneB\n", encoding="utf-8-sig")
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "querygraph.py"),
            "-i",
            str(genes),
            "--st",
            "Gene",
            "--tt",
            "Chemical",
            "--kg-policy",
            "auto",
            "--dry-run-plan",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    plan = json.loads(result.stdout)
    assert plan["source_entities_preview"][0] == "GeneA"


def test_drug_workflow_commands_merge_before_reasoning(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("TP53\nEGFR\n", encoding="utf-8")
    out_dir = tmp_path / "wf_drug"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "workflow.py"),
            "drug",
            "--genes",
            str(genes),
            "--kg-policy",
            "auto",
            "--question",
            "find drugs for these genes",
            "-o",
            str(out_dir),
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    commands = (out_dir / "commands.ps1").read_text(encoding="utf-8")
    assert "brick plan merge-csv" in commands
    assert commands.index("brick plan merge-csv") < commands.index("brick reason normalize")


def test_celltype_workflow_uses_produce_relation(tmp_path):
    marker_dir = tmp_path / "markers"
    marker_dir.mkdir()
    (marker_dir / "cluster_0.txt").write_text("GeneA\nGeneB\n", encoding="utf-8")
    out_dir = tmp_path / "wf_celltype"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "workflow.py"),
            "celltype",
            "--marker-dir",
            str(marker_dir),
            "--kg-policy",
            "auto",
            "-o",
            str(out_dir),
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    commands = (out_dir / "commands.ps1").read_text(encoding="utf-8")
    assert "-r produce" in commands
    assert "marker_of" not in commands


def test_enrichment_workflow_lets_policy_choose_relations(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("TP53\nEGFR\n", encoding="utf-8")
    out_dir = tmp_path / "wf_enrich"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "workflow.py"),
            "enrichment",
            "--genes",
            str(genes),
            "--target-types",
            "Chemical,Disease",
            "--kg-policy",
            "auto",
            "-o",
            str(out_dir),
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    commands = (out_dir / "commands.ps1").read_text(encoding="utf-8")
    assert "brick plan merge-csv" in commands
    assert commands.index("brick plan merge-csv") < commands.index("brick reason normalize")
    assert "pharmaco_effect_on,interact_with,regulated" in commands
    assert "causally_related_to,disease_relation,participate_in" in commands
    assert "marker_of" not in commands
