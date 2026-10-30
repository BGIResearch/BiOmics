import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from agent_checker import has_errors, has_failed, run_artifact_validators, validate_agent_plan


def _subprocess_env():
    env = os.environ.copy()
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = str(SCRIPTS) if not existing else str(SCRIPTS) + os.pathsep + existing
    env.setdefault("PYTHONIOENCODING", "utf-8")
    return env


def _agent_plan(tmp_path, command=None):
    workdir = tmp_path / "run"
    workdir.mkdir()
    return {
        "schema_version": "1.0",
        "task": "unit",
        "workdir": str(workdir),
        "execution": {"mode": "dry_run", "executor": "brick", "allow_shell": False},
        "steps": [
            {
                "id": "query_kg",
                "type": "kg_query",
                "command": command
                or [
                    "brick",
                    "qr",
                    "-i",
                    "genes.txt",
                    "--st",
                    "Gene",
                    "--tt",
                    "Chemical",
                    "--hops",
                    "2",
                    "-o",
                    str(workdir / "kg.csv"),
                ],
                "inputs": ["genes.txt"],
                "outputs": [str(workdir / "kg.csv")],
                "depends_on": [],
                "validators": [
                    "kg_policy_not_broad",
                    "max_hops_respected",
                    "kg_edge_metadata_columns_present",
                ],
                "retry_policy": {"max_retries": 1},
                "checkpoint": True,
                "allow_failure": False,
            }
        ],
    }


def test_agent_plan_accepts_scoped_multihop_query(tmp_path):
    issues = validate_agent_plan(_agent_plan(tmp_path), base_dir=tmp_path)
    assert not has_errors(issues)


def test_agent_plan_rejects_broad_kg_query(tmp_path):
    plan = _agent_plan(tmp_path, command=["brick", "qr", "-i", "genes.txt", "--st", "Gene", "-o", str(tmp_path / "run" / "kg.csv")])
    issues = validate_agent_plan(plan, base_dir=tmp_path)
    assert has_errors(issues)
    assert any(issue.validator == "kg_policy_not_broad" for issue in issues)


def test_co_mentioned_validators_use_edge_metadata_and_weak_cap(tmp_path):
    path = tmp_path / "evidence.csv"
    pd.DataFrame(
        [
            {
                "path.0.name": "GeneA",
                "path.1.relation": "co_mentioned",
                "path.1.info_source": "PubMed",
                "path.1.original_relation": "GeneA and DiseaseX appeared in the same abstract.",
                "path.1.condition": "disease context",
                "path.1.relation_confidence": 0.8,
                "path.2.name": "DiseaseX",
                "evidence_strength": "direct",
                "manual_review_required": False,
            }
        ]
    ).to_csv(path, index=False)

    metadata = run_artifact_validators(path, ["kg_edge_metadata_required_for_co_mentioned"])
    assert not has_failed(metadata)

    weak_cap = run_artifact_validators(path, ["co_mentioned_weak_cap"])
    assert has_failed(weak_cap)


def test_co_mentioned_requires_original_relation_and_info_source(tmp_path):
    path = tmp_path / "kg.csv"
    pd.DataFrame(
        [
            {
                "path.1.relation": "co_mentioned",
                "path.1.info_source": "",
                "path.1.original_relation": "",
                "path.1.condition": "",
                "path.1.relation_confidence": "",
            }
        ]
    ).to_csv(path, index=False)
    results = run_artifact_validators(path, ["kg_edge_metadata_required_for_co_mentioned"])
    assert has_failed(results)


def test_brickplan_can_emit_agent_plan_and_agent_check_plan_accepts_it(tmp_path):
    genes = tmp_path / "genes.txt"
    genes.write_text("TP53\nEGFR\n", encoding="utf-8")
    out_dir = tmp_path / "model_run"
    plan_path = out_dir / "agent_plan.json"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "brickplan.py"),
            "template",
            "--task",
            "query-rank",
            "-i",
            str(genes),
            "--source-type",
            "Gene",
            "--target-type",
            "Chemical",
            "-r",
            "pharmaco_effect_on",
            "--agent-plan",
            "--output-dir",
            str(out_dir),
            "-o",
            str(plan_path),
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    assert "steps" in plan
    assert plan["steps"][0]["type"] == "kg_query"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "agent.py"),
            "check-plan",
            "-p",
            str(plan_path),
            "--format",
            "json",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )
    assert json.loads(result.stdout)["valid"] is True
