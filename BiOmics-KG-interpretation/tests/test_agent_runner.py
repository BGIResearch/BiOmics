import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from agent_runner import run_plan


def _write_plan(tmp_path, command, output_name="out.csv"):
    workdir = tmp_path / "run"
    workdir.mkdir()
    plan = {
        "schema_version": "1.0",
        "task": "runner_unit",
        "workdir": str(workdir),
        "execution": {
            "mode": "dry_run",
            "executor": "brick",
            "allow_shell": False,
            "max_retries": 0,
        },
        "steps": [
            {
                "id": "query_kg",
                "type": "kg_query",
                "command": command,
                "inputs": [],
                "outputs": [str(workdir / output_name)],
                "depends_on": [],
                "validators": ["kg_policy_not_broad", "max_hops_respected", "file_exists"],
                "retry_policy": {"max_retries": 0},
                "checkpoint": True,
                "allow_failure": False,
            }
        ],
    }
    path = tmp_path / "agent_plan.json"
    path.write_text(json.dumps(plan), encoding="utf-8")
    return path, workdir


def test_agent_runner_dry_run_writes_manifest_without_executing(tmp_path):
    output_file = tmp_path / "run" / "out.csv"
    plan_path, _workdir = _write_plan(
        tmp_path,
        [
            "brick",
            "qr",
            "TP53",
            "--st",
            "Gene",
            "--tt",
            "Chemical",
            "-o",
            str(output_file),
        ],
    )
    summary = run_plan(plan_path, execute=False)
    assert summary["ok"] is True
    manifest = Path(summary["manifest"])
    assert manifest.exists()
    events = [json.loads(line) for line in manifest.read_text(encoding="utf-8").splitlines()]
    assert any(event["event"] == "step_dry_run" for event in events)
    assert not output_file.exists()


def test_agent_runner_records_plan_failure_for_unsafe_command(tmp_path):
    plan_path, _workdir = _write_plan(tmp_path, ["python", "-c", "print(1)"])
    summary = run_plan(plan_path, execute=False)
    assert summary["ok"] is False
    events = [json.loads(line) for line in Path(summary["manifest"]).read_text(encoding="utf-8").splitlines()]
    assert any(event["event"] == "plan_issue" and "only allows commands" in event["message"] for event in events)


def test_agent_runner_dry_run_honors_step_dependencies(tmp_path):
    workdir = tmp_path / "run"
    workdir.mkdir()
    plan = {
        "schema_version": "1.0",
        "task": "ordered",
        "workdir": str(workdir),
        "execution": {"mode": "dry_run", "executor": "brick", "allow_shell": False},
        "steps": [
            {
                "id": "query_kg",
                "type": "kg_query",
                "command": ["brick", "qr", "TP53", "--st", "Gene", "--tt", "Chemical", "-o", str(workdir / "kg.csv")],
                "inputs": [],
                "outputs": [str(workdir / "kg.csv")],
                "depends_on": [],
                "validators": ["kg_policy_not_broad", "max_hops_respected"],
                "retry_policy": {"max_retries": 0},
                "checkpoint": True,
                "allow_failure": False,
            },
            {
                "id": "rank",
                "type": "ranking",
                "command": ["brick", "rk", "-i", str(workdir / "kg.csv"), "--rank", "-o", str(workdir / "ranked.csv")],
                "inputs": [str(workdir / "kg.csv")],
                "outputs": [str(workdir / "ranked.csv")],
                "depends_on": ["query_kg"],
                "validators": ["file_exists"],
                "retry_policy": {"max_retries": 0},
                "checkpoint": False,
                "allow_failure": False,
            },
        ],
    }
    plan_path = tmp_path / "agent_plan.json"
    plan_path.write_text(json.dumps(plan), encoding="utf-8")
    summary = run_plan(plan_path, execute=False)
    assert summary["ok"] is True
    events = [json.loads(line) for line in Path(summary["manifest"]).read_text(encoding="utf-8").splitlines()]
    started = [event["step_id"] for event in events if event["event"] == "step_started"]
    assert started == ["query_kg", "rank"]
