import json
import os
import shlex
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from querygraph import _extract_entities_from_text, _load_entities_from_file


def _subprocess_env():
    env = os.environ.copy()
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = str(SCRIPTS) if not existing else str(SCRIPTS) + os.pathsep + existing
    env.setdefault("PYTHONIOENCODING", "utf-8")
    return env


def test_inline_parser_preserves_internal_spaces_and_hgvs_colon():
    assert _extract_entities_from_text("rs7756992 A G") == ["rs7756992 A G"]
    assert _extract_entities_from_text("NC_000010.11:g.112988738=") == [
        "NC_000010.11:g.112988738="
    ]


def test_entity_files_preserve_spaces_commas_colons_and_parentheses(tmp_path):
    entities = tmp_path / "entities.txt"
    entities.write_text(
        "NC_000010.11:g.112988738=\n"
        "rs7756992 A G\n"
        "Diabetes Mellitus, Type 2\n"
        "NM_017774.3(CDKAL1):c.371+30101A>G\n",
        encoding="utf-8",
    )

    assert _load_entities_from_file(str(entities)) == [
        "NC_000010.11:g.112988738=",
        "rs7756992 A G",
        "Diabetes Mellitus, Type 2",
        "NM_017774.3(CDKAL1):c.371+30101A>G",
    ]


def test_dry_run_shell_command_round_trips_multiword_entity():
    source = "rs7756992 A G"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "querygraph.py"),
            source,
            "--st",
            "Mutation",
            "--tt",
            "Disease",
            "--dry-run-plan",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )

    plan = json.loads(result.stdout)
    assert plan["source_entities_preview"] == [source]
    assert shlex.split(plan["command_sh"]) == plan["command"]


def test_cli_rejoins_unquoted_snp_words_before_parsing():
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "querygraph.py"),
            "NC_000010.11:g.112998590",
            "A",
            "G",
            "--st",
            "Mutation",
            "--tt",
            "Disease",
            "--dry-run-plan",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )

    plan = json.loads(result.stdout)
    assert plan["source_entities_preview"] == ["NC_000010.11:g.112998590 A G"]
    assert plan["command_sh"].startswith(
        "brick qr 'NC_000010.11:g.112998590 A G'"
    )


def test_dry_run_files_keep_exact_snp_and_disease_literals(tmp_path):
    sources = tmp_path / "snps.txt"
    targets = tmp_path / "diseases.txt"
    sources.write_text("NC_000010.11:g.112988738=\n", encoding="utf-8")
    targets.write_text("Diabetes Mellitus, Type 2\n", encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "querygraph.py"),
            "-i",
            str(sources),
            "--target-file",
            str(targets),
            "--st",
            "Mutation",
            "--tt",
            "Disease",
            "--dry-run-plan",
        ],
        cwd=str(ROOT),
        env=_subprocess_env(),
        text=True,
        check=True,
        capture_output=True,
    )

    plan = json.loads(result.stdout)
    assert plan["source_entities_preview"] == ["NC_000010.11:g.112988738="]
    assert plan["target_entities_preview"] == ["Diabetes Mellitus, Type 2"]
