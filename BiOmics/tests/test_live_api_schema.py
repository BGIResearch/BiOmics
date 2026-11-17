import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


pytestmark = pytest.mark.skipif(
    os.environ.get("RUN_BRICK_API_TESTS") != "1",
    reason="set RUN_BRICK_API_TESTS=1 to run live BiOmics-KG API smoke tests",
)


def _run_query(output, *args):
    env = os.environ.copy()
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = str(SCRIPTS) if not existing else str(SCRIPTS) + os.pathsep + existing
    env.setdefault("PYTHONIOENCODING", "utf-8")
    subprocess.run(
        [sys.executable, str(SCRIPTS / "querygraph.py"), *args, "-n", "1", "-o", str(output)],
        cwd=str(ROOT),
        env=env,
        text=True,
        check=True,
        capture_output=True,
    )


def test_live_api_accepts_info_source_alias_and_hops(tmp_path):
    onehop = tmp_path / "tp53_drugbank.csv"
    _run_query(
        onehop,
        "TP53",
        "--st",
        "Gene",
        "--tt",
        "Chemical",
        "-r",
        "pharmaco_effect_on",
        "--info-source",
        "DrugBank",
    )
    if not onehop.exists() or onehop.stat().st_size == 0:
        pytest.skip("live API returned no one-hop rows for the compact smoke query")
    df = pd.read_csv(onehop)
    if df.empty:
        pytest.skip("live API returned an empty one-hop table for the compact smoke query")
    reviewable_fields = {"path.1.original_relation", "path.1.info_source", "path.1.condition", "path.1.relation_confidence"}
    assert reviewable_fields.intersection(df.columns)

    multihop = tmp_path / "tp53_disease_2hop.csv"
    _run_query(
        multihop,
        "TP53",
        "--st",
        "Gene",
        "--tt",
        "Disease",
        "--kg-policy",
        "auto",
        "--hops",
        "2",
    )
    assert multihop.exists()
