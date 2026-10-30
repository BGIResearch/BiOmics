import math
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import rankgraph


def query_frame():
    return pd.DataFrame(
        {
            "path.0.id": ["S1", "S2", "S3"],
            "path.0.name": ["SOURCE1", "SOURCE2", "SOURCE3"],
            "path.0.type": ["SourceType"] * 3,
            "path.1.relation": ["related_to"] * 3,
            "path.2.id": ["T1", "T1", "T2"],
            "path.2.name": ["Target One", "Target One", "Target Two"],
            "path.2.type": ["TargetType"] * 3,
        }
    )


def reverse_rows(include_second=True):
    rows = [
        {
            "path.0.id": "T1",
            "path.0.name": "Target One",
            "path.1.relation": "related_to",
            "path.1.info_source": ["E1"],
            "path.2.id": "B1",
            "path.2.name": "Background One",
        },
        {
            "path.0.id": "T1",
            "path.0.name": "Target One",
            "path.1.relation": "related_to",
            "path.1.info_source": ["E2"],
            "path.2.id": "B2",
            "path.2.name": "Background Two",
        },
        {
            "path.0.id": "T1",
            "path.0.name": "Target One",
            "path.1.relation": "related_to",
            "path.1.info_source": ["E3"],
            "path.2.id": "B3",
            "path.2.name": "Background Three",
        },
    ]
    rows.append(dict(rows[0]))  # exact transport duplicate; must count once
    if include_second:
        rows.append(
            {
                "path.0.id": "T2",
                "path.0.name": "Target Two",
                "path.1.relation": "related_to",
                "path.1.info_source": ["E4"],
                "path.2.id": "B4",
                "path.2.name": "Background Four",
            }
        )
    return pd.DataFrame(rows)


class RelationBackgroundTests(unittest.TestCase):
    def test_query_subset_degraded_is_local_and_explicitly_labeled(self):
        with (
            patch.object(
                rankgraph,
                "query_cypher",
                side_effect=AssertionError("degraded mode must not call Cypher"),
            ),
            patch.object(
                rankgraph,
                "_query_relation_single_batch",
                side_effect=AssertionError("degraded mode must not call live relation API"),
            ),
        ):
            result = rankgraph.enrich(
                query_frame(),
                background_method="query_subset_degraded",
            )

        by_id = result.set_index("path.2.id")
        self.assertEqual(int(by_id.loc["T1", "path.2.background_count"]), 2)
        self.assertEqual(int(by_id.loc["T2", "path.2.background_count"]), 1)
        self.assertTrue(
            result["path.2.background_method"]
            .eq("query_subset_degraded_not_kg_wide")
            .all()
        )
        self.assertTrue(
            result["path.2.enrichment_status"]
            .eq("degraded_query_subset_not_kg_wide")
            .all()
        )
        self.assertTrue(
            result["path.2.enrich_fdr_status"]
            .eq("degraded_query_subset_all_targets")
            .all()
        )

    def test_kg_wide_targets_are_queried_individually(self):
        calls = []

        def fake_relation(**kwargs):
            calls.append(kwargs)
            self.assertEqual(len(kwargs["source_entity_set"]), 1)
            target_id = kwargs["source_entity_set"][0]
            target_name = "Target One" if target_id == "T1" else "Target Two"
            background_count = 2 if target_id == "T1" else 1
            return pd.DataFrame(
                [
                    {
                        "path.0.id": target_id,
                        "path.0.name": target_name,
                        "path.1.relation": "related_to",
                        "path.2.id": f"B-{target_id}-{index}",
                        "path.2.name": f"Background {index}",
                    }
                    for index in range(background_count)
                ]
            )

        with patch.object(rankgraph, "_query_relation_single_batch", side_effect=fake_relation):
            result = rankgraph.enrich(
                query_frame(),
                background_method="relation",
                relation_batch_size=10,
                relation_workers=1,
                relation_retries=0,
            )

        self.assertEqual(len(calls), 2)
        self.assertEqual([call["source_entity_set"] for call in calls], [["T1"], ["T2"]])
        self.assertTrue(all(call["limit_number"] == 10001 for call in calls))
        self.assertTrue(
            result["path.2.background_method"]
            .str.startswith("relation_id_target_relations_kg_wide")
            .all()
        )

    def test_each_target_uses_only_its_observed_relation_signature(self):
        frame = pd.DataFrame(
            {
                "path.0.id": ["S1", "S2"],
                "path.0.name": ["SOURCE1", "SOURCE2"],
                "path.0.type": ["SourceType", "SourceType"],
                "path.1.relation": ["relation_a", "relation_b"],
                "path.2.id": ["T1", "T2"],
                "path.2.name": ["Target One", "Target Two"],
                "path.2.type": ["TargetType", "TargetType"],
            }
        )
        calls = []

        def fake_relation(**kwargs):
            calls.append(kwargs)
            target_id = kwargs["source_entity_set"][0]
            target_name = "Target One" if target_id == "T1" else "Target Two"
            background_id = "B1" if target_id == "T1" else "B2"
            return pd.DataFrame(
                [
                    {
                        "path.0.id": target_id,
                        "path.0.name": target_name,
                        "path.1.relation": kwargs["relation"][0],
                        "path.2.id": background_id,
                        "path.2.name": "Background",
                    }
                ]
            )

        with patch.object(rankgraph, "_query_relation_single_batch", side_effect=fake_relation):
            result = rankgraph.enrich(
                frame,
                background_method="relation",
                relation_batch_size=10,
                relation_workers=1,
            )

        self.assertEqual(len(calls), 2)
        by_target = {call["source_entity_set"][0]: call["relation"] for call in calls}
        self.assertEqual(by_target, {"T1": ("relation_a",), "T2": ("relation_b",)})
        self.assertTrue(
            result["path.2.background_method"]
            .str.startswith("relation_id_target_relations")
            .all()
        )

    def test_transient_timeout_is_retried_and_audited(self):
        timeout = requests.Timeout("temporary KG timeout")

        calls = []

        def fake_relation(**kwargs):
            calls.append(kwargs)
            if len(calls) == 1:
                raise timeout
            target_id = kwargs["source_entity_set"][0]
            background_count = 2 if target_id == "T1" else 1
            return pd.DataFrame(
                [
                    {
                        "path.0.id": target_id,
                        "path.0.name": "Target One" if target_id == "T1" else "Target Two",
                        "path.1.relation": "related_to",
                        "path.2.id": f"B-{target_id}-{index}",
                        "path.2.name": f"Background {index}",
                    }
                    for index in range(background_count)
                ]
            )

        with (
            patch.object(
                rankgraph,
                "_query_relation_single_batch",
                side_effect=fake_relation,
            ) as relation_call,
            patch.object(rankgraph.time, "sleep"),
        ):
            result = rankgraph.enrich(
                query_frame(),
                background_method="relation",
                relation_batch_size=10,
                relation_workers=1,
                relation_retries=1,
            )

        self.assertEqual(relation_call.call_count, 3)
        by_id = result.set_index("path.2.id")
        self.assertIn("transport_retry_1", by_id.loc["T1", "path.2.background_method"])
        self.assertNotIn("transport_retry_1", by_id.loc["T2", "path.2.background_method"])

    def test_nontransient_http_error_is_not_retried(self):
        response = requests.Response()
        response.status_code = 400
        invalid_request = requests.HTTPError(response=response)
        with patch.object(
            rankgraph,
            "_query_relation_single_batch",
            side_effect=invalid_request,
        ) as relation_call:
            with self.assertRaises(requests.HTTPError):
                rankgraph.enrich(
                    query_frame(),
                    background_method="relation",
                    relation_batch_size=10,
                    relation_workers=1,
                    relation_retries=2,
                )

        self.assertEqual(relation_call.call_count, 1)

    def test_same_label_targets_remain_separate_by_stable_id(self):
        frame = pd.DataFrame(
            {
                "path.0.id": ["S1", "S2"],
                "path.0.name": ["SOURCE1", "SOURCE2"],
                "path.0.type": ["SourceType", "SourceType"],
                "path.1.relation": ["related_to", "related_to"],
                "path.2.id": ["T1", "T2"],
                "path.2.name": ["Shared Label", "Shared Label"],
                "path.2.type": ["TargetType", "TargetType"],
            }
        )
        with patch.object(
            rankgraph,
            "_query_relation_single_batch",
            return_value=reverse_rows(include_second=True),
        ):
            result = rankgraph.enrich(
                frame,
                background_method="relation",
                relation_batch_size=10,
                relation_workers=1,
            )

        self.assertEqual(set(result["path.2.id"]), {"T1", "T2"})
        self.assertEqual(len(result), 2)

    def test_exact_bound_screen_skips_targets_that_cannot_pass(self):
        calls = []

        def fake_relation(**kwargs):
            calls.append(kwargs)
            return reverse_rows(include_second=False)

        with patch.object(rankgraph, "_query_relation_single_batch", side_effect=fake_relation):
            result = rankgraph.enrich(
                query_frame(),
                background_method="relation",
                pvalue_cutoff=1e-4,
                relation_batch_size=10,
                relation_workers=1,
            )

        by_id = result.set_index("path.2.id")
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["source_entity_set"], ["T1"])
        self.assertEqual(int(by_id.loc["T1", "path.2.background_count"]), 3)
        self.assertEqual(by_id.loc["T1", "path.2.enrichment_status"], "exact")
        self.assertTrue(math.isnan(by_id.loc["T2", "path.2.enrich_pvalue"]))
        self.assertEqual(
            by_id.loc["T2", "path.2.enrichment_status"],
            "screened_exact_bound",
        )
        self.assertEqual(
            by_id.loc["T2", "path.2.background_method"],
            "not_queried_exact_bound",
        )

    def test_kg_wide_failure_bound_has_no_pseudo_pvalue(self):
        frame = query_frame().loc[lambda value: value["path.2.id"] == "T1"]
        source_n = frame["path.0.id"].nunique()
        threshold = rankgraph._background_failure_threshold(2, source_n, 1e-4)
        calls = []

        def fake_relation(**kwargs):
            calls.append(kwargs)
            return pd.DataFrame(
                [
                    {
                        "path.0.id": "T1",
                        "path.0.name": "Target One",
                        "path.1.relation": "related_to",
                        "path.2.id": f"B{i}",
                        "path.2.name": f"Background {i}",
                    }
                    for i in range(threshold)
                ]
            )

        with patch.object(rankgraph, "_query_relation_single_batch", side_effect=fake_relation):
            result = rankgraph.enrich(
                frame,
                background_method="relation",
                pvalue_cutoff=1e-4,
                relation_workers=1,
            )

        row = result.iloc[0]
        self.assertEqual(calls[0]["limit_number"], threshold)
        self.assertEqual(row["path.2.enrichment_status"], "screened_kg_wide_failure_bound")
        self.assertTrue(math.isnan(row["path.2.enrich_pvalue"]))

    def test_server_cap_is_excluded_as_hub_not_treated_as_exact(self):
        frame = query_frame().loc[lambda value: value["path.2.id"] == "T1"]

        def fake_relation(**kwargs):
            return pd.DataFrame(
                [
                    {
                        "path.0.id": "T1",
                        "path.0.name": "Target One",
                        "path.1.relation": "related_to",
                        "path.2.id": f"B{i}",
                        "path.2.name": f"Background {i}",
                    }
                    for i in range(5)
                ]
            )

        with patch.object(rankgraph, "_query_relation_single_batch", side_effect=fake_relation):
            result = rankgraph.enrich(
                frame,
                background_method="relation",
                pvalue_cutoff=1e-4,
                relation_workers=1,
                relation_server_row_cap=5,
            )

        row = result.iloc[0]
        self.assertEqual(row["path.2.enrichment_status"], "excluded_kg_wide_saturated_hub")
        self.assertTrue(math.isnan(row["path.2.enrich_pvalue"]))
        self.assertIn("server_cap_5", row["path.2.background_method"])

    def test_auto_falls_back_only_for_missing_cypher_endpoint(self):
        response = requests.Response()
        response.status_code = 404
        missing_endpoint = requests.HTTPError(response=response)

        with (
            patch.object(rankgraph, "query_cypher", side_effect=missing_endpoint),
            patch.object(
                rankgraph,
                "_query_relation_single_batch",
                return_value=reverse_rows(include_second=True),
            ),
        ):
            result = rankgraph.enrich(
                query_frame(),
                background_method="auto",
                relation_batch_size=10,
                relation_workers=1,
            )

        by_id = result.set_index("path.2.id")
        self.assertEqual(int(by_id.loc["T1", "path.2.background_count"]), 3)
        self.assertEqual(int(by_id.loc["T2", "path.2.background_count"]), 1)
        self.assertTrue(by_id["path.2.enrich_pvalue"].notna().all())
        self.assertTrue(by_id["path.2.enrich_fdr"].notna().all())

    def test_incomplete_live_background_is_a_hard_failure(self):
        frame = query_frame().loc[lambda value: value["path.2.id"] == "T1"]
        with patch.object(
            rankgraph,
            "_query_relation_single_batch",
            return_value=pd.DataFrame(),
        ):
            with self.assertRaisesRegex(RuntimeError, "background is incomplete"):
                rankgraph.enrich(
                    frame,
                    background_method="relation",
                    relation_batch_size=10,
                    relation_workers=1,
                )


if __name__ == "__main__":
    unittest.main()
