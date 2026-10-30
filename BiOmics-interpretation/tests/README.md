# Skill Tests

This folder contains lightweight tests for the BiOmics-BRICK skill.

## Files

- `test_agent_checker.py`
- `test_agent_runner.py`
- `test_drug_discovery.py`
- `test_embedgraph.py`
- `test_evidence_packet.py`
- `test_grace.py`
- `test_kg_policy.py`
- `test_live_api_schema.py`
- `test_querygraph_entity_parsing.py`
- `test_rankgraph_enrichment_fallback.py`
- `test_reasoning_evidence.py`
- `test_reasoning_literature.py`
- `test_task_contract.py`

## Run

From the skill directory:

```bash
pip install pytest
python -m pytest tests/
```

The live API schema test is optional and should only be run when KG API access
is intentionally enabled:

```bash
RUN_BRICK_API_TESTS=1 python -m pytest tests/test_live_api_schema.py
```
