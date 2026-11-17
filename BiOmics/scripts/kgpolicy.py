"""Controlled KG source/hop policy for BiOmics-BRICK."""

from __future__ import annotations

import argparse
import ast
import json
import re
import shlex
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

MAX_HOPS = 3
KG_POLICY_CHOICES = ("manual", "auto", "off", "direct", "source-scoped", "multihop")
MIN_AUTO_INTENT_CONFIDENCE = 0.25
DEFAULT_SCHEMA_PATH = Path(__file__).resolve().parents[1] / "references" / "kg_schema.json"


def load_kg_schema(schema_file: Optional[str] = None) -> Dict[str, Any]:
    path = Path(schema_file) if schema_file else DEFAULT_SCHEMA_PATH
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {
        "queryable_edge_types": [],
        "blocked_query_edge_types": ["no_relation"],
        "weak_evidence_edge_types": ["mentioned", "co_mentioned"],
        "relation_aliases": {},
        "relation_compatibility": {},
    }


_SCHEMA = load_kg_schema()
QUERYABLE_KG_RELATIONS = set(_SCHEMA.get("queryable_edge_types", []))
WEAK_EVIDENCE_RELATIONS = set(_SCHEMA.get("weak_evidence_edge_types", []))
BLOCKED_QUERY_RELATIONS = set(_SCHEMA.get("blocked_query_edge_types", ["no_relation"]))
WEAK_TRACE_RELATIONS = set(
    _SCHEMA.get("weak_trace_edge_types", list(WEAK_EVIDENCE_RELATIONS | BLOCKED_QUERY_RELATIONS))
)
STRONG_CURATED_RELATIONS = set(_SCHEMA.get("strong_curated_edge_types", []))
KNOWN_KG_RELATIONS = QUERYABLE_KG_RELATIONS | WEAK_EVIDENCE_RELATIONS
RELATION_ALIASES = dict(_SCHEMA.get("relation_aliases", {}))

INFO_SOURCE_CANONICAL_ALIASES = {
    "clinical trials": "ClinicalTrials",
    "clinicaltrials": "ClinicalTrials",
    "clinicaltrials.gov": "ClinicalTrials",
    "clinicaltrails": "ClinicalTrials",
    "gwas_cataog": "GWAS_Catalog",
    "gwas catalog": "GWAS_Catalog",
    "gwas_catalog": "GWAS_Catalog",
    "gwas-catalog": "GWAS_Catalog",
    "gwascatalog": "GWAS_Catalog",
    "ebi gwas catalog": "GWAS_Catalog",
    "nhgri-ebi gwas catalog": "GWAS_Catalog",
}

INFO_SOURCE_QUERY_ALIASES = {
    "ClinicalTrials": [
        "ClinicalTrials",
        "ClinicalTrails",
        "ClinicalTrials.gov",
    ],
    "GWAS_Catalog": [
        "GWAS_Catalog",
        "GWAS Catalog",
        "GWAS_cataog",
        "gwas_catalog",
        "GWASCatalog",
        "EBI GWAS Catalog",
    ],
}

DEFAULT_CATALOG: Dict[str, Any] = {
    "version": 1,
    "default_edge_fields": ["relation", "original_relation", "info_source", "condition", "relation_confidence"],
    "co_mentioned_policy": {
        "claim_allowed": False,
        "max_evidence_strength": "weak",
        "requires_fields": ["original_relation", "info_source", "condition", "relation_confidence"],
    },
    "source_quality_tiers": {
        "curated_clinical": ["ClinVar", "ClinVar_gene", "OMIM", "MEDGEN", "ORPHANET", "DrugBank", "DrugCentral", "DRUGCENTRAL", "IUPHAR", "PHARMGKB", "ClinicalTrials", "ClinicalTrails", "REPODB"],
        "curated_genetic_association": ["GWAS_Catalog", "GWAS_cataog"],
        "curated_biological": ["GOA", "GOA_mm", "GENEONTOLOGY", "Reactome", "Reactome_anno", "KEGG", "KEGG_anno", "KEGG_anno_mm", "RHEA", "CORUM", "IntAct", "INTACT", "BioGRID", "BIOGRID", "BioGRID_hs", "BioGRID_mm", "SIGNOR", "UNIPROT", "UNIPROT_human", "UNIPROT_rodents"],
        "single_cell_curated": ["CellMarker", "CellMarker_mm", "PanglaoDB", "CellTaxonomy", "CellChat", "CellTalkDB", "CellTalkDB_mm", "SCENIC", "SCENIC_mm"],
        "ontology": ["CL", "CLO", "CELLLINEONTOLOGY", "UBERON", "GO", "GO-PLUS", "HPO", "HUMAN PHENOTYPE ONTOLOGY", "MONDO", "DOID", "BTO", "CHEBI", "EFO", "PATO", "NCIT", "SNOMED_CT", "MEDDRA", "MESH", "MPATH", "NBO", "UPHENO", "UPHENO_ALL_WITH_RELATIONS", "XREF", "MA", "EMAPA", "FOODON"],
        "integrated_network": ["PrimeKG", "Bioportal", "AGRKB", "AGRKB_Disease", "AGRKB_EXPRESSION", "AGRKB_GENE", "BioSNAP_ChCh-Miner", "BioSNAP_ChChSe-Decagon", "BioSNAP_ChG-Miner", "BioSNAP_ChSe-Decagon", "BioSNAP_DCh-Miner", "BioSNAP_DG-Miner"],
        "predicted_or_text_mined": ["StringDB", "STITCH", "StingDB", "StingDB_mm", "SIDER", "PubMed text-mined"],
    },
    "intents": [
        {
            "intent": "cell_marker_annotation",
            "keywords": ["cell type", "celltype", "marker", "annotation", "cluster", "single cell"],
            "source_type": "Gene", "target_type": "Cell",
            "relations": ["produce"],
            "primary_info_sources": ["CellMarker", "CellMarker_mm", "PanglaoDB", "CellTaxonomy"],
            "fallback_info_sources": ["CL", "CLO", "BTO", "EFO", "NCIT"],
            "default_hops": 1, "allow_multihop": False,
            "avoid_relations": ["co_mentioned", "no_relation"],
        },
        {
            "intent": "cell_cell_communication",
            "keywords": ["communication", "ligand", "receptor", "cell-cell", "cellchat", "celltalk"],
            "source_type": "Cell", "target_type": "Cell",
            "relations": ["interact_with", "regulated", "functionally_related_to"],
            "primary_info_sources": ["CellChat", "CellTalkDB", "CellTalkDB_mm"],
            "fallback_info_sources": ["BioGRID_hs", "IntAct", "StringDB"],
            "default_hops": 2, "allow_multihop": True,
            "avoid_relations": ["co_mentioned", "no_relation"],
        },
        {
            "intent": "drug_target_or_repositioning",
            "keywords": ["drug", "chemical", "compound", "therapy", "treatment", "inhibitor", "agonist", "antagonist", "repurposing", "repositioning"],
            "source_type": "Gene", "target_type": "Chemical",
            "relations": ["pharmaco_effect_on", "interact_with", "regulated"],
            "primary_info_sources": ["DrugBank", "DrugCentral", "DRUGCENTRAL", "IUPHAR", "PHARMGKB", "DrugCentral_hs", "DrugCentral_mm"],
            "fallback_info_sources": ["STITCH", "BIOSNAP_ChG-Miner", "BIOSNAP_DG-Miner", "PrimeKG"],
            "default_hops": 1, "allow_multihop": True,
            "avoid_relations": ["co_mentioned", "no_relation"],
        },
        {
            "intent": "disease_gene_or_variant",
            "keywords": ["disease", "phenotype", "variant", "mutation", "causal", "pathogenic", "symptom", "gwas"],
            "source_type": "Gene", "target_type": "Disease",
            "relations": ["causally_related_to", "disease_relation", "participate_in", "is_a"],
            "primary_info_sources": ["ClinVar", "ClinVar_gene", "GWAS_Catalog", "OMIM", "HPO_g2p_anno", "HPO_anno", "MONDO", "DOID", "MEDGEN", "ORPHANET"],
            "fallback_info_sources": ["PrimeKG", "AGRKB_Disease", "PHARMGKB"],
            "default_hops": 1, "allow_multihop": True,
            "avoid_relations": ["co_mentioned", "no_relation"],
        },
        {
            "intent": "pathway_function_enrichment",
            "keywords": ["pathway", "enrichment", "function", "process", "reaction", "go", "reactome", "kegg"],
            "source_type": "Gene", "target_type": "Pathway",
            "relations": ["participate_in", "functionally_related_to", "regulated"],
            "primary_info_sources": ["GOA", "GOA_mm", "Reactome_anno", "Reactome", "KEGG_anno", "KEGG_anno_mm", "KEGG", "RHEA", "SMPDB", "PANTHER_PATHWAY"],
            "fallback_info_sources": ["GO", "GO-PLUS", "PW", "METACYC"],
            "default_hops": 1, "allow_multihop": True,
            "avoid_relations": ["co_mentioned", "no_relation"],
        },
        {
            "intent": "ppi_or_protein_complex",
            "keywords": ["ppi", "protein interaction", "interactome", "complex", "binding"],
            "source_type": "Protein", "target_type": "Protein",
            "relations": ["interact_with", "functionally_related_to"],
            "primary_info_sources": ["BioGRID_hs", "BioGRID_mm", "IntAct", "CORUM", "SIGNOR"],
            "fallback_info_sources": ["StringDB", "StingDB", "StingDB_mm", "PrimeKG"],
            "default_hops": 1, "allow_multihop": True,
            "avoid_relations": ["co_mentioned", "no_relation"],
        },
        {
            "intent": "regulatory_network",
            "keywords": ["regulate", "regulatory", "tf", "transcription factor", "scenic", "trrust", "phosphorylation", "signaling"],
            "source_type": "Gene", "target_type": "Gene",
            "relations": ["regulated", "causally_related_to", "interact_with"],
            "primary_info_sources": ["TRRUST", "SCENIC", "SCENIC_mm", "PhosphoSitePlus", "SIGNOR", "Reactome"],
            "fallback_info_sources": ["STITCH", "StringDB", "PrimeKG"],
            "default_hops": 1, "allow_multihop": True,
            "avoid_relations": ["co_mentioned", "no_relation"],
        },
        {
            "intent": "anatomy_or_trajectory",
            "keywords": ["tissue", "anatomy", "develop", "development", "trajectory", "differentiate", "differentiation", "spatial"],
            "source_type": "Cell", "target_type": "Cell",
            "relations": ["develope_into", "is_a", "has_part", "spatially_related_to"],
            "primary_info_sources": ["CL", "CLO", "UBERON", "EMAPA", "MA", "BTO"],
            "fallback_info_sources": ["EFO", "NCIT", "CellTaxonomy"],
            "default_hops": 1, "allow_multihop": True,
            "avoid_relations": ["co_mentioned", "no_relation"],
        },
    ],
}

QUALITY_WEIGHT = {
    "curated_clinical": 1.0, "curated_genetic_association": 0.95, "curated_biological": 0.92, "single_cell_curated": 0.90,
    "ontology": 0.82, "integrated_network": 0.68, "predicted_or_text_mined": 0.45, "unknown": 0.55,
}
TARGET_HINTS = {"cell":"cell_marker_annotation", "chemical":"drug_target_or_repositioning", "disease":"disease_gene_or_variant", "phenotype":"disease_gene_or_variant", "pathway":"pathway_function_enrichment", "process":"pathway_function_enrichment", "function":"pathway_function_enrichment", "cell_component":"pathway_function_enrichment", "protein":"ppi_or_protein_complex", "tissue":"anatomy_or_trajectory"}
RELATION_HINTS = {"marker_of":"cell_marker_annotation", "produce":"cell_marker_annotation", "pharmaco_effect_on":"drug_target_or_repositioning", "causally_related_to":"disease_gene_or_variant", "disease_relation":"disease_gene_or_variant", "participate_in":"pathway_function_enrichment", "functionally_related_to":"pathway_function_enrichment", "interact_with":"ppi_or_protein_complex", "regulated":"regulatory_network", "regulate":"regulatory_network", "develope_into":"anatomy_or_trajectory", "spatially_related_to":"anatomy_or_trajectory"}


def _default_path() -> Path:
    return Path(__file__).resolve().parents[1] / "references" / "kg_source_policy.json"


def load_policy_catalog(path: Optional[str] = None) -> Dict[str, Any]:
    if path:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    default = _default_path()
    if default.exists():
        try:
            return json.loads(default.read_text(encoding="utf-8"))
        except Exception:
            pass
    return DEFAULT_CATALOG


def split_csv(value: Optional[Any]) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        out: List[str] = []
        for item in value:
            out.extend(split_csv(item))
        return _dedupe(out)
    return [x.strip() for x in str(value).split(",") if x.strip()]


def _dedupe(items: Iterable[str]) -> List[str]:
    seen = set(); out = []
    for item in items:
        key = str(item).strip().casefold()
        if key and key not in seen:
            seen.add(key); out.append(str(item).strip())
    return out


def resolve_relation_aliases(
    relations: Optional[Any],
    schema_file: Optional[str] = None,
) -> tuple[List[str], List[str], Dict[str, List[str]]]:
    """Map concept-level relation names to KG relation names.

    Returns ``(resolved, unknown, alias_map)``. ``alias_map`` records only
    replacements, for example ``{"marker_of": ["produce"]}``.
    """
    schema = load_kg_schema(schema_file) if schema_file else _SCHEMA
    aliases = schema.get("relation_aliases", RELATION_ALIASES)
    queryable = set(schema.get("queryable_edge_types", []))
    weak = set(schema.get("weak_evidence_edge_types", []))
    blocked = set(schema.get("blocked_query_edge_types", ["no_relation"]))
    known = queryable | weak

    resolved: List[str] = []
    alias_map: Dict[str, List[str]] = {}
    for relation in split_csv(relations):
        key = _norm(relation).replace("-", "_")
        if key in aliases:
            replacements = list(aliases.get(key) or [])
            alias_map[relation] = replacements
            resolved.extend(replacements)
        else:
            resolved.append(relation)
    resolved = _dedupe([
        r for r in resolved
        if _norm(r).replace("-", "_") not in blocked
    ])
    unknown = [
        r for r in resolved
        if _norm(r).replace("-", "_") not in known
    ]
    return resolved, unknown, alias_map


def _norm(value: Any) -> str:
    return "" if value is None else str(value).strip().casefold()


def split_info_source_tokens(value: Optional[Any]) -> List[str]:
    """Split edge info_source metadata while preserving source names."""
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        raw_items = list(value)
    else:
        text = str(value).strip()
        if not text:
            return []
        raw_items: list[Any]
        try:
            parsed = ast.literal_eval(text)
        except Exception:
            parsed = None
        if isinstance(parsed, (list, tuple, set)):
            raw_items = list(parsed)
        else:
            raw_items = re.split(r"[,;|]", text)
    out: List[str] = []
    for item in raw_items:
        token = str(item).strip().strip("\"'")
        if not token or token == "Undef":
            continue
        out.append(token)
    return out


def canonicalize_info_source(value: Optional[Any]) -> str:
    """Canonical display name for known KG edge info_source variants."""
    tokens = split_info_source_tokens(value)
    if not tokens:
        return ""
    canonical: List[str] = []
    for token in tokens:
        key = _norm(token).replace("_", " ").replace("-", " ")
        key_compact = key.replace(" ", "_")
        mapped = INFO_SOURCE_CANONICAL_ALIASES.get(key) or INFO_SOURCE_CANONICAL_ALIASES.get(key_compact)
        canonical.append(mapped or token)
    return "|".join(_dedupe(canonical))


def canonicalize_info_sources(values: Optional[Any]) -> List[str]:
    return _dedupe([canonicalize_info_source(token) or token for token in split_csv(values)])


def expand_info_source_filter(values: Optional[Any]) -> List[str]:
    """Expand canonical source names to legacy KG source spellings for querying."""
    expanded: List[str] = []
    for token in split_csv(values):
        canonical = canonicalize_info_source(token) or token
        expanded.extend(INFO_SOURCE_QUERY_ALIASES.get(canonical, [canonical]))
    return _dedupe(expanded)


def _intent_cfg(catalog: Dict[str, Any], name: str) -> Dict[str, Any]:
    for item in catalog.get("intents", []):
        if _norm(item.get("intent")) == _norm(name):
            return item
    return {}


def infer_intent(question: Optional[str] = None, source_type: Optional[str] = None, target_type: Optional[str] = None, relation: Optional[Any] = None, task: Optional[str] = None, catalog: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    catalog = catalog or load_policy_catalog()
    text = " ".join(str(x) for x in [question or "", task or "", source_type or "", target_type or "", relation or ""] if x).casefold()
    scores: Dict[str, float] = {}; reasons: Dict[str, List[str]] = {}
    def add(intent: str, score: float, reason: str) -> None:
        if intent:
            scores[intent] = scores.get(intent, 0.0) + score
            reasons.setdefault(intent, []).append(reason)
    if task:
        tl = _norm(task)
        if "cell" in tl: add("cell_marker_annotation", 3, f"task={task}")
        if "drug" in tl: add("drug_target_or_repositioning", 3, f"task={task}")
        if "enrich" in tl and (not target_type or _norm(target_type) in {"pathway", "process", "function", "cell_component"}):
            add("pathway_function_enrichment", 3, f"task={task}")
        if "trajectory" in tl: add("anatomy_or_trajectory", 3, f"task={task}")
    add(TARGET_HINTS.get(_norm(target_type), ""), 2.0, f"target_type={target_type}")
    for rel in split_csv(relation):
        add(RELATION_HINTS.get(_norm(rel), ""), 2.5, f"relation={rel}")
    for item in catalog.get("intents", []):
        name = item.get("intent", "")
        if _norm(source_type) and _norm(source_type) == _norm(item.get("source_type")): add(name, 0.3, f"source_type={source_type}")
        if _norm(target_type) and _norm(target_type) == _norm(item.get("target_type")): add(name, 0.6, f"target_type matches")
        for kw in item.get("keywords", []):
            if str(kw).casefold() in text: add(name, 1.0, f"keyword={kw}")
    if not scores:
        return {"intent": "generic_kg_lookup", "confidence": 0.0, "reason": "No strong intent match."}
    best, score = sorted(scores.items(), key=lambda x: (-x[1], x[0]))[0]
    return {"intent": best, "confidence": round(min(0.99, score/5), 3), "reason": "; ".join(reasons.get(best, [])[:5])}


def source_quality_tier(info_source: Any, catalog: Optional[Dict[str, Any]] = None) -> str:
    catalog = catalog or load_policy_catalog()
    tokens = split_info_source_tokens(info_source)
    canonical_tokens = [canonicalize_info_source(token) or token for token in tokens]
    if not canonical_tokens:
        return "unknown"
    text = "|".join(_norm(token) for token in canonical_tokens if token)
    pieces = [_norm(token) for token in canonical_tokens if token]
    for tier, names in catalog.get("source_quality_tiers", {}).items():
        names_l = {_norm(x) for x in names}
        if text in names_l or any(p in names_l for p in pieces):
            return tier
    for tier, names in catalog.get("source_quality_tiers", {}).items():
        for name in names:
            nl = _norm(name)
            if nl and (nl in text or any(nl in p for p in pieces)):
                return tier
    return "unknown"


def source_quality_weight(info_source: Any, catalog: Optional[Dict[str, Any]] = None) -> float:
    return QUALITY_WEIGHT.get(source_quality_tier(info_source, catalog), QUALITY_WEIGHT["unknown"])


def select_kg_policy(question: Optional[str] = None, source_type: Optional[str] = None, target_type: Optional[str] = None, relation: Optional[Any] = None, info_source: Optional[Any] = None, hops: Optional[int] = None, task: Optional[str] = None, kg_policy: str = "auto", policy_file: Optional[str] = None, intent: Optional[str] = None, schema_file: Optional[str] = None, **legacy_kwargs: Any) -> Dict[str, Any]:
    # Compatibility with earlier experimental interface.
    if "mode" in legacy_kwargs and kg_policy == "auto":
        kg_policy = legacy_kwargs["mode"]
    if "user_info_sources" in legacy_kwargs and info_source is None:
        info_source = legacy_kwargs["user_info_sources"]
    if "requested_hops" in legacy_kwargs and hops is None:
        hops = legacy_kwargs["requested_hops"]
    if kg_policy == "manual":
        kg_policy = "direct"
    if kg_policy not in KG_POLICY_CHOICES:
        raise ValueError(f"Unsupported kg_policy={kg_policy!r}")
    catalog = load_policy_catalog(policy_file)
    inferred = infer_intent(question, source_type, target_type, relation, task, catalog)
    cfg = _intent_cfg(catalog, intent or inferred["intent"])
    raw_relations = split_csv(relation)
    explicit_relations, unknown_relations, relation_aliases = resolve_relation_aliases(relation, schema_file=schema_file)
    explicit_sources = canonicalize_info_sources(info_source)
    catalog_relations, catalog_unknown, catalog_aliases = resolve_relation_aliases(cfg.get("relations", []), schema_file=schema_file)
    relations_explicitly_dropped = bool(raw_relations and relation_aliases and not explicit_relations)
    relations = explicit_relations if (explicit_relations or relations_explicitly_dropped) else catalog_relations
    unknown_relations = unknown_relations or catalog_unknown
    relation_aliases = relation_aliases or catalog_aliases
    sources = explicit_sources or list(cfg.get("primary_info_sources", cfg.get("info_sources", [])))
    fallback = list(cfg.get("fallback_info_sources", []))
    use_kg = kg_policy != "off"
    req_hops = int(hops) if hops is not None else int(cfg.get("default_hops", 1) or 1)
    if kg_policy == "multihop": req_hops = max(2, req_hops)
    max_hops = max(1, min(MAX_HOPS, req_hops))
    confidence = float(inferred.get("confidence", 0.0) or 0.0)
    has_user_scope = bool(target_type or relation or task or intent or info_source)
    if kg_policy == "auto" and confidence < MIN_AUTO_INTENT_CONFIDENCE and not has_user_scope:
        return {
            "use_kg": use_kg,
            "kg_policy": "manual_required",
            "intent": "ambiguous",
            "intent_confidence": confidence,
            "source_type": source_type,
            "target_type": target_type,
            "relations": [],
            "primary_info_sources": [],
            "info_sources": [],
            "fallback_info_sources": [],
            "source_quality_tiers": {},
            "max_hops": 1,
            "hops": 1,
            "allow_multihop": False,
            "avoid_relations": ["co_mentioned", "no_relation"],
            "required_edge_fields": catalog.get("default_edge_fields", []),
            "co_mentioned_policy": catalog.get("co_mentioned_policy", {}),
            "broad_query": True,
            "broad_query_allowed": False,
            "review_required": True,
            "reason": "Question/context is too ambiguous for automatic KG source selection.",
            "warnings": [
                "Low-confidence KG routing. Add --tt, -r, --task, --intent, or --info-source.",
            ],
            "unknown_relations": [],
            "relation_aliases": {},
        }
    if kg_policy == "direct" and not explicit_sources:
        sources = []
    if kg_policy in {"source-scoped", "multihop"} and not sources:
        sources = fallback
    broad = not relations and not sources and not target_type
    reason_bits = [inferred.get("reason", "")]
    if sources: reason_bits.append("edge info_source scoped")
    if max_hops > 1: reason_bits.append(f"limited multi-hop enabled ({max_hops})")
    if broad: reason_bits.append("broad query risk")
    warnings = []
    if broad:
        warnings.append("Broad query risk: add --tt, -r, --info-source, or explicit target.")
    if relation_aliases:
        warnings.append(f"Resolved relation aliases: {relation_aliases}")
    if unknown_relations:
        warnings.append(f"Unknown KG relation(s): {unknown_relations}")
    return {
        "use_kg": use_kg, "kg_policy": kg_policy, "intent": intent or inferred["intent"], "intent_confidence": confidence,
        "source_type": source_type or cfg.get("source_type"), "target_type": target_type or cfg.get("target_type"),
        "relations": relations, "primary_info_sources": sources, "info_sources": sources, "fallback_info_sources": fallback,
        "source_quality_tiers": {s: source_quality_tier(s, catalog) for s in sources}, "max_hops": max_hops, "hops": max_hops,
        "allow_multihop": bool(cfg.get("allow_multihop", max_hops > 1)), "avoid_relations": list(cfg.get("avoid_relations", ["co_mentioned", "no_relation"])),
        "required_edge_fields": catalog.get("default_edge_fields", []), "co_mentioned_policy": catalog.get("co_mentioned_policy", {}),
        "broad_query": broad, "broad_query_allowed": False, "review_required": True,
        "reason": "; ".join(x for x in reason_bits if x) or "Use explicit filters or conservative direct lookup.",
        "warnings": warnings,
        "unknown_relations": unknown_relations,
        "relation_aliases": relation_aliases,
    }


def query_options_from_policy(policy: Dict[str, Any], include_relations: bool = True, include_hops: bool = True) -> List[str]:
    opts: List[str] = []
    if include_relations and policy.get("relations"):
        opts += ["-r", ",".join(policy["relations"])]
    if policy.get("primary_info_sources"):
        opts += ["--info-source", ",".join(policy["primary_info_sources"])]
    if include_hops and int(policy.get("max_hops", 1) or 1) > 1:
        opts += ["--hops", str(policy["max_hops"])]
    return opts


def render_query_command(source: Optional[str] = None, input_file: Optional[str] = None, output: Optional[str] = None, source_type: Optional[str] = None, target_type: Optional[str] = None, limit: Optional[Any] = None, policy: Optional[Dict[str, Any]] = None) -> List[str]:
    cmd = ["brick", "qr"]
    if input_file: cmd += ["-i", input_file]
    elif source: cmd.append(source)
    if source_type or (policy and policy.get("source_type")): cmd += ["--st", source_type or policy.get("source_type")]
    if target_type or (policy and policy.get("target_type")): cmd += ["--tt", target_type or policy.get("target_type")]
    if policy: cmd += query_options_from_policy(policy)
    if limit is not None: cmd += ["-n", str(limit)]
    if output: cmd += ["-o", output]
    return [str(x) for x in cmd]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="brick policy", description="Plan controlled KG source/hop query policies")
    sub = parser.add_subparsers(dest="command", required=True)
    route = sub.add_parser("route", aliases=["route-query"], help="Infer KG policy from question/query context")
    route.add_argument("-q", "--question")
    route.add_argument("--task")
    route.add_argument("--intent")
    route.add_argument("--source-type", "--st", dest="source_type")
    route.add_argument("--target-type", "--tt", dest="target_type")
    route.add_argument("-r", "--relation")
    route.add_argument("--info-source")
    route.add_argument("--hops", type=int, choices=range(1, MAX_HOPS+1))
    route.add_argument("--kg-policy", "--mode", dest="kg_policy", choices=KG_POLICY_CHOICES, default="auto")
    route.add_argument("--policy-file")
    route.add_argument("--kg-schema", default=None, help="Path to KG schema JSON. Defaults to references/kg_schema.json.")
    route.add_argument("--format", choices=["json", "table", "sh"], default="json")
    route.add_argument("-o", "--output")
    route.add_argument("--command", dest="emit_command", action="store_true")
    route.add_argument("--source")
    route.add_argument("--input-file")
    route.add_argument("--output-csv")
    route.add_argument("-n", "--limit", default="500")
    listed = sub.add_parser("list", help="List known KG policy intents")
    listed.add_argument("--policy-file")
    listed.add_argument("--format", choices=["json", "table"], default="table")
    export = sub.add_parser("export-default", aliases=["export"], help="Export the built-in KG policy catalog JSON")
    export.add_argument("-o", "--output")
    return parser


def _print_table(policy: Dict[str, Any]) -> None:
    for key in ["use_kg", "kg_policy", "intent", "intent_confidence", "source_type", "target_type", "relations", "primary_info_sources", "fallback_info_sources", "max_hops", "broad_query", "reason"]:
        print(f"{key}: {policy.get(key)}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command in {"export-default", "export"}:
            text = json.dumps(DEFAULT_CATALOG, indent=2, ensure_ascii=False) + "\n"
            if args.output:
                Path(args.output).write_text(text, encoding="utf-8")
                print(f"Wrote {args.output}")
            else:
                print(text, end="")
            return 0
        if args.command == "list":
            data = load_policy_catalog(args.policy_file)
            intents = data.get("intents", [])
            if args.format == "json": print(json.dumps(intents, indent=2, ensure_ascii=False))
            else:
                for item in intents:
                    print(f"{item.get('intent')}: relations={','.join(item.get('relations', [])[:5])} info_sources={','.join(item.get('primary_info_sources', [])[:5])} hops={item.get('default_hops', 1)}")
            return 0
        policy = select_kg_policy(args.question, args.source_type, args.target_type, args.relation, args.info_source, args.hops, args.task, args.kg_policy, args.policy_file, args.intent, schema_file=args.kg_schema)
        result: Dict[str, Any] = {"policy": policy}
        if getattr(args, "emit_command", False):
            cmd = render_query_command(args.source, args.input_file, args.output_csv, args.source_type, args.target_type, args.limit, policy)
            result["command"] = cmd; result["command_sh"] = " ".join(shlex.quote(x) for x in cmd)
        if args.output:
            Path(args.output).write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8"); print(f"Wrote {args.output}"); return 0
        if args.format == "json": print(json.dumps(result, indent=2, ensure_ascii=False))
        elif args.format == "sh": print(result.get("command_sh") or " ".join(shlex.quote(x) for x in query_options_from_policy(policy)))
        else:
            _print_table(policy)
            if "command_sh" in result: print("command: " + result["command_sh"])
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


def _cli_main() -> None:
    raise SystemExit(main())


if __name__ == "__main__":
    raise SystemExit(main())
