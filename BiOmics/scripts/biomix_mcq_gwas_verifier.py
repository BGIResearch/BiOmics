#!/usr/bin/env python3
"""Biomix MCQ GWAS Catalog option verifier.

Generation reads only Biomix MCQ questions/options. It does not read gold
labels, Biomni outputs, consensus outputs, or previous answer files. Each option
is checked as a hypothesis against GWAS Catalog associations for the two traits
named in the question; low-evidence cases fall back to deepseek-v4-pro.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import requests
from langchain_core.messages import HumanMessage, SystemMessage
from langchain_openai import ChatOpenAI


GWAS_BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2"
CHAT_BASE_URL = "https://dcsapi.dcs.cloud/api/aigress/unified/v1"
CHAT_API_KEY_ENV = "BRICK_API_KEY"
DEFAULT_METHOD = "deepseek_v4_pro_D_gwas_sourceaware_20260707"
_THREAD_LOCAL = threading.local()
_CACHE_LOCK = threading.Lock()


def default_benchmark_root() -> Path:
    repo_from_skill = Path(__file__).resolve().parents[3]
    for candidate in [Path.cwd() / "biomics-benchmark", Path.cwd(), repo_from_skill]:
        if (candidate / "dataset" / "biomix" / "mcq_biomix_pure_question.json").exists():
            return candidate.resolve()
    return repo_from_skill.resolve()


def read_json(path: Path, default: Any = None) -> Any:
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def load_local_env_var(key: str) -> str:
    value = os.environ.get(key)
    if value:
        return value.strip()
    candidates: list[Path] = []
    local_appdata = os.environ.get("LOCALAPPDATA")
    if local_appdata:
        candidates.append(Path(local_appdata) / "hermes" / ".env")
    candidates.append(Path.home() / "AppData" / "Local" / "hermes" / ".env")
    for env_path in candidates:
        if not env_path.exists():
            continue
        for line in env_path.read_text(encoding="utf-8-sig").splitlines():
            if line.startswith(f"{key}="):
                value = line.split("=", 1)[1].strip().strip('"').strip("'")
                if value:
                    os.environ[key] = value
                    return value
    raise RuntimeError(f"Missing API key env var: {key}")


def normalize(text: Any) -> str:
    return re.sub(r"\s+", " ", str(text or "").strip().strip("\"'`.,;:()[]{}")).casefold()


def dedupe(values: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = re.sub(r"\s+", " ", str(value or "").strip())
        key = text.casefold()
        if text and key not in seen:
            seen.add(key)
            out.append(text)
    return out


def parse_mcq(question: str) -> dict[str, Any]:
    match = re.search(
        r"which\s+(?P<kind>Gene|Variant|Mutation|Protein|Disease|Chemical)\s+is\s+associated\s+with\s+"
        r"(?P<targets>.+?)\.\s+Given list is:\s*(?P<opts>.+)$",
        question,
        flags=re.I | re.S,
    )
    if not match:
        return {"kind": "", "targets": [], "options": []}
    options = [item.strip(" .;:,[]()\"'`") for item in match.group("opts").split(",") if item.strip(" .;:,[]()\"'`")]
    target_text = match.group("targets").strip()
    if " and " in target_text:
        t1, t2 = target_text.rsplit(" and ", 1)
        targets = [t1.strip(), t2.strip()]
    else:
        targets = [target_text]
    return {
        "kind": match.group("kind"),
        "targets": targets,
        "options": options,
    }


def trait_aliases(trait: str) -> list[str]:
    t = re.sub(r"\s+", " ", trait.strip())
    low = t.casefold()
    aliases = [t]
    replacements = [
        (" cancer", " carcinoma"),
        (" carcinoma", " cancer"),
        ("tumour", "tumor"),
        ("tumor", "cancer"),
    ]
    for old, new in replacements:
        if old.strip() in low:
            aliases.append(re.sub(re.escape(old), new, t, flags=re.I))
    if "male breast cancer" in low:
        aliases.extend(["male breast carcinoma", "breast carcinoma", "breast cancer"])
    if "triple-receptor negative breast cancer" in low or "triple receptor negative breast cancer" in low:
        aliases.extend(["triple negative breast cancer", "estrogen-receptor negative breast cancer", "breast carcinoma"])
    if "estrogen-receptor negative breast cancer" in low or "estrogen receptor negative breast cancer" in low:
        aliases.extend(["estrogen receptor-negative breast cancer", "triple negative breast cancer", "breast carcinoma"])
    if "kuhnt-junius" in low:
        aliases.extend(["age related macular degeneration", "age-related macular degeneration", "macular degeneration"])
    if "sjogren" in low:
        aliases.extend(["Sjogren syndrome", "Sjögren syndrome", "Sjogren's syndrome"])
    if "non-alcoholic" in low:
        aliases.append(t.replace("non-alcoholic", "nonalcoholic"))
    if "head and neck" in low:
        aliases.extend(["head and neck carcinoma", "head and neck neoplasm"])
    if "oral cavity cancer" in low:
        aliases.extend(["oral cavity carcinoma", "mouth cancer", "oral cancer"])
    if "oropharynx cancer" in low:
        aliases.extend(["oropharyngeal cancer", "oropharyngeal carcinoma"])
    if "cervix carcinoma" in low:
        aliases.extend(["cervical carcinoma", "cervical cancer"])
    if "liver cirrhosis" in low:
        aliases.extend(["cirrhosis", "hepatic cirrhosis"])
    if "membranous glomerulonephritis" in low:
        aliases.extend(["membranous nephropathy", "glomerulonephritis"])
    if "collagenous colitis" in low:
        aliases.extend(["microscopic colitis", "colitis"])
    return dedupe(aliases)


def association_text(assoc: dict[str, Any]) -> str:
    parts: list[str] = []
    for trait in assoc.get("efo_traits") or []:
        parts.append(str(trait.get("efo_trait", "")))
        parts.append(str(trait.get("efo_id", "")))
    for trait in assoc.get("reported_trait") or []:
        parts.append(str(trait))
    parts.extend(str(gene) for gene in assoc.get("mapped_genes") or [])
    for allele in assoc.get("snp_allele") or []:
        parts.append(str(allele.get("rs_id", "")))
    return normalize(" ".join(parts))


def trait_supported(assoc: dict[str, Any], trait: str) -> bool:
    text = association_text(assoc)
    aliases = trait_aliases(trait)
    for alias in aliases:
        alias_norm = normalize(alias)
        if alias_norm and alias_norm in text:
            return True
    tokens = [tok for tok in re.findall(r"[a-z0-9]+", normalize(trait)) if len(tok) > 3]
    return bool(tokens) and sum(tok in text for tok in tokens) >= max(1, min(2, len(tokens)))


def gwas_session() -> requests.Session:
    session = getattr(_THREAD_LOCAL, "gwas_session", None)
    if session is None:
        session = requests.Session()
        session.headers.update({"Accept": "application/json", "User-Agent": "BiOmics-BRICK-benchmark/strict-mcq"})
        setattr(_THREAD_LOCAL, "gwas_session", session)
    return session


def load_cache(path: Path) -> dict[str, Any]:
    data = read_json(path, default={})
    return data if isinstance(data, dict) else {}


def save_cache(path: Path, cache: dict[str, Any]) -> None:
    with _CACHE_LOCK:
        write_json(path, cache)


def query_gwas(args: argparse.Namespace, cache: dict[str, Any], option: str, option_kind: str, trait_alias: str) -> list[dict[str, Any]]:
    option = option.strip()
    trait_alias = trait_alias.strip()
    key = json.dumps([option_kind, option, trait_alias], ensure_ascii=False)
    with _CACHE_LOCK:
        if key in cache:
            return list(cache[key])

    params: dict[str, Any] = {"efo_trait": trait_alias, "size": args.gwas_size}
    if option_kind.casefold() in {"variant", "mutation"} or option.casefold().startswith("rs"):
        params["rs_id"] = option
    else:
        params["mapped_gene"] = option

    last_error = ""
    associations: list[dict[str, Any]] = []
    for attempt in range(args.gwas_retries + 1):
        try:
            response = gwas_session().get(f"{GWAS_BASE_URL}/associations", params=params, timeout=args.gwas_timeout)
            response.raise_for_status()
            payload = response.json()
            associations = payload.get("_embedded", {}).get("associations", []) or []
            break
        except Exception as exc:
            last_error = repr(exc)
            if attempt < args.gwas_retries:
                time.sleep(0.5 * (attempt + 1))
    if last_error and not associations:
        associations = [{"_error": last_error}]

    with _CACHE_LOCK:
        cache[key] = associations
    return associations


def assoc_strength(assoc: dict[str, Any]) -> float:
    if assoc.get("_error"):
        return 0.0
    p_value = assoc.get("p_value")
    try:
        p_score = min(max(-math.log10(float(p_value)), 0.0), 12.0) / 12.0
    except Exception:
        p_score = 0.25
    mapped_bonus = 0.25 if assoc.get("mapped_genes") else 0.0
    return 1.0 + p_score + mapped_bonus


def score_option(args: argparse.Namespace, cache: dict[str, Any], option: str, option_kind: str, targets: list[str]) -> dict[str, Any]:
    target_details: list[dict[str, Any]] = []
    total = 0.0
    covered = 0
    for target in targets:
        best = {"score": 0.0, "alias": "", "hits": []}
        for alias in trait_aliases(target)[: args.max_trait_aliases]:
            hits = [
                assoc for assoc in query_gwas(args, cache, option, option_kind, alias)
                if trait_supported(assoc, target)
            ]
            if hits:
                score = sum(assoc_strength(hit) for hit in hits[: args.max_hits_per_alias])
                if score > best["score"]:
                    best = {
                        "score": round(score, 4),
                        "alias": alias,
                        "hits": [
                            {
                                "efo_traits": hit.get("efo_traits", []),
                                "reported_trait": hit.get("reported_trait", []),
                                "mapped_genes": hit.get("mapped_genes", []),
                                "snp_allele": hit.get("snp_allele", []),
                                "p_value": hit.get("p_value"),
                                "pubmed_id": hit.get("pubmed_id"),
                                "accession_id": hit.get("accession_id"),
                            }
                            for hit in hits[:3]
                        ],
                    }
        if best["score"] > 0:
            covered += 1
        total += float(best["score"])
        target_details.append({"target": target, **best})
    both_bonus = args.both_target_bonus if covered == len(targets) and targets else 0.0
    return {
        "option": option,
        "coverage": covered,
        "score": round(total + both_bonus, 4),
        "raw_score": round(total, 4),
        "both_target_bonus": both_bonus,
        "targets": target_details,
    }


def get_llm(args: argparse.Namespace) -> ChatOpenAI:
    key = f"llm_{args.model}_{args.max_tokens}_{args.temperature}"
    llm = getattr(_THREAD_LOCAL, key, None)
    if llm is None:
        llm = ChatOpenAI(
            model=args.model,
            base_url=args.base_url,
            api_key=load_local_env_var(args.api_key_env),
            temperature=args.temperature,
            timeout=args.request_timeout,
            max_retries=1,
            max_completion_tokens=args.max_tokens,
            extra_body={"thinking": {"type": "disabled"}},
        )
        setattr(_THREAD_LOCAL, key, llm)
    return llm


def clean_choice(raw: str, options: list[str]) -> str:
    text = re.sub(r"\s+", " ", str(raw or "").strip().strip("\"'`"))
    low = normalize(text)
    for option in options:
        if normalize(option) == low:
            return option
    hits = [option for option in options if re.search(rf"(?<![A-Za-z0-9_]){re.escape(option)}(?![A-Za-z0-9_])", text)]
    return hits[0] if len(set(hits)) == 1 else ""


def native_fallback(args: argparse.Namespace, question: str, options: list[str]) -> dict[str, Any]:
    messages = [
        SystemMessage(content="You are a biomedical GWAS/disease-association MCQ solver. Return only one option exactly as listed."),
        HumanMessage(
            content=(
                "Select the option most specifically associated with BOTH traits in the question. "
                "Avoid generic hub genes unless they are specifically linked to both traits.\n"
                f"Question: {question}\nOptions: {', '.join(options)}\nReturn only the exact option text."
            )
        ),
    ]
    try:
        raw = str(get_llm(args).invoke(messages).content or "")
        return {"answer": clean_choice(raw, options), "raw": raw, "error": ""}
    except Exception as exc:
        return {"answer": "", "raw": "", "error": repr(exc)}


def choose_answer(args: argparse.Namespace, option_scores: list[dict[str, Any]], fallback_answer: str) -> tuple[str, dict[str, Any]]:
    ranked = sorted(option_scores, key=lambda x: (x["coverage"], x["score"]), reverse=True)
    if not ranked:
        return fallback_answer, {"source": "native_fallback_no_gwas_scores"}
    best = ranked[0]
    second = ranked[1] if len(ranked) > 1 else {"coverage": 0, "score": 0.0}
    coverage_margin = int(best["coverage"]) - int(second.get("coverage", 0))
    score_margin = float(best["score"]) - float(second.get("score", 0.0))
    if best["coverage"] >= args.min_gwas_coverage and (
        coverage_margin >= args.min_coverage_margin or score_margin >= args.min_score_margin
    ):
        return best["option"], {
            "source": "gwas_catalog_option_verifier",
            "best": best,
            "second": second,
            "coverage_margin": coverage_margin,
            "score_margin": round(score_margin, 4),
        }
    if fallback_answer:
        return fallback_answer, {
            "source": "deepseek_native_fallback_low_gwas_margin",
            "best": best,
            "second": second,
            "coverage_margin": coverage_margin,
            "score_margin": round(score_margin, 4),
        }
    return best["option"], {
        "source": "gwas_best_effort_no_native",
        "best": best,
        "second": second,
        "coverage_margin": coverage_margin,
        "score_margin": round(score_margin, 4),
    }


def run_one(args: argparse.Namespace, cache: dict[str, Any], row: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    qid = int(row["id"])
    question = str(row["question"])
    parsed = parse_mcq(question)
    options = parsed["options"]
    targets = parsed["targets"]
    kind = parsed["kind"]
    option_scores = [score_option(args, cache, option, kind, targets) for option in options]
    fallback = native_fallback(args, question, options) if args.native_fallback else {"answer": "", "raw": "", "error": ""}
    answer, decision = choose_answer(args, option_scores, fallback.get("answer", ""))
    detail = {
        "id": qid,
        "question": question,
        "kind": kind,
        "targets": targets,
        "options": options,
        "answer": answer,
        "decision": decision,
        "option_scores": option_scores,
        "native_fallback": fallback,
        "gold_access_during_generation": "none",
        "forbidden_sources": ["gold", "Biomni", "previous answer files", "consensus files"],
        "evidence_source": "GWAS Catalog REST API v2",
        "model": args.model,
    }
    return {"id": qid, "answer": answer}, detail


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run strict Biomix MCQ GWAS option verification.")
    parser.add_argument("--benchmark-root", type=Path, default=default_benchmark_root())
    parser.add_argument("--method", default=DEFAULT_METHOD)
    parser.add_argument("--model", default="deepseek-v4-pro")
    parser.add_argument("--base-url", default=CHAT_BASE_URL)
    parser.add_argument("--api-key-env", default=CHAT_API_KEY_ENV)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--native-fallback", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--request-timeout", type=float, default=45)
    parser.add_argument("--max-tokens", type=int, default=128)
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--gwas-timeout", type=float, default=20)
    parser.add_argument("--gwas-retries", type=int, default=1)
    parser.add_argument("--gwas-size", type=int, default=20)
    parser.add_argument("--max-trait-aliases", type=int, default=8)
    parser.add_argument("--max-hits-per-alias", type=int, default=5)
    parser.add_argument("--both-target-bonus", type=float, default=4.0)
    parser.add_argument("--min-gwas-coverage", type=int, default=1)
    parser.add_argument("--min-coverage-margin", type=int, default=1)
    parser.add_argument("--min-score-margin", type=float, default=1.0)
    parser.add_argument("--cache-file", type=Path, default=None)
    return parser.parse_args()


def _cli_main() -> int:
    args = parse_args()
    args.benchmark_root = args.benchmark_root.resolve()
    input_path = args.benchmark_root / "dataset" / "biomix" / "mcq_biomix_pure_question.json"
    answer_dir = args.benchmark_root / "agent_answer" / "biomix" / args.method
    detail_dir = args.benchmark_root / "run_agent_pipeline" / "brickkg" / "summaries" / args.method
    output_path = answer_dir / "mcq_biomix_pure_question.json"
    detail_path = detail_dir / "mcq_biomix_pure_question.details.json"
    manifest_path = detail_dir / "gwas_option_verifier_manifest.json"
    args.cache_file = args.cache_file or detail_dir / "gwas_cache.json"
    if output_path.exists() and not args.overwrite:
        print(f"exists: {output_path}", flush=True)
        return 0
    rows = read_json(input_path, default=[])
    if args.limit is not None:
        rows = rows[: args.limit]
    cache = load_cache(args.cache_file)
    print(f"gwas-mcq: loaded={len(rows)} cache={len(cache)} model={args.model}", flush=True)

    results: dict[int, dict[str, Any]] = {}
    details: dict[str, Any] = {}
    failures = 0
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {pool.submit(run_one, args, cache, row): int(row["id"]) for row in rows}
        for idx, future in enumerate(as_completed(futures), start=1):
            qid = futures[future]
            try:
                result, detail = future.result()
                results[qid] = result
                details[str(qid)] = detail
                decision_source = detail.get("decision", {}).get("source")
                print(f"[OK] id={qid} answer={result.get('answer')} source={decision_source}", flush=True)
            except Exception as exc:
                failures += 1
                results[qid] = {"id": qid, "answer": ""}
                details[str(qid)] = {"id": qid, "error": repr(exc)}
                print(f"[ERR] id={qid} {exc!r}", file=sys.stderr, flush=True)
            if idx % 20 == 0 or idx == len(rows):
                ordered = [results[int(row["id"])] for row in rows if int(row["id"]) in results]
                write_json(output_path, ordered)
                write_json(detail_path, details)
                save_cache(args.cache_file, cache)
                print(f"checkpoint {idx}/{len(rows)} failed={failures} cache={len(cache)}", flush=True)
    write_json(
        manifest_path,
        {
            "method": args.method,
            "model": args.model,
            "mode": "strict_question_only_gwas_option_verifier",
            "question_file": str(input_path),
            "answers": str(output_path),
            "details": str(detail_path),
            "cache_file": str(args.cache_file),
            "evidence_source": "GWAS Catalog REST API v2",
            "evidence_source_url": GWAS_BASE_URL,
            "gold_access_during_generation": "none",
            "forbidden_sources": ["gold", "Biomni", "previous answer files", "consensus files"],
            "native_fallback": args.native_fallback,
            "thresholds": {
                "min_gwas_coverage": args.min_gwas_coverage,
                "min_coverage_margin": args.min_coverage_margin,
                "min_score_margin": args.min_score_margin,
            },
        },
    )
    print(f"answers: {answer_dir}", flush=True)
    print(f"details: {detail_dir}", flush=True)
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(_cli_main())
