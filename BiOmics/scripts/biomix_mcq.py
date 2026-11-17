#!/usr/bin/env python3
"""Biomix MCQ helpers for question-only KG option scoring."""

from __future__ import annotations

import re
from collections.abc import Iterable


def _clean(text: str) -> str:
    return re.sub(r"\s+", " ", str(text or "").strip())


def _add(out: list[str], value: str) -> None:
    value = _clean(value)
    key = value.casefold()
    if value and key not in {item.casefold() for item in out}:
        out.append(value)


def expand_disease_target(term: str) -> list[str]:
    """Return conservative disease-name aliases derived only from question text.

    The aliases are generic lexical/biomedical surface variants. They are used
    for KG target matching and evidence scoring, never for looking at gold
    answers.
    """

    base = _clean(term)
    if not base:
        return []
    aliases: list[str] = []
    _add(aliases, base)

    lower = base.casefold()
    no_apostrophe = re.sub(r"[']", "", base)
    if no_apostrophe != base:
        _add(aliases, no_apostrophe)

    hyphen_space = base.replace("-", " ")
    if hyphen_space != base:
        _add(aliases, hyphen_space)
    space_hyphen = re.sub(r"\s+", "-", base)
    if space_hyphen != base and " " in base:
        _add(aliases, space_hyphen)

    swaps = [
        (" carcinoma", " cancer"),
        (" cancer", " carcinoma"),
        (" cancer", " neoplasm"),
        (" tumor", " cancer"),
        (" tumour", " cancer"),
        (" neoplasm", " cancer"),
        (" syndrome", " syndromes"),
        (" syndromes", " syndrome"),
    ]
    for old, new in swaps:
        if lower.endswith(old):
            _add(aliases, base[: -len(old)] + new)

    if lower.startswith("cervix "):
        rest = base.split(" ", 1)[1]
        _add(aliases, f"cervical {rest}")
        if rest.casefold() == "carcinoma":
            _add(aliases, "cervical cancer")
            _add(aliases, "cervical carcinoma")
    if lower.startswith("lung small cell"):
        suffix = base.split("lung small cell", 1)[1].strip()
        _add(aliases, f"small cell lung {suffix}")
        _add(aliases, "small cell lung cancer")
        _add(aliases, "small cell lung carcinoma")
    if lower.startswith("oral cavity "):
        rest = base.split(" ", 2)[2]
        _add(aliases, f"oral {rest}")
    if lower.startswith("gastric fundus "):
        rest = base.split(" ", 2)[2]
        _add(aliases, f"gastric {rest}")
        _add(aliases, f"stomach {rest}")
    if lower == "liver cirrhosis":
        _add(aliases, "cirrhosis of liver")
        _add(aliases, "hepatic cirrhosis")
        _add(aliases, "cirrhosis")
    if lower == "membranous glomerulonephritis":
        _add(aliases, "membranous nephropathy")
    if "sjogren" in lower:
        _add(aliases, "Sjogren syndrome")
        _add(aliases, "Sjogren's syndrome")
    if lower == "triple-receptor negative breast cancer":
        _add(aliases, "triple-negative breast cancer")
        _add(aliases, "triple negative breast cancer")
        _add(aliases, "triple receptor-negative breast cancer")
    if lower == "male breast cancer":
        _add(aliases, "breast cancer")
        _add(aliases, "male breast carcinoma")
    if lower == "androgenic alopecia":
        _add(aliases, "androgenetic alopecia")
    if lower == "dengue hemorrhagic fever":
        _add(aliases, "dengue haemorrhagic fever")
        _add(aliases, "dengue fever")
    if lower == "intracranial aneurysm":
        _add(aliases, "cerebral aneurysm")
        _add(aliases, "brain aneurysm")
    if lower == "cocaine dependence":
        _add(aliases, "cocaine use disorder")
        _add(aliases, "cocaine addiction")
    if lower == "atrophic gastritis":
        _add(aliases, "chronic atrophic gastritis")
    if lower == "myeloproliferative neoplasm":
        _add(aliases, "myeloproliferative neoplasms")
        _add(aliases, "myeloproliferative disorder")
    if lower == "head and neck cancer":
        _add(aliases, "head and neck carcinoma")
        _add(aliases, "head and neck neoplasm")
        _add(aliases, "head and neck squamous cell carcinoma")

    return aliases[:12]


def build_target_alias_groups(targets: Iterable[str]) -> dict[str, list[str]]:
    return {_clean(target): expand_disease_target(target) for target in targets if _clean(target)}


def flatten_target_alias_groups(groups: dict[str, list[str]]) -> list[str]:
    out: list[str] = []
    for aliases in groups.values():
        for alias in aliases:
            _add(out, alias)
    return out
