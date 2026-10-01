from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable, Mapping

import pandas as pd

MASK_LEVELS = ("none", "taxonomy", "strict")
DEFAULT_TOKEN = "[TAXON]"
COMPOUND_MIN_LENGTH = 6

_WORD = re.compile(r"[A-Za-z]+(?:-[A-Za-z]+)*")
_ABBREVIATED = re.compile(r"\b([A-Z])\.\s?(?:[a-z]\.\s?)?([a-z]+)((?:\s+[a-z]+)?)")
_NEXT_WORD = re.compile(r"[\s-]*([A-Za-z]+)")


def text_column(level: str) -> str:
    if level not in MASK_LEVELS:
        raise ValueError(f"unknown mask level {level!r}; expected one of {MASK_LEVELS}")
    return "text" if level == "none" else f"text_{level}"


def derived_terms(name: str) -> set[str]:
    """Adjective and noun forms of a family-group name: Felidae -> felid, feline; Lutrinae -> lutrine."""
    low = name.lower()
    if low.endswith("idae"):
        return {low[:-4] + "id", low[:-4] + "ine"}
    if low.endswith("inae"):
        return {low[:-4] + "ine"}
    if low.endswith("formia"):
        return {low[:-2]}
    return set()


def plural_forms(noun: str) -> set[str]:
    forms = {noun, noun + "s"}
    if noun.endswith(("s", "x", "ch", "sh")):
        forms.add(noun + "es")
    if noun.endswith("f"):
        forms.add(noun[:-1] + "ves")
    if noun.endswith("fe"):
        forms.add(noun[:-2] + "ves")
    if len(noun) > 2 and noun.endswith("y") and noun[-2] not in "aeiou":
        forms.add(noun[:-1] + "ies")
    return forms


def head_nouns(common_names: Iterable[object]) -> set[str]:
    """Last word of each common name ("Short-eared Dog" -> dog), plus the last part of a hyphenated one."""
    nouns: set[str] = set()
    for cell in common_names:
        if not isinstance(cell, str):
            continue
        for name in cell.split("|"):
            words = _WORD.findall(name)
            if not words:
                continue
            last = words[-1].lower()
            nouns.add(last)
            nouns.add(last.split("-")[-1])
    return {n for n in nouns if len(n) >= 3}


@dataclass(frozen=True)
class MaskVocabulary:
    names: frozenset[str]              # genera and higher-rank names with their derived forms; masked wherever they occur
    compound_names: frozenset[str]     # long genus names; a word containing one (Propoecilogale) is masked whole
    epithets: frozenset[str]           # current species epithets; masked wherever they occur
    binomial_epithets: frozenset[str]  # current and synonym epithets; masked only after a genus or an abbreviated genus
    genus_initials: frozenset[str]
    head_nouns: frozenset[str]         # strict level only, plural forms included
    keep_terms: frozenset[str] = frozenset()
    keep_when_followed_by: Mapping[str, frozenset[str]] = field(default_factory=dict)


def build_vocabulary(
    species: pd.DataFrame,
    synonyms: pd.DataFrame | None = None,
    extra_terms: Iterable[str] = (),
    keep_terms: Iterable[str] = (),
    keep_when_followed_by: Mapping[str, Iterable[str]] | None = None,
) -> MaskVocabulary:
    genera = {str(g).lower() for g in species["genus"].dropna()}
    nouns: set[str] = set()
    for col in ("common_name", "other_common_names"):
        if col in species.columns:
            nouns |= head_nouns(species[col])
    epithets = {str(e).lower() for e in species["specific_epithet"].dropna()}
    binomial_epithets = set(epithets)
    names = set(genera)
    for col in ("order", "suborder", "family", "subfamily"):
        if col in species.columns:
            for value in species[col].dropna().astype(str):
                names.add(value.lower())
                names |= derived_terms(value)
    if synonyms is not None:
        for root in synonyms["root_name"].dropna().astype(str):
            if re.fullmatch(r"[A-Za-z]{3,}", root):
                binomial_epithets.add(root.lower())
        for combination in synonyms["normalized_original_combination"].dropna().astype(str):
            words = combination.split()
            # a historical genus spelled like a current common name (Hyena, Coati, Serval) is left to the strict level
            if words and re.fullmatch(r"[A-Z][a-z]{2,}", words[0]) and words[0].lower() not in nouns:
                genera.add(words[0].lower())
            binomial_epithets.update(w for w in words[1:] if re.fullmatch(r"[a-z]{3,}", w))
        names |= genera
    names |= {str(t).lower() for t in extra_terms}
    return MaskVocabulary(
        names=frozenset(names),
        compound_names=frozenset(g for g in genera if len(g) >= COMPOUND_MIN_LENGTH),
        epithets=frozenset(epithets),
        binomial_epithets=frozenset(binomial_epithets),
        genus_initials=frozenset(g[0].upper() for g in genera),
        head_nouns=frozenset(form for noun in nouns for form in plural_forms(noun)),
        keep_terms=frozenset(str(t).lower() for t in keep_terms),
        keep_when_followed_by={str(k).lower(): frozenset(str(v).lower() for v in vs) for k, vs in (keep_when_followed_by or {}).items()},
    )


def _alternation(terms: Iterable[str]) -> str:
    return "|".join(re.escape(t) for t in sorted(terms, key=lambda t: (-len(t), t)))


class Masker:
    def __init__(self, vocabulary: MaskVocabulary, token: str = DEFAULT_TOKEN):
        self.vocabulary = vocabulary
        self.token = token
        v = vocabulary
        self._names = re.compile(rf"\b(?P<name>(?i:{_alternation(v.names)}))(?P<plural>e?s)?\b(?P<tail>(?:\s+[a-z]+){{0,2}})") if v.names else None
        self._compounds = re.compile(rf"\b[A-Za-z]*(?i:{_alternation(v.compound_names)})[A-Za-z]*\b") if v.compound_names else None
        self._epithets = re.compile(rf"\b(?i:{_alternation(v.epithets)})s?\b") if v.epithets else None
        self._nouns = re.compile(rf"\b(?i:{_alternation(v.head_nouns)})\b") if v.head_nouns else None
        self._collapse = re.compile(rf"{re.escape(token)}(?:[\s-]+{re.escape(token)})+")

    def _kept(self, term: str, text: str, end: int) -> bool:
        term = term.lower()
        if term in self.vocabulary.keep_terms:
            return True
        followers = self.vocabulary.keep_when_followed_by.get(term)
        if followers:
            nxt = _NEXT_WORD.match(text, end)
            return bool(nxt and nxt.group(1).lower() in followers)
        return False

    def mask(self, text: str, level: str, counts: Counter | None = None) -> str:
        if level not in MASK_LEVELS:
            raise ValueError(f"unknown mask level {level!r}; expected one of {MASK_LEVELS}")
        if level == "none":
            return text
        v, token = self.vocabulary, self.token

        def note(kind: str, surface: str) -> None:
            if counts is not None:
                counts[(kind, surface.lower())] += 1

        def abbreviated(m: re.Match) -> str:
            initial, epithet, extra = m.group(1), m.group(2), m.group(3)
            if initial not in v.genus_initials or epithet not in v.binomial_epithets:
                return m.group(0)
            handed_back = "" if extra.strip() in v.binomial_epithets else extra
            note("abbreviated", m.group(0)[: len(m.group(0)) - len(handed_back)])
            return token + handed_back

        def name(m: re.Match) -> str:
            surface = m.group("name") + (m.group("plural") or "")
            if self._kept(surface, m.string, m.end("plural") if m.group("plural") else m.end("name")):
                return m.group(0)
            tail = m.group("tail")
            words = re.findall(r"\s+[a-z]+", tail)
            used = 0
            while used < len(words) and words[used].strip() in v.binomial_epithets:
                used += 1
            consumed = "".join(words[:used])
            note("binomial" if used else "name", surface + consumed)
            return token + tail[len(consumed):]

        def simple(kind: str):
            def replace(m: re.Match) -> str:
                if self._kept(m.group(0), m.string, m.end()):
                    return m.group(0)
                note(kind, m.group(0))
                return token
            return replace

        out = _ABBREVIATED.sub(abbreviated, text)
        # words after a name are matched as a possible epithet and handed back if they are not one,
        # so repeat until a name hiding in such a tail has been masked too
        while self._names:
            again = self._names.sub(name, out)
            if again == out:
                break
            out = again
        if self._compounds:
            out = self._compounds.sub(simple("compound"), out)
        if self._epithets:
            out = self._epithets.sub(simple("epithet"), out)
        if level == "strict" and self._nouns:
            out = self._nouns.sub(simple("head_noun"), out)
        return self._collapse.sub(token, out)

    def count(self, masked: str) -> int:
        return masked.count(self.token)


def term_leaks(text: str, term: str) -> bool:
    """True if `term` survives in `text`: as a whole word (plural allowed), or anywhere inside a word when it is long."""
    term = str(term).lower()
    if len(term) >= COMPOUND_MIN_LENGTH:
        return term in text.lower()
    return re.search(rf"\b{re.escape(term)}(?:e?s)?\b", text, re.IGNORECASE) is not None


def leak_counts(texts: Iterable[str], terms: Iterable[str]) -> int:
    return int(sum(term_leaks(text, term) for text, term in zip(texts, terms)))
