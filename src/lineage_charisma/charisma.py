from __future__ import annotations

import json
from typing import Any, Mapping, Sequence

import numpy as np

# docs/preregistration_amendment_secondary.md, section B and Change 1. The prompt text below is the registered text;
# do not edit it without a dated Change entry.
DIMENSIONS = ("warmth", "anthropomorphism", "aesthetic", "human_relevance")
SCORE_RANGE = (0, 3)
VALIDATION_SIZE = 50
AGREEMENT_FLOOR = 0.5
MAX_RETRIES = 5

SYSTEM_PROMPT = (
    "You are annotating encyclopedia descriptions of mammal species for a research study. Score only the language of the text you "
    "are given, not what you know about the animal. Some names in the text have been replaced by the token [TAXON]; do not try to "
    "guess them. Respond with JSON only."
)

RUBRIC = """warmth (affective tone): 0 = neutral, technical, no evaluative or emotional words; 1 = one or two mildly positive or emotive words; 2 = several affective or evaluative phrases, a noticeably warm or dramatic tone in places; 3 = affect-laden throughout, written to make the reader feel something about the animal.

anthropomorphism (mental-state language): 0 = no mental states, intentions or personality attributed; 1 = one passing attribution (e.g. "curious", "shy"); 2 = several attributions of feelings, intentions or personality; 3 = the animal is described mainly through human-like mental life, motives, emotions or social roles.

aesthetic (beauty, cuteness, impressiveness): 0 = appearance described only in measurements and neutral colour terms; 1 = one aesthetic or impressive term (e.g. "striking", "large and powerful"); 2 = several, with appearance framed as beautiful, cute or impressive; 3 = aesthetic or impressive framing is central to the description.

human_relevance (danger, usefulness, cultural role): 0 = no mention of humans, human uses, danger or culture; 1 = one brief mention (e.g. hunted for fur, a single cultural reference); 2 = several mentions, or one developed passage, on danger, economic use, conflict or cultural role; 3 = human relations (danger, use, symbolism, conflict) are a major theme of the description."""

USER_TEMPLATE = (
    "Score the description below from 0 to 3 on each of four dimensions, using these definitions.\n\n"
    + RUBRIC
    + '\n\nReturn exactly this JSON object and nothing else, with an integer from 0 to 3 for each key:\n'
    '{"warmth": <0-3>, "anthropomorphism": <0-3>, "aesthetic": <0-3>, "human_relevance": <0-3>}\n\n'
    'Description:\n"""\n{TEXT}\n"""'
)

OUTPUT_SCHEMA = {
    "type": "object",
    "properties": {d: {"type": "integer", "enum": [0, 1, 2, 3]} for d in DIMENSIONS},
    "required": list(DIMENSIONS),
    "additionalProperties": False,
}


class InvalidResponse(ValueError):
    """A model response that does not meet the registered format; it is logged and the request retried."""


def user_message(text: str) -> str:
    return USER_TEMPLATE.replace("{TEXT}", text)


def validate_scores(raw: str | Mapping[str, Any]) -> dict[str, int]:
    """Parse and check one response: a JSON object with exactly the four keys, each an integer from 0 to 3."""
    if isinstance(raw, str):
        try:
            data = json.loads(raw)
        except json.JSONDecodeError as e:
            raise InvalidResponse(f"not valid JSON: {e.msg}") from None
    else:
        data = raw
    if not isinstance(data, dict):
        raise InvalidResponse(f"expected a JSON object, got {type(data).__name__}")
    if set(data) != set(DIMENSIONS):
        raise InvalidResponse(f"keys {sorted(data)} are not exactly {list(DIMENSIONS)}")
    out = {}
    lo, hi = SCORE_RANGE
    for d in DIMENSIONS:
        v = data[d]
        if isinstance(v, bool) or not isinstance(v, int):
            raise InvalidResponse(f"{d} = {v!r} is not an integer")
        if not lo <= v <= hi:
            raise InvalidResponse(f"{d} = {v} is outside {lo}..{hi}")
        out[d] = v
    return out


def validation_sample(order: Sequence[str], seed: int, size: int = VALIDATION_SIZE) -> list[str]:
    """The registered hand-labelling sample: rng(seed).choice(len(order), size, replace=False), indexing the pilot order."""
    idx = np.random.default_rng(seed).choice(len(order), size=size, replace=False)
    return [order[i] for i in idx]


def charisma_distance(scores: np.ndarray) -> np.ndarray:
    """Euclidean distance between species on the retained dimensions, each standardized to mean 0, SD 1 across species.

    `scores` is species x retained dimensions (mean of the two runs). With one dimension this is the absolute difference
    of the standardized scores.
    """
    scores = np.asarray(scores, dtype=np.float64)
    if scores.ndim == 1:
        scores = scores[:, None]
    if scores.shape[1] == 0:
        raise ValueError("no retained dimensions: the charisma-content matrix is undefined")
    sd = scores.std(axis=0)
    if np.any(sd == 0):
        raise ValueError("a retained dimension has no variance across species")
    z = (scores - scores.mean(axis=0)) / sd
    diff = z[:, None, :] - z[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    np.fill_diagonal(dist, 0.0)
    return (dist + dist.T) / 2.0
