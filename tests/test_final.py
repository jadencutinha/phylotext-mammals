from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name: str):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), ROOT / "scripts" / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def verify():
    return load("34_verify_numbers.py")


@pytest.fixture(scope="module")
def figures():
    return load("35_final_figures.py")


def test_quoted_number_matches_at_its_own_precision(verify):
    assert verify.matches("0.033", 0.03270586)
    assert verify.matches("0.0285", 0.0285)
    assert not verify.matches("0.442", 0.44147)          # the recorded erratum
    assert verify.matches("0.441", 0.44147)
    assert verify.matches("−0.028", -0.02821)            # typographic minus
    assert verify.matches("+0.056", 0.05633)
    assert not verify.matches("0.034", -0.034)


def test_every_erratum_names_a_checked_quantity(verify):
    checked = {(desc, file) for desc, file, _, _ in verify.claims(ROOT / "reports" / "tables")}
    diagnostic = {("taxonomy mask, different-order mean text distance", f) for f in
                  ("reports/week6_report.md", "docs/preregistration_amendment_secondary.md")}
    assert set(verify.ERRATA) <= checked | diagnostic


def test_p_text_floor_and_precision(figures):
    assert figures.p_text(0.0001, 9999) == "p ≤ 0.0001"
    assert figures.p_text(0.001, 999) == "p ≤ 0.0010"
    assert figures.p_text(0.0285, 9999) == "p = 0.0285"


def test_charisma_figure_is_skipped_without_results(figures, tmp_path):
    (tmp_path / "secondary").mkdir()
    assert figures.fig_charisma(tmp_path, tmp_path) is False
    assert not list(tmp_path.glob("fig10*"))
