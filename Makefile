PYTHON ?= .venv/bin/python
FORCE ?=
CLADE ?=
ifneq ($(CLADE),)
export LC_CLADE := $(CLADE)
endif
ARGS = $(if $(FORCE),--force,)
OFFLINE ?=
NET = $(if $(OFFLINE),--offline,)
FULL = --clade mammalia
ACCEPT_RESIDUAL ?=
H2_ARGS ?=

.PHONY: setup-locked all pilot fullscale secondary final charisma zoo lock week2 week3 week4 week5 setup mdd tree reconcile prune wiki corpus embed matrices matrices-sensitivity sanity h1 ecology attention diagnostics report h2 exploratory report5 test posterior posterior-h1 clean-interim clean

# Every result, in the order it was produced. With OFFLINE=1 the API stages read their caches only and fail on a miss.
# Needs the posterior archive (make posterior) for the pilot H2 and full-scale posterior trees. Analysis B (make charisma) needs
# ANTHROPIC_API_KEY and committed hand labels, and zoo signage (make zoo) needs data/raw/zoo_signage/signage.csv, so neither is in `all`;
# the final stage reports them as pending or not run when their tables are absent.
all: pilot fullscale secondary final

pilot: week2 week3 week4 posterior-h1 week5

week2: mdd tree reconcile prune wiki

week3: corpus embed matrices sanity

week4: h1 ecology attention diagnostics report

# h2 reads the posterior archive (make posterior, about 1.2 GB); H2_ARGS=--skip-posterior leaves those trees out
week5: h2 exploratory report5

setup:
	python3 -m venv .venv
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -e '.[dev]'

# the exact environment that produced the results (Python 3.12.10, macOS arm64; see requirements.lock)
setup-locked:
	python3.12 -m venv .venv
	$(PYTHON) -m pip install -r requirements.lock
	$(PYTHON) -m pip install --no-deps -e .

mdd:
	$(PYTHON) scripts/01_fetch_mdd.py $(ARGS)

tree: mdd
	$(PYTHON) scripts/02_fetch_tree.py $(ARGS)

reconcile: tree
	$(PYTHON) scripts/03_reconcile_names.py $(ARGS) $(NET)

prune: reconcile
	$(PYTHON) scripts/04_prune_tree.py $(ARGS)

wiki: reconcile
	$(PYTHON) scripts/05_fetch_wikipedia.py $(ARGS) $(NET)

corpus: wiki
	$(PYTHON) scripts/05b_build_corpus.py $(ARGS)

embed: corpus
	$(PYTHON) scripts/06_embed.py $(ARGS)

matrices: embed prune
	$(PYTHON) scripts/07_build_matrices.py $(ARGS)

matrices-sensitivity: embed prune
	$(PYTHON) scripts/07_build_matrices.py --sensitivity $(ARGS)

sanity: matrices
	$(PYTHON) scripts/08_sanity_report.py

h1: matrices
	$(PYTHON) scripts/09_run_h1.py

ecology: matrices
	$(PYTHON) scripts/10_build_ecology.py $(ARGS) $(if $(ACCEPT_RESIDUAL),--accept-residual,)

attention: matrices
	$(PYTHON) scripts/11_build_attention.py $(NET)

diagnostics: ecology attention
	$(PYTHON) scripts/12_confound_diagnostics.py

report: h1 diagnostics
	$(PYTHON) scripts/13_week4_report.py

posterior:
	$(PYTHON) scripts/02_fetch_tree.py --posterior

posterior-h1: posterior matrices
	$(PYTHON) scripts/14_posterior_h1.py

h2: ecology attention
	$(PYTHON) scripts/15_run_h2.py $(H2_ARGS)

exploratory: ecology attention
	$(PYTHON) scripts/16_exploratory.py

report5: h2 exploratory
	$(PYTHON) scripts/17_week5_report.py

# Full scale, all mammals (amendment 3). The frozen species list is never overwritten (19 and 07 stop if it would change).
fullscale: posterior
	$(PYTHON) scripts/01_fetch_mdd.py $(FULL) $(ARGS)
	$(PYTHON) scripts/02_fetch_tree.py $(FULL) $(ARGS)
	$(PYTHON) scripts/03_reconcile_names.py $(FULL) $(ARGS) $(NET)
	$(PYTHON) scripts/03b_residual_rules.py $(FULL)
	$(PYTHON) scripts/04_prune_tree.py $(FULL) $(ARGS)
	$(PYTHON) scripts/05_fetch_wikipedia.py $(FULL) $(ARGS) $(NET)
	$(PYTHON) scripts/05b_build_corpus.py $(FULL) $(ARGS)
	$(PYTHON) scripts/19_freeze_species.py $(FULL)
	$(PYTHON) scripts/06_embed.py $(FULL) $(ARGS)
	$(PYTHON) scripts/07_build_matrices.py $(FULL) $(ARGS)
	$(PYTHON) scripts/10_build_ecology.py $(FULL) $(ARGS)
	$(PYTHON) scripts/11_build_attention.py $(FULL) $(NET)
	$(PYTHON) scripts/20_coverage_by_order.py $(FULL)
	$(PYTHON) scripts/21_replication.py
	$(PYTHON) scripts/22_fullscale_tests.py $(FULL)
	$(PYTHON) scripts/23_posterior_fullscale.py $(FULL) --trees-only
	$(PYTHON) scripts/23_posterior_fullscale.py $(FULL)
	$(PYTHON) scripts/24_exploratory_fullscale.py $(FULL)
	$(PYTHON) scripts/25_week6_figures.py $(FULL)

# Amendment 4: post-hoc D, then secondary A (ADW pages come from data/raw/adw/ once fetched; OFFLINE=1 never requests)
secondary:
	$(PYTHON) scripts/26_length_matched.py $(FULL)
	$(PYTHON) scripts/29_fetch_adw.py $(NET)
	$(PYTHON) scripts/30_register_comparison.py

# Secondary B: refuses to call the model until the 50 hand labels are committed; reads ANTHROPIC_API_KEY; responses are cached.
# The labelling sheet itself was written once by scripts/27_annotation_sheet.py, which never overwrites a sheet that exists.
charisma:
	$(PYTHON) scripts/31_annotate.py --run 1
	$(PYTHON) scripts/31_annotate.py --run 2
	$(PYTHON) scripts/32_charisma_analysis.py

# Secondary C: needs the transcribed signs
zoo:
	$(PYTHON) scripts/28_zoo_signage.py

# Week 8 lock: results index, number check against the tables, publication figures
final: lock
	$(PYTHON) scripts/35_final_figures.py

lock:
	$(PYTHON) scripts/33_results_index.py
	$(PYTHON) scripts/34_verify_numbers.py

test:
	$(PYTHON) -m pytest

clean-interim:
	find data/interim -mindepth 1 -maxdepth 1 ! -name api_cache ! -name wiki -exec rm -rf {} +

clean:
	rm -rf data/interim/* data/processed/*
