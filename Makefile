PYTHON ?= .venv/bin/python
FORCE ?=
ARGS = $(if $(FORCE),--force,)
ACCEPT_RESIDUAL ?=

.PHONY: all week2 week3 week4 setup mdd tree reconcile prune wiki corpus embed matrices matrices-sensitivity sanity h1 ecology attention diagnostics report test posterior posterior-h1 clean-interim clean

all: week2 week3 week4

week2: mdd tree reconcile prune wiki

week3: corpus embed matrices sanity

week4: h1 ecology attention diagnostics report

setup:
	python3 -m venv .venv
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -e '.[dev]'

mdd:
	$(PYTHON) scripts/01_fetch_mdd.py $(ARGS)

tree: mdd
	$(PYTHON) scripts/02_fetch_tree.py $(ARGS)

reconcile: tree
	$(PYTHON) scripts/03_reconcile_names.py $(ARGS)

prune: reconcile
	$(PYTHON) scripts/04_prune_tree.py $(ARGS)

wiki: reconcile
	$(PYTHON) scripts/05_fetch_wikipedia.py $(ARGS)

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
	$(PYTHON) scripts/11_build_attention.py

diagnostics: ecology attention
	$(PYTHON) scripts/12_confound_diagnostics.py

report: h1 diagnostics
	$(PYTHON) scripts/13_week4_report.py

posterior:
	$(PYTHON) scripts/02_fetch_tree.py --posterior

posterior-h1: posterior matrices
	$(PYTHON) scripts/14_posterior_h1.py

test:
	$(PYTHON) -m pytest

clean-interim:
	find data/interim -mindepth 1 -maxdepth 1 ! -name api_cache ! -name wiki -exec rm -rf {} +

clean:
	rm -rf data/interim/* data/processed/*
