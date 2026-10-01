PYTHON ?= .venv/bin/python
FORCE ?=
ARGS = $(if $(FORCE),--force,)

.PHONY: all week2 week3 setup mdd tree reconcile prune wiki corpus embed matrices matrices-sensitivity sanity test posterior clean-interim clean

all: week2 week3

week2: mdd tree reconcile prune wiki

week3: corpus embed matrices sanity

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

posterior:
	$(PYTHON) scripts/02_fetch_tree.py --posterior

test:
	$(PYTHON) -m pytest

clean-interim:
	find data/interim -mindepth 1 -maxdepth 1 ! -name api_cache ! -name wiki -exec rm -rf {} +

clean:
	rm -rf data/interim/* data/processed/*
