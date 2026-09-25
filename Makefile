PYTHON ?= .venv/bin/python
FORCE ?=
ARGS = $(if $(FORCE),--force,)

.PHONY: all setup mdd tree reconcile prune wiki test posterior clean-interim clean

all: mdd tree reconcile prune wiki

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

posterior:
	$(PYTHON) scripts/02_fetch_tree.py --posterior

test:
	$(PYTHON) -m pytest

clean-interim:
	find data/interim -mindepth 1 -maxdepth 1 ! -name api_cache ! -name wiki -exec rm -rf {} +

clean:
	rm -rf data/interim/* data/processed/*
