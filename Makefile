
SHELL = /bin/bash
ACTIVATE = source .venv/bin/activate 

.PHONY: install_packages build upload check isort tags clean data

.venv:
	python3 -m venv $@

install_packages: .venv
	${ACTIVATE} && pip install -e .

requirements.txt:
	${ACTIVATE} && pip freeze > $@

build:
	${ACTIVATE} && python3 -m build

upload: build
	${ACTIVATE} && python3 -m twine upload -r pypi -u __token__ dist/*

check:
	ruff check src tests

isort:
	isort src tests

tags:
	ctags-universal --recurse src tests

clean:
	rm -rf src/*.egg-info/ **/__pycache__/ build/ dist/ report.xml

data:
	${MAKE} -C src/gemsembler/data/
