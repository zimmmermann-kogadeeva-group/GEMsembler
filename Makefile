
SHELL = /bin/bash
ACTIVATE = source .venv/bin/activate 

.venv:
	python3 -m venv $@

.PHONY: install_packages
install_packages: .venv requirements.txt
	${ACTIVATE} && pip install -r requirements.txt

requirements.txt:
	${ACTIVATE} && pip freeze > $@

.PHONY: clean
clean:
	rm -rf example/*.pkl *.egg-info/ **/__pycache__/ build/ report.xml
