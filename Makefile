
SHELL = /bin/bash
ACTIVATE = source .venv/bin/activate 

.PHONY: install_packages build upload check isort tags clean

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

PAPER_URL = https://raw.githubusercontent.com/SystemsBioinformatics/pub-data/master/reconstruction-tools-assessment/pipeline/reconstructions/
VMH_URL = https://www.vmh.life/files/reconstructions/AGORA/
LP_FILES = CA1.xml.gz \
		   CA2.xml.gz \
		   MS2.sbml.gz \
		   protein_fasta.faa.gz \
		   WCFS1_agora.xml.gz \
		   WCFS1.fasta.gz

.PHONY: lp_data
lp_data: $(addprefix src/gemsembler/data/LP/LP_,${LP_FILES})

src/gemsembler/data/LP/LP_CA1.xml.gz:
	mkdir -p $(dir $@); wget -qO - ${PAPER_URL}/CA/LPL/CA1.xml | gzip -c > $@

src/gemsembler/data/LP/LP_CA2.xml.gz:
	mkdir -p $(dir $@); wget -qO - ${PAPER_URL}/CA/LPL/CA2.xml | gzip -c > $@

src/gemsembler/data/LP/LP_MS2.sbml.gz:
	mkdir -p $(dir $@); wget -qO - ${PAPER_URL}/MS/LPL/MS2.sbml | gzip -c > $@

src/gemsembler/data/LP/LP_protein_fasta.faa.gz:
	mkdir -p $(dir $@); wget -qO - ${PAPER_URL}/inputs/LPL/protein_fasta.faa | gzip -c > $@

src/gemsembler/data/LP/LP_iLP728.xml.gz:
	mkdir -p $(dir $@); wget -qO - ${PAPER_URL}/manually_curated_models/LPL/iLP728.xml | gzip -c > $@

src/gemsembler/data/LP/LP_WCFS1_agora.xml.gz:
	mkdir -p $(dir $@); wget -qO - ${VMH_URL}/1.03/reconstructions/sbml/Lactobacillus_plantarum_WCFS1.xml | gzip -c > $@ 

src/gemsembler/data/LP/LP_WCFS1.fasta.gz:
	mkdir -p $(dir $@); \
	wget -q - ${VMH_URL}/genomes/AGORA-Genomes.zip -O agora_genomes.zip ; \
	unzip -p agora_genomes.zip Lactobacillus_plantarum_WCFS1.fasta | gzip -c > $@ && \
	rm agora_genomes.zip
