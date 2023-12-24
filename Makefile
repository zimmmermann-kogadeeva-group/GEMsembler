
SHELL = /bin/bash
ACTIVATE = source .venv/bin/activate 

.venv:
	python3 -m venv $@

.PHONY: install_packages
install_packages: .venv requirements.txt
	${ACTIVATE} && pip install -r requirements.txt

requirements.txt:
	${ACTIVATE} && pip freeze > $@

.PHONY: build
build:
	${ACTIVATE} && python3 -m build

.PHONY: upload
upload: build
	${ACTIVATE} && python3 -m twine upload -r pypi -u __token__ dist/*

.PHONY: clean
clean:
	rm -rf example/*.pkl *.egg-info/ **/__pycache__/ build/ dist/ report.xml

PAPER_URL = https://raw.githubusercontent.com/SystemsBioinformatics/pub-data/master/reconstruction-tools-assessment/pipeline/reconstructions/
VMH_URL = https://www.vmh.life/files/reconstructions/AGORA/
LP_FILES = CA1.xml.gz \
		   CA2.xml.gz \
		   MS2.sbml.gz \
		   protein_fasta.faa.gz \
		   WCFS1_agora.xml.gz \
		   WCFS1.fasta.gz

.PHONY: lp_data
lp_data: $(addprefix example/LP_,${LP_FILES})

example/LP_CA1.xml.gz:
	wget -qO - ${PAPER_URL}/CA/LPL/CA1.xml | gzip -c > $@

example/LP_CA2.xml.gz:
	wget -qO - ${PAPER_URL}/CA/LPL/CA2.xml | gzip -c > $@

example/LP_MS2.sbml.gz:
	wget -qO - ${PAPER_URL}/MS/LPL/MS2.sbml | gzip -c > $@

example/LP_protein_fasta.faa.gz:
	wget -qO - ${PAPER_URL}/inputs/LPL/protein_fasta.faa | gzip -c > $@

example/LP_iLP728.xml.gz:
	wget -qO - ${PAPER_URL}/manually_curated_models/LPL/iLP728.xml | gzip -c > $@

example/LP_WCFS1_agora.xml.gz:
	wget -qO - ${VMH_URL}/1.03/reconstructions/sbml/Lactobacillus_plantarum_WCFS1.xml | gzip -c > $@ 

example/LP_WCFS1.fasta.gz:
	wget -q - ${VMH_URL}/genomes/AGORA-Genomes.zip -O agora_genomes.zip ; \
	unzip -p agora_genomes.zip Lactobacillus_plantarum_WCFS1.fasta | gzip -c > $@ && \
	rm agora_genomes.zip
