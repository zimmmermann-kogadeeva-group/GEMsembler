
PAPER_URL = https://raw.githubusercontent.com/SystemsBioinformatics/pub-data/master/reconstruction-tools-assessment/pipeline/reconstructions/
VMH_URL = https://www.vmh.life/files/reconstructions/AGORA/
LP_FILES = CA1.xml.gz \
	   CA2.xml.gz \
	   MS2.sbml.gz \
	   protein_fasta.faa.gz \
	   WCFS1_agora.xml.gz \
	   WCFS1.fasta.gz

all: $(addprefix LP/LP_,${LP_FILES})

LP:
	mkdir -p $@

# model from paper - carveme with universal template
LP/LP_CA1.xml.gz: LP
	wget -qO - ${PAPER_URL}/CA/LPL/CA1.xml | gzip -c > $@

# model from paper - carveme with universal gram-positive template
LP/LP_CA2.xml.gz: LP
	wget -qO - ${PAPER_URL}/CA/LPL/CA2.xml | gzip -c > $@

# model from paper (built with modelseed with gram-positive template)
LP/LP_MS2.sbml.gz: LP
	wget -qO - ${PAPER_URL}/MS/LPL/MS2.sbml | gzip -c > $@

# genome from paper
LP/LP_protein_fasta.faa.gz: LP
	wget -qO - ${PAPER_URL}/inputs/LPL/protein_fasta.faa | gzip -c > $@

# curated model from paper (similar to carveme)
LP/LP_iLP728.xml.gz: LP
	wget -qO - ${PAPER_URL}/manually_curated_models/LPL/iLP728.xml | gzip -c > $@

# agora model
LP/LP_WCFS1_agora.xml.gz: LP
	wget -qO - ${VMH_URL}/1.03/reconstructions/sbml/Lactobacillus_plantarum_WCFS1.xml | gzip -c > $@ 

# genome for agora model
LP/LP_WCFS1.fasta.gz: LP
	wget -q - ${VMH_URL}/genomes/AGORA-Genomes.zip -O agora_genomes.zip ; \
	unzip -p agora_genomes.zip Lactobacillus_plantarum_WCFS1.fasta | gzip -c > $@ && \
	rm agora_genomes.zip
