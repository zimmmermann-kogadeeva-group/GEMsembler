
all: \
	genome.fna \
	genome_agora.fna \
	agora_model.xml
	# LP/genome.faa \

# LP

# This protein fasta file does not match the one that was downloaded couple of
# months ago
LP_GENOME_FAA_URL='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/203/855/GCF_000203855.3_ASM20385v3/GCF_000203855.3_ASM20385v3_protein.faa.gz'
LP/LP_genome.faa.gz:
	curl ${LP_GENOME_FAA_URL} -o $@

LP_GENOME_FNA_URL='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/203/855/GCF_000203855.3_ASM20385v3/GCF_000203855.3_ASM20385v3_genomic.fna.gz'
LP/LP_genome.fna.gz:
	curl ${LP_GENOME_FNA_URL} -o $@

LP_GENOME_AGORA_FNA_URL='https://www.vmh.life/files/reconstructions/AGORA/genomes/AGORA-Genomes.zip'
LP/LP_genome_agora.fna.gz:
	wget ${LP_GENOME_AGORA_FNA_URL} && \
	unzip -p AGORA-Genomes.zip Lactobacillus_plantarum_WCFS1.fasta | gzip > $@ && \
	rm AGORA-Genomes.zip

LP_AGORA_URL='https://www.vmh.life/files/reconstructions/AGORA2/version2.01/sbml_files/individual_reconstructions/Lactobacillus_plantarum_WCFS1.xml'
LP/LP_agora_model.xml:
	curl ${LP_AGORA_URL} -o $@

LP/LP_agora_model.xml.gz: agora_model.xml
	gzip $^
