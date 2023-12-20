# GEMsembler

GEMsembler tool for assembling and comparing several types of Genome-Scale Metabolic
Models. 

**THIS IS A BETA VERSION! BUGS CAN BE EXPECTED**

### Instalation

Instalation from git via pip
```
pip install git+ssh://git@git.embl.de/grp-zimmermann-kogadeeva/GEMsembler.git
```

You also have to install BLAST in advance

### Usage

Input models have to be COBRApy readable files. And models need to be particular type. Curently models made by CarveMe (carveme), ModelSEED (modelseed), gapseq (gapseq) and models downloaded from AGORA VMH database (agora) are supported. Custom type is comming soon. Genomes, from which the models are built will allow to convert and assemble genes as well.

```
input_data = {
    "curated_LP": {"path_to_model": "./iLP728_revision_data_met_C_c.xml", "model_type": "carveme", "path_to_genome": "./protein_fasta.faa"},
    "cauniv_LP": {"path_to_model": "./CA1.xml", "model_type": "carveme", "path_to_genome": "./protein_fasta.faa"},
    "cagram_LP": {"path_to_model": "./CA2.xml", "model_type": "carveme", "path_to_genome": "./protein_fasta.faa"},
    "msgram_LP": {"path_to_model": "./MS2.sbml", "model_type": "modelseed", "path_to_genome": "./protein_fasta.faa"},
    "agora_LP": {"path_to_model": "./Lactobacillus_plantarum_WCFS1_agora.xml", "model_type": "agora", "path_to_genome": "./Lactobacillus_plantarum_WCFS1.fasta"}
}
```

Fisrt stage is the creation of gathered models, a class, that performs conversion and contains results of all stages.

```
gathered = GatheredModels()
for i, v in input_data.items():
    gathered.add_model(i, **v)
gathered.run()

```

Second stage is actual assembly of supermodel from the in formation in gathered models. User has to provide output folder. And for gene conversion user hast provide either final genes in fasta. Then all gene will be converted to ids in these files. Or if user provides NCBI assembly ID for his organism of interest, corresponding genome will be downloaded autpmaticly and all genes will be converted to the locus tags of the organis.
```
supermodel_lp = gathered.assemble_supermodel("./gemsembler_output/", assembly_id = "GCF_000203855.3")
```

After supermodel is assembled different comparison methods can be run
```
supermodel_lp.at_least_in(2)
```

And results of comparison can be extracted as typical COBRApy models
```
core2 = get_model_of_interest(supermodel_lp, "core2", "./gemsembler_output/LP_core2_output_model.xml")
```
