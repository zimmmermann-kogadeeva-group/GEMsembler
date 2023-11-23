from .gathering import GatheredModels


if __name__ == "__main__":
    pass
"""       

    # Genes convertions
    input_genomes_names = {
        "agora": "Bacteroides_uniformis_ATCC_8492.fasta",
        "carveme": "GCF_000154205.1_ASM15420v1_protein.faa",
        "gapseq": "GCF_000154205.1_ASM15420v1_genomic.fna",
        "modelseed": "GCF_000154205.1_ASM15420v1_protein.faa",
    }
    gapseq_genes_names = {"gapseq": "Gapseq_model_genes.fasta"}
    ncbi_cds_name = "GCF_000154205.1_ASM15420v1_cds_from_genomic.fna"
    ncbi_protein_name = "GCF_000154205.1_ASM15420v1_protein.faa"
    feature_table_name = "GCF_000154205.1_ASM15420v1_feature_table.txt"
    out_nt_fasta = "NCBI_cds_old_locus_tag.fasta"
    out_aa_fasta = "NCBI_proteins_old_locus_tag.fasta"
    genes.runGenesConversion(
        model_type_list,
        curated_models,
        input_genomes_names,
        gapseq_genes_names,
        ncbi_cds_name,
        ncbi_protein_name,
        feature_table_name,
        out_nt_fasta,
        out_aa_fasta,
        do_old_locus_tag=True,
    )
    # creating supermodel
    supermodel = creation.runSupermodelCreation(
        model_type_list,
        final_m,
        final_m_not_sel,
        final_r,
        final_r_not_sel,
        curated_models,
        bigg_all_m,
        bigg_all_r,
        allmet_converted,
        allreact_converted,
        additional_p_m,
        periplasmic_r,
    )
    # getting core and different types of intersections in supermodel
    comparison.runComparison(supermodel)
    core_model = anticreation.getModelOfInterest(
        supermodel, "core4", name="BU_core_model.xml"
    )
    union_model = anticreation.getModelOfInterest(
        supermodel, "core1", name="BU_union_model.xml"
    )
    Yes_a_No_cmg_model = anticreation.getModelOfInterest(
        supermodel, "Yes_a_No_cgm", name="BU_Yes_a_No_cmg_model.xml"
    )
    carveme_model = anticreation.getModelOfInterest(
        supermodel, "carveme", name="BU_carveme_out_model.xml"
    ) 
        
        
"""
