from .gathering import GatheredModels


if __name__ == "__main__":
    pass
"""       
    # checking ids from models that supposed to be with BiGG but still may be old (or potentialy wrond but that i don't remember)
    # allmet_checked, allmet_not_pass = conversion.runNoneConversionChecking(
    #     models_NOTto_convert, curated_models, ConversionStrategies, "metabolites"
    # )
    # allreact_checked, allreact_not_pass = conversion.runNoneConversionChecking(
    #     models_NOTto_convert, curated_models, ConversionStrategies, "reactions"
    # )

    # checking reaction equations for models with BiGG
    allreact_checked_struct, allreact_not_pass_struct = structural.runStructuralCheck(
        models_NOTto_convert,
        allreact_checked,
        allreact_not_pass,
        curated_models,
        bigg_db_network,
    )

    # converting reactions via reactions equations from BiGG database via metabolites that were converted 1-1 or 1-n and selecting ones that were converted with equations and were converted uniquely
    structural_r_info, structural_r_sel = structural.runStructuralConversion(
        models_to_convert,
        allmet_selected.get("one_to_one"),
        allmet_selected,
        allreact_selected,
        curated_models,
        bigg_db_network,
        models_wo_periplasmic,
        allmet_selected.get("one_to_many"),
    )
    # checking consistency in reactions uniquely converted via equation for different models with same database
    (
        struct_r_consistent,
        struct_r_consist,
        struct_r_not_consist,
    ) = selection.checkDBConsistency(
        models_same_db, structural_r_sel, "reactions", write_files=False, do_stat=False
    )
    # selecting which consistent structural reactions are converted 1-1 and n-1
    struct_r_uniq, struct_r_not_uniq = selection.checkFromOneFromMany(
        models_to_convert, struct_r_consistent
    )
    # getting suggestions for 1-n metabolites from first run of structural conversion and tring to structural with n-1 metabolites. Get 1-1 metabolites plus suggestions
    met_struct = structural.runSuggestionsMet(
        models_to_convert,
        structural_r_info,
        struct_r_uniq,
        allmet_selected,
        models_same_db,
        curated_models,
        bigg_db_network,
    )
    # running structural reaction conversion second time with update metabolites (1-1 plus suggestions)
    struct_final_r_info, struct_final_r_sel = structural.runStructuralConversion(
        models_to_convert,
        met_struct.get("one_one_sugg_met"),
        allmet_selected,
        allreact_selected,
        curated_models,
        bigg_db_network,
        models_wo_periplasmic,
    )
    # checking consistency in 2d run of reactions uniquely converted via equation for different models with same database
    (
        struct_final_r_consistent,
        struct_final_r_consist,
        struct_final_r_not_consist,
    ) = selection.checkDBConsistency(
        models_same_db,
        struct_final_r_sel,
        "reactions",
        write_files=False,
        do_stat=False,
    )
    # selecting which consistent structural reactions from run 2 are converted 1-1 and n-1
    struct_final_r_uniq, struct_final_r_not_uniq = selection.checkFromOneFromMany(
        models_to_convert, struct_final_r_consistent
    )
    # getting metabolites and reactions that became preiplasmic for models without periplasmic compartment originally
    periplasmic_m, periplasmic_r = structural.getSuggestionPeriplasmic(
        models_wo_periplasmic,
        struct_final_r_uniq,
        struct_final_r_info,
        bigg_db_network,
        curated_models,
    )
    # getting in standart format final metabolites and reactions for supermodel creation (converted)
    final_r = deepcopy(struct_final_r_uniq)
    final_r_not_uniq = {}
    for (
        typ
    ) in (
        models_to_convert
    ):  # adding n-1 (duplicated) reactions structural reactions if they are indeed originaly duplicated in models
        true_dupl = list(
            set(struct_final_r_not_uniq.get(typ).keys())
            & set(duplicated_reactions.get(typ)[0]["ID"].tolist())
        )
        if true_dupl:
            final_r.get(typ).update(
                {td: struct_final_r_not_uniq.get(typ).get(td) for td in true_dupl}
            )
        false_dupl = list(
            set(struct_final_r_not_uniq.get(typ).keys())
            - set(duplicated_reactions.get(typ)[0]["ID"].tolist())
        )
        if false_dupl:
            final_r_not_uniq.update(
                {
                    typ: {
                        fd: struct_final_r_not_uniq.get(typ).get(fd)
                        for fd in false_dupl
                    }
                }
            )
    # getting in standard format not selected for not converted in supermodel
    final_m_not_sel = selection.runNotSelectedMet(
        models_to_convert, met_struct, allmet_selected
    )
    final_r_not_sel = selection.runNotSelectedR(
        models_to_convert,
        final_r,
        struct_final_r_not_consist,
        final_r_not_uniq,
        struct_final_r_info,
        curated_models,
    )
    final_m = deepcopy(met_struct.get("one_one_sugg_met"))
    additional_p_m = {}
    # dealing with periplasmic metabolites: replacing if final and creating additional periplasmic metabolites if original metabolite works with boths compartments
    for typ in models_wo_periplasmic:
        additional_p_m.update({typ: {}})
        for orig_id in periplasmic_m.get(typ).keys():
            if periplasmic_m.get(typ).get(orig_id)[4] == "replace":
                final_m.get(typ)[orig_id] = [
                    final_m.get(typ).get(orig_id)[0],
                    [periplasmic_m.get(typ).get(orig_id)[1]],
                ]
            else:
                additional_p_m.get(typ).update(
                    {
                        orig_id: [
                            final_m.get(typ).get(orig_id)[0],
                            [periplasmic_m.get(typ).get(orig_id)[1]],
                        ]
                    }
                )
    # combining final metabolites and reactions with ones coming from BiGG models (wo conversion)
    for typ in models_NOTto_convert:
        final_m.update({typ: allmet_checked.get(typ)})
        final_r.update({typ: allreact_checked_struct.get(typ)})
        m_notsel = general.findKeysByValue(
            allmet_not_pass.get(typ), "not_found_in_new_and_old_bigg", operator.eq
        )
        if m_notsel:
            final_m_not_sel.update(
                {typ: {m: [allmet_not_pass.get(typ).get(m)[0], m] for m in m_notsel}}
            )
        else:
            final_m_not_sel.update({typ: {}})
        r_notsel = [
            k
            for k, v in allreact_not_pass_struct.get(typ).items()
            if v[1] == "not_found_in_new_and_old_bigg"
        ]
        if r_notsel:
            final_r_not_sel.update(
                {
                    typ: {
                        r: [allreact_not_pass_struct.get(typ).get(r)[0], [r]]
                        for r in r_notsel
                    }
                }
            )
        else:
            final_r_not_sel.update({typ: {}})
        for kg, vg in allreact_not_pass_struct.get(typ).items():
            if vg[1] == "Growth_reaction":
                final_r_not_sel.get(typ).update(
                    {kg: [allreact_not_pass_struct.get(typ).get(kg)[0], ["Biomass"]]}
                )
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
