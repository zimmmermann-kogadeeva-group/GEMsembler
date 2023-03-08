import copy
import operator
import os
from os.path import join, exists
import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model
import curation
import BiGGnetwork
import conversion
import general
import selection
import structural
import creation
import dill
from copy import deepcopy

if __name__ == '__main__':
    # region Open conversion tables
    fileDir = os.path.dirname(os.path.realpath('__file__'))  # getting directory of the script for pathes to files
    seed_orig = pd.read_csv(join(fileDir, "../Data/seed_to_bigg.tsv"),
                            sep="\t")
    seed_orig_m = seed_orig[seed_orig["type"] == "m"][["seed_ids", "bigg_ids"]]
    seed_orig_m.columns = ["old", "new"]
    seed_orig_r = seed_orig[seed_orig["type"] == "r"][["seed_ids", "bigg_ids"]]
    seed_orig_r.columns = ["old", "new"]
    seed_addit = pd.read_csv(join(fileDir, "../Data/seed_to_bigg_metanetx.tsv"),
                             sep="\t")
    seed_addit_m = seed_addit[seed_addit["type"] == "m"][["seed_ids", "bigg_ids"]]
    seed_addit_m.columns = ["old", "new"]
    seed_addit_r = seed_addit[seed_addit["type"] == "r"][["seed_ids", "bigg_ids"]]
    seed_addit_r.columns = ["old", "new"]
    kegg_bigg = pd.read_csv(join(fileDir, "../Data/kegg_to_bigg_metanetx.tsv"),
                            sep="\t")
    kegg_bigg_m = kegg_bigg[kegg_bigg["type"] == "m"][["kegg_ids", "bigg_ids"]]
    kegg_bigg_m.columns = ["old", "new"]
    kegg_bigg_r = kegg_bigg[kegg_bigg["type"] == "r"][["kegg_ids", "bigg_ids"]]
    kegg_bigg_r.columns = ["old", "new"]
    old_new_bigg = pd.read_csv(join(fileDir, "../Data/old_to_new_bigg.tsv"),
                               sep="\t")
    old_new_bigg_m = old_new_bigg[old_new_bigg["type"] == "m"][["old_bigg_ids", "bigg_ids"]]
    old_new_bigg_m.columns = ["old", "new"]
    old_new_bigg_m["new"] = old_new_bigg_m["new"].str[:-2]
    old_new_bigg_r = old_new_bigg[old_new_bigg["type"] == "r"][["old_bigg_ids", "bigg_ids"]]
    old_new_bigg_r.columns = ["old", "new"]
    old_new_bigg_r["new"] = old_new_bigg_r[
        "new"]  # .str.replace("_[cep]", "", regex = True) - not sure, whether it  is reasnable
    bigg_all_m = pd.read_csv(join(fileDir, "../Data/bigg_models_metabolites.txt"), sep="\t")
    bigg_m = list(set(bigg_all_m["universal_bigg_id"]))
    bigg_all_r = pd.read_csv(join(fileDir, "../Data/bigg_models_reactions.txt"), sep="\t")
    bigg_all_r["universal_bigg_id"] = bigg_all_r["bigg_id"]
    bigg_r = list(set(bigg_all_r["universal_bigg_id"]))

    # endregion
    bigg_network_dill_file = join(fileDir, "../Scripts/bigg_network.pkl")
    if exists(bigg_network_dill_file):
        bigg_db_network = dill.load(open(bigg_network_dill_file, "rb"))
    else:
        bigg_db_network = BiGGnetwork.getBiGGnetwork(bigg_all_r)
        dill.dump(bigg_db_network, open(bigg_network_dill_file, "wb"))

    # region Open models
    model_type_list = ["carveme", "gapseq", "modelseed", "agora"]
    models_same_db = {"modelseed": ["gapseq", "modelseed"]}
    models_to_curate = ["modelseed"]
    models_to_convert = model_type_list[1:]
    models_NOTto_convert = model_type_list[:1]
    models_wo_periplasmic = ["modelseed", "agora"]

    fileDir = os.path.dirname(os.path.realpath('__file__'))  # getting directory of the script for paths to files
    name_carv_hom = join(fileDir, "../Data/BU_carveme_hom.xml")
    name_gapseq = join(fileDir, "../Data/BU_gapseq.xml")
    name_modelseed = join(fileDir, "../Data/BU_modelSEED.sbml")
    name_agora = join(fileDir, "../Data/BU_agora.xml")
    names = {"carveme": name_carv_hom, "gapseq": name_gapseq, "modelseed": name_modelseed, "agora": name_agora}
    all_models_dill_file = join(fileDir, "../Scripts/all_models.pkl")
    if exists(all_models_dill_file):
        all_models = dill.load(open(all_models_dill_file, "rb"))
    else:
        all_models = {}
        for typ in model_type_list:
            nam = names.get(typ)
            model = read_sbml_model(nam)
            all_models.update({typ: model})
        dill.dump(all_models, open(all_models_dill_file, "wb"))

    # endregion

    # region Perform curation
    curated_models = copy.deepcopy(all_models)
    for cur in models_to_curate:
        curated_models[cur] = curation.removeBtypeExchange(curated_models[cur])

    duplicated_reactions = curation.checkDuplicatedReactions(model_type_list, curated_models)
    # endregion

    # region Perform conversion
    GapseqConv = conversion.ConversionForGapseq(seed_orig_m, seed_orig_r, seed_addit_m, seed_addit_r, bigg_m, bigg_r)
    ModelseedConv = conversion.ConversionForModelseed(seed_orig_m, seed_orig_r, seed_addit_m, seed_addit_r, bigg_m,
                                                      bigg_r)
    AgoraConv = conversion.ConversionForAgora(old_new_bigg_m, old_new_bigg_r, kegg_bigg_m, kegg_bigg_r, bigg_m, bigg_r)
    CarvemeConv = conversion.ConversionForCarveMe(old_new_bigg_m, old_new_bigg_r, bigg_m, bigg_r)
    ConversionStrategies = {"gapseq": GapseqConv, "modelseed": ModelseedConv, "agora": AgoraConv,
                            "carveme": CarvemeConv}
    # CarvemeComp = conversion.CompartmentsForCarveme()
    # GapseqComp = conversion.CompartmentsForGapseq()
    # ModelseedComp = conversion.CompartmentsForModelseed()
    # AgoraComp = conversion.CompartmentsForAgora()
    # CompartmentsStrategies = {"carveme": CarvemeComp, "gapseq": GapseqComp, "modelseed": ModelseedComp,
    #                           "agora": AgoraComp}
    allmet_converted_dill_file = join(fileDir, "../Scripts/allmet_converted.pkl")
    allreact_converted_dill_file = join(fileDir, "../Scripts/allreact_converted.pkl")
    if exists(allmet_converted_dill_file) & exists(allreact_converted_dill_file):
        allmet_converted = dill.load(open(allmet_converted_dill_file, "rb"))
        allreact_converted = dill.load(open(allreact_converted_dill_file, "rb"))
    else:
        allmet_converted = conversion.runConversionForALLmodels(models_to_convert, curated_models,
                                                                ConversionStrategies, "metabolites")
        allreact_converted = conversion.runConversionForALLmodels(models_to_convert, curated_models,
                                                                  ConversionStrategies, "reactions")
        dill.dump(allmet_converted, open(allmet_converted_dill_file, "wb"))
        dill.dump(allreact_converted, open(allreact_converted_dill_file, "wb"))
    allmet_checked, allmet_not_pass = conversion.runNoneConversionChecking(models_NOTto_convert, curated_models,
                                                                           ConversionStrategies, "metabolites")
    allreact_checked, allreact_not_pass = conversion.runNoneConversionChecking(models_NOTto_convert, curated_models,
                                                                               ConversionStrategies, "reactions")
    allreact_checked_struct, allreact_not_pass_struct = structural.runStructuralCheck(models_NOTto_convert,
                                                                                      allreact_checked,
                                                                                      allreact_not_pass, curated_models,
                                                                                      bigg_db_network)
    allmet_selected = selection.runSelection(models_to_convert, allmet_converted, "metabolites",
                                             models_same_db)
    allreact_selected = selection.runSelection(models_to_convert, allreact_converted, "reactions",
                                               models_same_db)
    structural_r_info, structural_r_sel = structural.runStructuralConversion(models_to_convert,
                                                                             allmet_selected.get("one_to_one"),
                                                                             allmet_selected, allreact_selected,
                                                                             curated_models,
                                                                             bigg_db_network, models_wo_periplasmic,
                                                                             allmet_selected.get("one_to_many"))
    struct_r_consistent, struct_r_consist, struct_r_not_consist = selection.checkDBConsistency(models_same_db,
                                                                                               structural_r_sel,
                                                                                               "reactions",
                                                                                               write_files=False,
                                                                                               do_stat=False)
    struct_r_uniq, struct_r_not_uniq = selection.checkFromOneFromMany(models_to_convert, struct_r_consistent)
    met_struct = structural.runSuggestionsMet(models_to_convert, structural_r_info, struct_r_uniq, allmet_selected,
                                              models_same_db, curated_models, bigg_db_network)
    struct_final_r_info, struct_final_r_sel = structural.runStructuralConversion(models_to_convert,
                                                                                 met_struct.get("one_one_sugg_met"),
                                                                                 allmet_selected, allreact_selected,
                                                                                 curated_models, bigg_db_network,
                                                                                 models_wo_periplasmic)
    struct_final_r_consistent, struct_final_r_consist, struct_final_r_not_consist = selection.checkDBConsistency(
        models_same_db, struct_final_r_sel,
        "reactions", write_files=False, do_stat=False)
    struct_final_r_uniq, struct_final_r_not_uniq = selection.checkFromOneFromMany(models_to_convert,
                                                                                  struct_final_r_consistent)
    periplasmic_m, periplasmic_r = structural.getSuggestionPeriplasmic(models_wo_periplasmic, struct_final_r_uniq,
                                                                       struct_final_r_info, bigg_db_network,
                                                                       curated_models)
    final_r = deepcopy(struct_final_r_uniq)
    final_r_not_uniq = {}
    for typ in models_to_convert:
        true_dupl = list(
            set(struct_final_r_not_uniq.get(typ).keys()) & set(duplicated_reactions.get(typ)[0]["ID"].tolist()))
        if true_dupl:
            final_r.get(typ).update({td: struct_final_r_not_uniq.get(typ).get(td) for td in true_dupl})
        false_dupl = list(
            set(struct_final_r_not_uniq.get(typ).keys()) - set(duplicated_reactions.get(typ)[0]["ID"].tolist()))
        if false_dupl:
            final_r_not_uniq.update({typ: {fd: struct_final_r_not_uniq.get(typ).get(fd) for fd in false_dupl}})
    final_m_not_sel = selection.runNotSelectedMet(models_to_convert, met_struct, allmet_selected)
    final_r_not_sel = selection.runNotSelectedR(models_to_convert, final_r, struct_final_r_not_consist,
                                                final_r_not_uniq, struct_final_r_info, curated_models)
    final_m = deepcopy(met_struct.get("one_one_sugg_met"))
    additional_p_m = {}
    for typ in models_wo_periplasmic:
        additional_p_m.update({typ: {}})
        for orig_id in periplasmic_m.get(typ).keys():
            if periplasmic_m.get(typ).get(orig_id)[4] == "replace":
                final_m.get(typ)[orig_id] = [final_m.get(typ).get(orig_id)[0], [periplasmic_m.get(typ).get(orig_id)[1]]]
            else:
                additional_p_m.get(typ).update({orig_id: [final_m.get(typ).get(orig_id)[0], [periplasmic_m.get(typ).get(orig_id)[1]]]})

    for typ in models_NOTto_convert:
        final_m.update({typ: allmet_checked.get(typ)})
        final_r.update({typ: allreact_checked_struct.get(typ)})
        m_notsel = general.findKeysByValue(allmet_not_pass.get(typ), "not_found_in_new_and_old_bigg", operator.eq)
        if m_notsel:
            final_m_not_sel.update({typ: {m: [allmet_not_pass.get(typ).get(m)[0], m] for m in m_notsel}})
        else:
            final_m_not_sel.update({typ: {}})
        r_notsel = general.findKeysByValue(allreact_not_pass_struct.get(typ), "not_found_in_new_and_old_bigg",
                                           operator.eq)
        if r_notsel:
            final_r_not_sel.update({typ: {r: [allreact_not_pass_struct.get(typ).get(r)[0], r] for r in r_notsel}})
        else:
            final_r_not_sel.update({typ: {}})
    supermodel = creation.runSupermodelCreation(model_type_list, final_m, final_m_not_sel, final_r, final_r_not_sel,
                                                curated_models, bigg_all_m, bigg_all_r, additional_p_m)
