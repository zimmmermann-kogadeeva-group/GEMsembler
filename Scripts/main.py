import copy
import os
from os.path import join, exists
import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model
import curation
import BiGGnetwork
import conversion
import selection
import structural
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
    models_to_curate = ["modelseed"]
    curated_models = copy.deepcopy(all_models)
    for cur in models_to_curate:
        curated_models[cur] = curation.removeBtypeExchange(curated_models[cur])
    # endregion

    # region Perform conversion
    GapseqConv = conversion.ConversionForGapseq(seed_orig_m, seed_orig_r, seed_addit_m, seed_addit_r, bigg_m, bigg_r)
    ModelseedConv = conversion.ConversionForModelseed(seed_orig_m, seed_orig_r, seed_addit_m, seed_addit_r, bigg_m,
                                                      bigg_r)
    AgoraConv = conversion.ConversionForAgora(old_new_bigg_m, old_new_bigg_r, kegg_bigg_m, kegg_bigg_r, bigg_m, bigg_r)
    ConversionStrategies = {"gapseq": GapseqConv, "modelseed": ModelseedConv, "agora": AgoraConv}
    CarvemeComp = conversion.CompartmentsForCarveme()
    GapseqComp = conversion.CompartmentsForGapseq()
    ModelseedComp = conversion.CompartmentsForModelseed()
    AgoraComp = conversion.CompartmentsForAgora()
    CompartmentsStrategies = {"carveme": CarvemeComp, "gapseq": GapseqComp, "modelseed": ModelseedComp,
                              "agora": AgoraComp}
    models_to_convert = model_type_list[1:]
    allmet_converted_dill_file = join(fileDir, "../Scripts/allmet_converted.pkl")
    allreact_converted_dill_file = join(fileDir, "../Scripts/allreact_converted.pkl")
    if exists(allmet_converted_dill_file) & exists(allreact_converted_dill_file):
        allmet_converted = dill.load(open(allmet_converted_dill_file, "rb"))
        allreact_converted = dill.load(open(allreact_converted_dill_file, "rb"))
    else:
        allmet_converted = conversion.runConversionForALLmodels(models_to_convert, curated_models,
                                                                CompartmentsStrategies,
                                                                ConversionStrategies, "metabolites")
        allreact_converted = conversion.runConversionForALLmodels(models_to_convert, curated_models,
                                                                  CompartmentsStrategies,
                                                                  ConversionStrategies, "reactions")
        dill.dump(allmet_converted, open(allmet_converted_dill_file, "wb"))
        dill.dump(allreact_converted, open(allreact_converted_dill_file, "wb"))
    allmet_selected = selection.runSelection(models_to_convert, allmet_converted, "metabolites",
                                             models_same_db)
    allreact_selected = selection.runSelection(models_to_convert, allreact_converted, "reactions",
                                               models_same_db)
    structural_r_all, structural_r_one = structural.runAdditionalConversion(models_to_convert, allmet_selected,
                                                                           allreact_selected, curated_models,
                                                                           bigg_db_network, "reactions")
    '''
    # check structural reaction consistency
    tmp_r_one_consist, struct_r_one_not_consist = selection.checkDBConsistency(models_same_db, structural_r_one,
                                                                                "reactions", write_files=False,
                                                                               do_stat=False)
    struct_r_one_consist = deepcopy(structural_r_one)
    for same_models in models_same_db.values():
        for model in same_models:
            struct_r_one_consist[model] = tmp_r_one_consist[model]
    # check structural reactions many_to_one
    for typ in models_to_convert:
        struct_r_one_consist[typ] = {k: [v[0], v[1][0]] for k, v in struct_r_one_consist[typ].items()}
    structural_r_one_one, structural_r_many_one = selection.checkFromOneFromMany(models_to_convert, struct_r_one_consist)



    # make suggestions based on one_t_one consistent reactions for one_to_many metabolites
    suggestions_one_many_m, suggestions_one_many_m_sel = selection.getSuggestionForMetabolites(models_to_convert,
                                                                                               structural_r_one_one,
                                                                                               structural_r_all)
    # unite suggestions coming for models from the same database
    met_same_db = set()
    for same_models in models_same_db.values():
        for model in same_models:
            met_same_db = met_same_db | set(suggestions_one_many_m_sel.get(model).keys())
    unite_suggestions_one_many_m_sel = {}
    for typ in models_to_convert:
        unite_suggestions_one_many_m_sel[typ] = {
            key: [allmet_selected.get("intermediate_data").get("highest").get(typ).get(key)[0], val] for key, val in
            suggestions_one_many_m_sel.get(typ).items()}
        if typ in allmet_selected.get("intermediate_data").get("consistent").keys():
            for met_to_add in met_same_db:
                if (met_to_add not in list(suggestions_one_many_m_sel.get(typ).keys())) & (
                        met_to_add in list(allmet_selected.get("intermediate_data").get("consistent").get(typ).keys())):
                    unite_suggestions_one_many_m_sel.get(typ).update({met_to_add: allmet_selected.get(
                        "intermediate_data").get("consistent").get(typ).get(met_to_add)})
    # check metabolite suggestions for consistency
    tmp_m_sug_consist, m_sug_not_consist = selection.checkDBConsistency(models_same_db,
                                                                        unite_suggestions_one_many_m_sel,
                                                                         "metabolites", write_files=False,
                                                                        do_stat=False)
    # check metabolites suggestions for many_to_one
    met_suggestions = deepcopy(unite_suggestions_one_many_m_sel)
    for same_models in models_same_db.values():
        for model in same_models:
            met_suggestions[model] = tmp_m_sug_consist[model]
    tmp_first_structural_met = {}
    for typ in models_to_convert:
        met_suggestions[typ] = {k: [v[0], v[1][0], {"selection_type": "structural_suggestions_for_one_many"}] for k, v
                                in met_suggestions[typ].items()}
        tmp_first_structural_met.update({typ: met_suggestions[typ]})
        tmp_first_structural_met.get(typ).update(
            {key: [val[0], val[1], {"selection_type": "selection_one_one"}] for key, val in
             allmet_selected.get("one_to_one").get(typ).items()})
    first_structural_met, first_structural_met_many_one = selection.checkFromOneFromMany(models_to_convert,
                                                                                         tmp_first_structural_met)



    # make suggestions for many_to_one metabolite
    many_to_one_suggestions = {}
    for typ in models_to_convert:
        many_to_one_suggestions.update({typ:
                                            selection.ManyOneMetFromStructuralMet(
                                                allmet_selected.get("many_to_one").get(typ), curated_models.get(typ),
                                                first_structural_met.get(typ), bigg_db_network.get("reactions"))})
   # unite suggestions coming for models from the same database
    met_same_db_mo = set()
    for same_models in models_same_db.values():
        for model in same_models:
            met_same_db_mo = met_same_db_mo | set(many_to_one_suggestions.get(model).keys())
    unite_suggestions_many_one_m_sel = {}
    for typ in models_to_convert:
        unite_suggestions_many_one_m_sel[typ] = {k: [v[0], [v[1]]] for k, v in many_to_one_suggestions.get(typ).items()}
        if typ in allmet_selected.get("intermediate_data").get("consistent").keys():
            for met_to_add in met_same_db_mo:
                if (met_to_add not in list(many_to_one_suggestions.get(typ).keys())) & (
                        met_to_add in list(allmet_selected.get("intermediate_data").get("consistent").get(typ).keys())):
                    unite_suggestions_many_one_m_sel.get(typ).update({met_to_add: allmet_selected.get(
                        "intermediate_data").get("consistent").get(typ).get(met_to_add)})
    # check metabolite suggestions for consistency
    tmp_mo_m_sug_consist, mo_m_sug_not_consist = selection.checkDBConsistency(models_same_db,
                                                                              unite_suggestions_many_one_m_sel,
                                                                         "metabolites", write_files=False,
                                                                              do_stat=False)
    # check metabolites suggestions for many_to_one
    mo_met_suggestions = deepcopy(unite_suggestions_many_one_m_sel)
    for same_models in models_same_db.values():
        for model in same_models:
            mo_met_suggestions[model] = tmp_mo_m_sug_consist[model]
    tmp_final_met = {}
    for typ in models_to_convert:
        mo_met_suggestions[typ] = {k: [v[0], v[1][0], {"selection_type": "structural_suggestions_for_one_many"}] for k, v
                                in mo_met_suggestions[typ].items()}
        tmp_final_met.update({typ: mo_met_suggestions[typ]})
        tmp_final_met.get(typ).update(first_structural_met.get(typ))
    final_metabolites, final_metabolites_many_one = selection.checkFromOneFromMany(models_to_convert,
                                                                                   tmp_final_met)
    test_r, test_r_one = selection.test(models_to_convert,final_metabolites, allmet_selected, allreact_selected, curated_models, bigg_db_network, "reactions")
    # check consistency for rtest
    tmp_r_consist, struct_r_not_consist = selection.checkDBConsistency(models_same_db, test_r_one,
                                                                                "reactions", write_files=False,
                                                                       do_stat=False)
    struct_r_consist = deepcopy(test_r_one)
    for same_models in models_same_db.values():
        for model in same_models:
            struct_r_consist[model] = tmp_r_consist[model]
    # check structural reactions many_to_one
    for typ in models_to_convert:
        struct_r_consist[typ] = {k: [v[0], v[1][0]] for k, v in struct_r_consist[typ].items()}
    test_structural_r_one_one_one, test_structural_r_many_one = selection.checkFromOneFromMany(models_to_convert, struct_r_consist)
    '''
