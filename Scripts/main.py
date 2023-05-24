import copy
import operator
import os
from os.path import join, exists
import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model

import anticreation
import curation
import BiGGnetwork
import conversion
import drawing
import general
import selection
import structural
import creation
import comparison
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

    # Conversion strategies
    GapseqConv = conversion.ConversionForGapseq(seed_orig_m, seed_orig_r, seed_addit_m, seed_addit_r, bigg_m, bigg_r)
    ModelseedConv = conversion.ConversionForModelseed(seed_orig_m, seed_orig_r, seed_addit_m, seed_addit_r, bigg_m,
                                                      bigg_r)
    AgoraConv = conversion.ConversionForAgora(old_new_bigg_m, old_new_bigg_r, kegg_bigg_m, kegg_bigg_r, bigg_m, bigg_r)
    CarvemeConv = conversion.ConversionForCarveMe(old_new_bigg_m, old_new_bigg_r, bigg_m, bigg_r)
    ConversionStrategies = {"gapseq": GapseqConv, "modelseed": ModelseedConv, "agora": AgoraConv,
                            "carveme": CarvemeConv}
    allmet_converted_dill_file = join(fileDir, "../Scripts/allmet_converted.pkl")
    allreact_converted_dill_file = join(fileDir, "../Scripts/allreact_converted.pkl")
    if exists(allmet_converted_dill_file) & exists(allreact_converted_dill_file):
        allmet_converted = dill.load(open(allmet_converted_dill_file, "rb"))
        allreact_converted = dill.load(open(allreact_converted_dill_file, "rb"))
    else:
        # converting ids for models not with BiGG ids
        allmet_converted = conversion.runConversionForALLmodels(models_to_convert, curated_models,
                                                                ConversionStrategies, "metabolites")
        allreact_converted = conversion.runConversionForALLmodels(models_to_convert, curated_models,
                                                                  ConversionStrategies, "reactions")
        dill.dump(allmet_converted, open(allmet_converted_dill_file, "wb"))
        dill.dump(allreact_converted, open(allreact_converted_dill_file, "wb"))
    # checking ids from models that supposed to be with BiGG but still may be old (or potentialy wrond but that i don't remember)
    allmet_checked, allmet_not_pass = conversion.runNoneConversionChecking(models_NOTto_convert, curated_models,
                                                                           ConversionStrategies, "metabolites")
    allreact_checked, allreact_not_pass = conversion.runNoneConversionChecking(models_NOTto_convert, curated_models,
                                                                               ConversionStrategies, "reactions")
    # checking reaction equations for models with BiGG
    allreact_checked_struct, allreact_not_pass_struct = structural.runStructuralCheck(models_NOTto_convert,
                                                                                      allreact_checked,
                                                                                      allreact_not_pass, curated_models,
                                                                                      bigg_db_network)
    # geting select ids that are converted 1-1, 1-n, n-1, n-n from converted ids (and check for consistency in models that are deferent but use ids from the same databes (modelseed in this case))
    allmet_selected = selection.runSelection(models_to_convert, allmet_converted, "metabolites",
                                             models_same_db)
    allreact_selected = selection.runSelection(models_to_convert, allreact_converted, "reactions",
                                               models_same_db)
    # converting reactions via reactions equations from BiGG database via metabolites that were converted 1-1 or 1-n and selecting ones that were converted with equations and were converted uniquely
    structural_r_info, structural_r_sel = structural.runStructuralConversion(models_to_convert,
                                                                             allmet_selected.get("one_to_one"),
                                                                             allmet_selected, allreact_selected,
                                                                             curated_models,
                                                                             bigg_db_network, models_wo_periplasmic,
                                                                             allmet_selected.get("one_to_many"))
    # checking consistency in reactions uniquely converted via equation for different models with same database
    struct_r_consistent, struct_r_consist, struct_r_not_consist = selection.checkDBConsistency(models_same_db,
                                                                                               structural_r_sel,
                                                                                               "reactions",
                                                                                               write_files=False,
                                                                                               do_stat=False)
    # selecting which consistent structural reactions are converted 1-1 and n-1
    struct_r_uniq, struct_r_not_uniq = selection.checkFromOneFromMany(models_to_convert, struct_r_consistent)
    # getting suggestions for 1-n metabolites from first run of structural conversion and tring to structural with n-1 metabolites. Get 1-1 metabolites plus suggestions
    met_struct = structural.runSuggestionsMet(models_to_convert, structural_r_info, struct_r_uniq, allmet_selected,
                                              models_same_db, curated_models, bigg_db_network)
    # running structural reaction conversion second time with update metabolites (1-1 plus suggestions)
    struct_final_r_info, struct_final_r_sel = structural.runStructuralConversion(models_to_convert,
                                                                                 met_struct.get("one_one_sugg_met"),
                                                                                 allmet_selected, allreact_selected,
                                                                                 curated_models, bigg_db_network,
                                                                                 models_wo_periplasmic)
    # checking consistency in 2d run of reactions uniquely converted via equation for different models with same database
    struct_final_r_consistent, struct_final_r_consist, struct_final_r_not_consist = selection.checkDBConsistency(
        models_same_db, struct_final_r_sel,
        "reactions", write_files=False, do_stat=False)
    # selecting which consistent structural reactions from run 2 are converted 1-1 and n-1
    struct_final_r_uniq, struct_final_r_not_uniq = selection.checkFromOneFromMany(models_to_convert,
                                                                                  struct_final_r_consistent)
    # getting metabolites and reactions that became preiplasmic for models without periplasmic compartment originally
    periplasmic_m, periplasmic_r = structural.getSuggestionPeriplasmic(models_wo_periplasmic, struct_final_r_uniq,
                                                                       struct_final_r_info, bigg_db_network,
                                                                       curated_models)
    # getting in standart format final metabolites and reactions for supermodel creation (converted)
    final_r = deepcopy(struct_final_r_uniq)
    final_r_not_uniq = {}
    for typ in models_to_convert:  # adding n-1 (duplicated) reactions structural reactions if they are indeed originaly duplicated in models
        true_dupl = list(
            set(struct_final_r_not_uniq.get(typ).keys()) & set(duplicated_reactions.get(typ)[0]["ID"].tolist()))
        if true_dupl:
            final_r.get(typ).update({td: struct_final_r_not_uniq.get(typ).get(td) for td in true_dupl})
        false_dupl = list(
            set(struct_final_r_not_uniq.get(typ).keys()) - set(duplicated_reactions.get(typ)[0]["ID"].tolist()))
        if false_dupl:
            final_r_not_uniq.update({typ: {fd: struct_final_r_not_uniq.get(typ).get(fd) for fd in false_dupl}})
    # getting in standard format not selected for not converted in supermodel
    final_m_not_sel = selection.runNotSelectedMet(models_to_convert, met_struct, allmet_selected)
    final_r_not_sel = selection.runNotSelectedR(models_to_convert, final_r, struct_final_r_not_consist,
                                                final_r_not_uniq, struct_final_r_info, curated_models)
    final_m = deepcopy(met_struct.get("one_one_sugg_met"))
    additional_p_m = {}
    # dealing with periplasmic metabolites: replacing if final and creating additional periplasmic metabolites if original metabolite works with boths compartments
    for typ in models_wo_periplasmic:
        additional_p_m.update({typ: {}})
        for orig_id in periplasmic_m.get(typ).keys():
            if periplasmic_m.get(typ).get(orig_id)[4] == "replace":
                final_m.get(typ)[orig_id] = [final_m.get(typ).get(orig_id)[0], [periplasmic_m.get(typ).get(orig_id)[1]]]
            else:
                additional_p_m.get(typ).update(
                    {orig_id: [final_m.get(typ).get(orig_id)[0], [periplasmic_m.get(typ).get(orig_id)[1]]]})
    # combining final metabolites and reactions with ones coming from BiGG models (wo conversion)
    for typ in models_NOTto_convert:
        final_m.update({typ: allmet_checked.get(typ)})
        final_r.update({typ: allreact_checked_struct.get(typ)})
        m_notsel = general.findKeysByValue(allmet_not_pass.get(typ), "not_found_in_new_and_old_bigg", operator.eq)
        if m_notsel:
            final_m_not_sel.update({typ: {m: [allmet_not_pass.get(typ).get(m)[0], m] for m in m_notsel}})
        else:
            final_m_not_sel.update({typ: {}})
        r_notsel = [k for k, v in allreact_not_pass_struct.get(typ).items() if v[1] == "not_found_in_new_and_old_bigg"]
        if r_notsel:
            final_r_not_sel.update({typ: {r: [allreact_not_pass_struct.get(typ).get(r)[0], [r]] for r in r_notsel}})
        else:
            final_r_not_sel.update({typ: {}})
        for kg, vg in allreact_not_pass_struct.get(typ).items():
            if vg[1] == "Growth_reaction":
                final_r_not_sel.get(typ).update({kg: [allreact_not_pass_struct.get(typ).get(kg)[0], ["Biomass"]]})
    # creating supermodel
    supermodel = creation.runSupermodelCreation(model_type_list, final_m, final_m_not_sel, final_r, final_r_not_sel,
                                                curated_models, bigg_all_m, bigg_all_r, additional_p_m, periplasmic_r)
    # getting core and different types of intersections in supermodel
    comparison.runComparioson(supermodel)
    core_model = anticreation.getModelOfInterest(supermodel, "core4", name="BU_core_model.xml")
    union_model = anticreation.getModelOfInterest(supermodel, "union1", name="BU_union_model.xml")
    Yes_a_No_cmg_model = anticreation.getModelOfInterest(supermodel, "Yes_a_No_cgm", name="BU_Yes_a_No_cmg_model.xml")
    carveme_model = anticreation.getModelOfInterest(supermodel, "carveme", name="BU_carveme_out_model.xml")
    # plotting some test pathways and core (intersection)
    # colorBrewer = {"reds": ["#feedde", "#fdbe85", "#fd8d3c", "#d94701"],
    #                "blues": ["#eff3ff", "#bdd7e7", "#6baed6", "#2171b5"],
    #                "purples": ["#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3"],
    #                "greens": ["#edf8e9", "#bae4b3", "#74c476", "#238b45"]}
    # aminoacids = ["ala__L", "arg__L", "asn__L", "asp__L", "cys__L", "gln__L", "glu__L", "gly", "his__L", "ile__L",
    #               "leu__L",
    #               "lys__L", "met__L", "phe__L", "pro__L", "ser__L", "thr__L", "trp__L", "tyr__L", "val__L"]
    # met_not_int = ["h", "h2o", "h2", "oh1", "o2", "co2", "coa", "ppi", "pi", "amp", "adp", "atp", "cmp", "cdp", "ctp",
    #                "gmp", "gdp", "gtp", "ump", "udp", "utp",
    #                "nad", "nadh", "nadp", "nadph", "dadp", "damp", "nh3", "nh4", "fadh2", "fad",
    #                "ac", "accoa", "h2s", "HC00250"]
    #
    # glycolisys = {
    #     "metabolites": ["glc__D_c", "g6p_c", "f6p_c", "fdp_c", "g3p_c", "13dpg_c", "3pg_c", "2pg_c", "pep_c", "pyr_c"],
    #     "reactions": ["HEX1", "PGI", "PFK", "FBA", "GAPD", "PGK", "PGM", "ENO", "PYK"]}
    # glycol = drawing.drawOnePathway(supermodel, glycolisys, met_not_int, colorBrewer, "glycolysis", aminoacids,
    #                                     directed=False, surrounding=True)
    # cys_syn1 = {
    #     "metabolites": ["glc__D_c", "g6p_c", "f6p_c", "fdp_c", "g3p_c", "13dpg_c", "3pg_c", "2pg_c", "pep_c", "pyr_c",
    #                     "oaa_c", "asp__L_c",
    #                     "4pasp_c", "aspsa_c",
    #                     "hom__L_c", "achms_c", "hcys__L_c",
    #                     "cyst__L_c", "cys__L_c"
    #                     ],
    #     "reactions": ["HEX1", "PGI", "PFK", "FBA", "GAPD", "PGK", "PGM", "ENO", "PYK",
    #                   "PC", "ASPTA",
    #                   "ASPK", "ASAD",
    #                   "HSDxi", "HSDy", "HSERTA", "AHSERL2",  # 2 first R - alternatives in NAD / NADP
    #                   "CYSTS", "CYSTGL",
    #                   ]}
    # cys_syn = {"metabolites": ["glc__D_c", "g6p_c", "f6p_c", "fdp_c", "g3p_c", "13dpg_c", "3pg_c",
    #                            "3php_c", "pser__L_c", "ser__L_c",
    #                            "acser_c", "cys__L_c"],
    #            "reactions": ["HEX1", "PGI", "PFK", "FBA", "GAPD", "PGK",
    #                          "PGCD", "PSERT", "PSP_L",
    #                          "SERAT", "CYSS"]
    #            }
    # cys0 = drawing.drawOnePathway(supermodel, cys_syn, met_not_int, colorBrewer, "cys_main")
    # cys_alt = drawing.drawTwoPathways(supermodel, cys_syn, cys_syn1, met_not_int, colorBrewer, "cys_altern_paths")
    # TCA = {"PYK": [("pep_c", "pyr_c")],
    #        "PPC": [("pep_c", "oaa_c")], "PPCK": [("pep_c", "oaa_c")], "PEPCK_re": [("pep_c", "oaa_c")],
    #        "PC": [("pyr_c", "oaa_c")],
    #        "PDH": [("pyr_c", "accoa_c")], "PFL": [("pyr_c", "accoa_c")],
    #        "CS": [("accoa_c", "cit_c"), ("oaa_c", "cit_c")], "ACONT": [("cit_c", "icit_c")],
    #        "ACONTa": [("cit_c", "acon_C_c")], "ACONTb": [("acon_C_c", "icit_c")],
    #        "ICL": [("icit_c", "succ_c"), ("icit_c", "glx_c")],
    #        "ICDHyr": [("icit_c", "akg_c")], "ICDHx": [("icit_c", "akg_c")], "ICITRED": [("icit_c", "osuc_c")],
    #        "OSUCCL": [("osuc_c", "akg_c")],
    #        "AKGDH": [("akg_c", "succoa_c")], "OOR2r": [("akg_c", "succoa_c"), ("fdxo_42_c", "fdxr_42_c")],
    #        "AKGDa": [("akg_c", "sdhlam_c"), ("lpam_c", "sdhlam_c")],
    #        "AKGDb": [("sdhlam_c", "succoa_c"), ("sdhlam_c", "dhlam_c")], "PDHcr": [("dhlam_c", "lpam_c")],
    #        "SUCOAS": [("succoa_c", "succ_c")],
    #        "SUCDi": [("succ_c", "fum_c")], "FRD7": [("fum_c", "succ_c")],
    #        "FUM": [("fum_c", "mal__L_c")], "MALS": [("glx_c", "mal__L_c")],
    #        "MDH": [("mal__L_c", "oaa_c")], "MDH2": [("mal__L_c", "oaa_c")], "MDH3": [("mal__L_c", "oaa_c")]}
    # tca_g = drawing.drawTCA(supermodel, TCA, met_not_int, colorBrewer, "TCA")
    # core = drawing.drawCore(supermodel, met_not_int, colorBrewer, "new_core")
    # print(
    #     f"BU core consist of {len(supermodel.reactions.core4)} reactions and {len(supermodel.metabolites.core4)} metabolites")
    # union = drawing.drawCore(supermodel, met_not_int, colorBrewer, "new_union", union=True)
    # print(
    #     f"BU supermodel consist of {len(supermodel.reactions.converted)} reactions and {len(supermodel.metabolites.converted)} metabolites")
    # biomass = drawing.drawBiomass(supermodel, "biomass", colorBrewer=colorBrewer)
    # biomass_diff = drawing.drawBiomass(supermodel, "biomass_diff", only_difference=True, colorBrewer=colorBrewer)
    # biomass_notconv = drawing.drawBiomass(supermodel, "biomass_notconverted", not_converted=True)
