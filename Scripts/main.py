import copy
import os
from os.path import join, exists
import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model
import conversion
import selection
import curation
import BiGGnetwork
import dill

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
    allmet_converted = conversion.runConversionForALLmodels(models_to_convert, curated_models, CompartmentsStrategies,
                                                            ConversionStrategies, "metabolites")
    allreact_converted = conversion.runConversionForALLmodels(models_to_convert, curated_models, CompartmentsStrategies,
                                                              ConversionStrategies, "reactions")
    allmet_selected = selection.selectBasedOnConversionQuality(models_to_convert, allmet_converted, "metabolites",
                                                                models_same_db)
    allreact_selected = selection.selectBasedOnConversionQuality(models_to_convert, allreact_converted, "reactions",
                                                                  models_same_db)
    structural_r = selection.runAdditionalConversion(models_to_convert, allmet_selected, allreact_selected, curated_models,
                                                     bigg_db_network, "reactions")
