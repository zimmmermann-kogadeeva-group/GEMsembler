import os
from os.path import join
import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model
import conversion
import selection

if __name__ == '__main__':
    # region Open conversion tables
    fileDir = os.path.dirname(os.path.realpath('__file__'))  # getting directory of the script for pathes to files
    seed_orig = pd.read_csv(join(fileDir, "../../../Databases/metabolic_ids_convertion/seed_to_bigg.tsv"),
                            sep="\t")
    seed_orig_m = seed_orig[seed_orig["type"] == "m"][["seed_ids", "bigg_ids"]]
    seed_orig_m.columns = ["old", "new"]
    seed_orig_r = seed_orig[seed_orig["type"] == "r"][["seed_ids", "bigg_ids"]]
    seed_orig_r.columns = ["old", "new"]
    seed_addit = pd.read_csv(join(fileDir, "../../../Databases/metabolic_ids_convertion/seed_to_bigg_metanetx.tsv"),
                            sep="\t")
    seed_addit_m = seed_addit[seed_addit["type"] == "m"][["seed_ids", "bigg_ids"]]
    seed_addit_m.columns = ["old", "new"]
    seed_addit_r = seed_addit[seed_addit["type"] == "r"][["seed_ids", "bigg_ids"]]
    seed_addit_r.columns = ["old", "new"]
    kegg_bigg = pd.read_csv(join(fileDir, "../../../Databases/metabolic_ids_convertion/kegg_to_bigg_metanetx.tsv"),
                            sep="\t")
    kegg_bigg_m = kegg_bigg[kegg_bigg["type"] == "m"][["kegg_ids", "bigg_ids"]]
    kegg_bigg_m.columns = ["old", "new"]
    kegg_bigg_r = kegg_bigg[kegg_bigg["type"] == "r"][["kegg_ids", "bigg_ids"]]
    kegg_bigg_r.columns = ["old", "new"]
    old_new_bigg = pd.read_csv(join(fileDir, "../../../Databases/metabolic_ids_convertion/old_to_new_bigg.tsv"),
                               sep="\t")
    old_new_bigg_m = old_new_bigg[old_new_bigg["type"] == "m"][["old_bigg_ids", "bigg_ids"]]
    old_new_bigg_m.columns = ["old", "new"]
    old_new_bigg_m["new"] = old_new_bigg_m["new"].str[:-2]
    old_new_bigg_r = old_new_bigg[old_new_bigg["type"] == "r"][["old_bigg_ids", "bigg_ids"]]
    old_new_bigg_r.columns = ["old", "new"]
    old_new_bigg_r["new"] = old_new_bigg_r[
        "new"]  # .str.replace("_[cep]", "", regex = True) - not sure, whether it  is reasnable
    bigg_all_m = pd.read_csv(join(fileDir, "../../../Databases/BiGG/bigg_models_metabolites.txt"), sep="\t")
    bigg_m = list(set(bigg_all_m["universal_bigg_id"]))
    bigg_all_r = pd.read_csv(join(fileDir, "../../../Databases/BiGG/bigg_models_reactions.txt"), sep="\t")
    bigg_all_r["universal_bigg_id"] = bigg_all_r["bigg_id"].str.replace("_[cep]", "", regex=True)
    bigg_r = list(set(bigg_all_r["universal_bigg_id"]))

    # endregion

    # region Open models
    model_type_list = ["carveme", "gapseq", "modelseed", "agora"]
    models_same_db = {"modelseed": ["gapseq", "modelseed"]}
    fileDir = os.path.dirname(os.path.realpath('__file__'))  # getting directory of the script for paths to files
    name_carv_hom = join(fileDir, "../Data/BU_carveme_hom.xml")
    name_gapseq = join(fileDir, "../Data/BU_gapseq.xml")
    name_modelseed = join(fileDir, "../Data/BU_modelSEED.sbml")
    name_agora = join(fileDir, "../Data/BU_agora.xml")
    names = {"carveme": name_carv_hom, "gapseq": name_gapseq, "modelseed": name_modelseed, "agora": name_agora}
    all_models = {}
    for typ in model_type_list:
        nam = names.get(typ)
        model = read_sbml_model(nam)
        all_models.update({typ: model})

    # endregion

    # region Perform conversion
    GapseqConv = conversion.ConversionForGapseq(seed_orig_m, seed_orig_r, seed_addit_m, seed_addit_r, bigg_m, bigg_r)
    ModelseedConv = conversion.ConversionForModelseed(seed_orig_m, seed_orig_r, seed_addit_m, seed_addit_r, bigg_m, bigg_r)
    AgoraConv = conversion.ConversionForAgora(old_new_bigg_m, old_new_bigg_r, kegg_bigg_m, kegg_bigg_r, bigg_m, bigg_r)
    ConversionStrategies = {"gapseq": GapseqConv, "modelseed": ModelseedConv, "agora": AgoraConv}
    CarvemeComp = conversion.CompartmentsForCarveme()
    GapseqComp = conversion.CompartmentsForGapseq()
    ModelseedComp = conversion.CompartmentsForModelseed()
    AgoraComp = conversion.CompartmentsForAgora()
    CompartmentsStrategies = {"carveme": CarvemeComp, "gapseq": GapseqComp, "modelseed": ModelseedComp,
                              "agora": AgoraComp}

    models_to_convert = model_type_list[1:]
    allmet_converted = conversion.runConversionForALLmodels(models_to_convert, all_models, CompartmentsStrategies,
                                                              ConversionStrategies, "metabolites")
    allreact_converted = conversion.runConversionForALLmodels(models_to_convert, all_models, CompartmentsStrategies,
                                                                ConversionStrategies, "reactions")
    highest_m = selection.getHighestConversion(models_to_convert, allmet_converted)
    highest_r = selection.getHighestConversion(models_to_convert, allreact_converted)
    consist_m = selection.checkSameConversion(models_same_db, highest_m, "metabolites")
    consist_r = selection.checkSameConversion(models_same_db, highest_r, "reactions")