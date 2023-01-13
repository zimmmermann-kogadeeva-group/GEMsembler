import cobra
from cobra.io import read_sbml_model
from os.path import join, dirname, realpath
import pandas as pd
import bisect


def checkDuplicatedReactions(model: cobra.core.model.Model):
    data = pd.DataFrame(columns=["ID", "Reactants", "Products", "GPR"])
    for r in model.reactions:
        reactants = []
        products = []
        for react in r.reactants:
            bisect.insort(reactants, react.id)
        for pro in r.products:
            bisect.insort(products, pro.id)
        if (r.lower_bound == 0) and (r.upper_bound > 0):
            data = pd.concat([data, pd.DataFrame([[r.id, str(reactants), str(products), r.gene_reaction_rule]],
                                                 columns=data.columns)],
                             ignore_index=True)
        if (r.lower_bound < 0) and (r.upper_bound == 0):
            data = pd.concat([data, pd.DataFrame([[r.id, str(products), str(reactants), r.gene_reaction_rule]],
                                                 columns=data.columns)],
                             ignore_index=True)
        if (r.lower_bound < 0) and (r.upper_bound > 0):
            data = pd.concat([data, pd.DataFrame([[r.id, str(reactants), str(products), r.gene_reaction_rule]],
                                                 columns=data.columns)],
                             ignore_index=True)
            data = pd.concat([data, pd.DataFrame([[r.id, str(products), str(reactants), r.gene_reaction_rule]],
                                                 columns=data.columns)],
                             ignore_index=True)
    struct_duplicated = data[data.duplicated(subset=['Reactants', 'Products'], keep=False)]
    struct_duplicated = struct_duplicated.sort_values(["Reactants", "Products"])
    gpr_duplicated = data[data.duplicated(subset=['Reactants', 'Products', "GPR"], keep=False)]
    gpr_duplicated = gpr_duplicated.sort_values(["Reactants", "Products", "GPR"])
    return struct_duplicated, gpr_duplicated


if __name__ == '__main__':
    # region Open models
    model_type_list = ["carveme", "gapseq", "modelseed", "agora"]
    fileDir = dirname(realpath('__file__'))  # getting directory of the script for paths to files
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
    for typ in model_type_list:
        struct_duplicated, gpr_duplicated = checkDuplicatedReactions(all_models.get(typ))
        print(struct_duplicated)
        print(gpr_duplicated)
