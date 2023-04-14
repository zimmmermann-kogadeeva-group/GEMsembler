import cobra
from cobra.io import read_sbml_model
from os.path import join, dirname, realpath
import pandas as pd


def checkDuplicatedReactions(model_types: [str], all_models: dict,):
    """ Getting tables with reactions with the same equation/same equation + same GPR for each model. """
    duplicated = {}
    for typ in model_types:
        data = pd.DataFrame(columns=["ID", "Reactants", "Products", "GPR"])
        for r in all_models.get(typ).reactions:
            reactants = sorted([react.id for react in r.reactants])
            products = sorted([pro.id for pro in r.products])
            if (r.lower_bound == 0) and (r.upper_bound > 0):
                data = pd.concat(
                    [data, pd.DataFrame([[r.id, " ".join(reactants), " ".join(products), r.gene_reaction_rule]],
                                        columns=data.columns)],
                    ignore_index=True)
            if (r.lower_bound < 0) and (r.upper_bound == 0):
                data = pd.concat(
                    [data, pd.DataFrame([[r.id, " ".join(products), " ".join(reactants), r.gene_reaction_rule]],
                                        columns=data.columns)],
                    ignore_index=True)
            if (r.lower_bound < 0) and (r.upper_bound > 0):
                data = pd.concat(
                    [data, pd.DataFrame([[r.id, " ".join(reactants), " ".join(products), r.gene_reaction_rule]],
                                        columns=data.columns)],
                    ignore_index=True)
                data = pd.concat(
                    [data, pd.DataFrame([[r.id, " ".join(products), " ".join(reactants), r.gene_reaction_rule]],
                                        columns=data.columns)],
                    ignore_index=True)
        struct_duplicated = data[data.duplicated(subset=['Reactants', 'Products'], keep=False)]
        struct_duplicated = struct_duplicated.sort_values(["Reactants", "Products"])
        gpr_duplicated = data[data.duplicated(subset=['Reactants', 'Products', "GPR"], keep=False)]
        gpr_duplicated = gpr_duplicated.sort_values(["Reactants", "Products", "GPR"])
        duplicated.update({typ: [struct_duplicated, gpr_duplicated]})
    return duplicated


def removeBtypeExchange(curated_model: cobra.core.model.Model):
    """ Removing boundary metabolites / exchange reactions with _b at the end if another (_e) compartment exist.
     Replacing _b with normal compartment if id with normal compartment doesn't exist separately. """
    all_met_ids = [m.id for m in curated_model.metabolites]
    b_met_ids = [m.id for m in curated_model.metabolites if m.id.endswith("_b")]
    all_r_ids = [r.id for r in curated_model.reactions]
    b_r_ids = [r.id for r in curated_model.reactions if r.id.endswith("_b")]
    met_to_remove = []
    r_to_remove = []
    for b_met in b_met_ids:
        b_met_e = b_met.removesuffix("_b") + "_" + curated_model.metabolites.get_by_id(b_met).compartment
        if b_met_e in all_met_ids:
            met_to_remove.append(curated_model.metabolites.get_by_id(b_met))
        else:
            curated_model.metabolites.get_by_id(b_met).id = b_met_e
            curated_model.repair()
    for b_r in b_r_ids:
        b_r_e = b_r.removesuffix("_b") + "_" + list(curated_model.reactions.get_by_id(b_r).compartments)[0]
        # TODO: make an approach in case compartments in specific filed is not the same as compartments at the end of id
        if b_r_e in all_r_ids:
            r_to_remove.append(curated_model.reactions.get_by_id(b_r))
        else:
            curated_model.reactions.get_by_id(b_r).id = b_r_e
            curated_model.repair()
    curated_model.remove_metabolites(met_to_remove)
    curated_model.remove_reactions(r_to_remove)
    return curated_model


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
    duplicated = checkDuplicatedReactions(model_type_list, all_models)
    for typ in model_type_list:
        print(duplicated.get(typ)[0])
        print(duplicated.get(typ)[1])
