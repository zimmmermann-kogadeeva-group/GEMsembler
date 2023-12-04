import warnings
from collections import defaultdict
from pathlib import PosixPath
from zipfile import Path
from .genes import makeNewGPR, uniteGPR
import pandas as pd


class NewObject:
    """ New object class - one metabolite or reaction for supermodel. """

    def __init__(
        self,
        new_id: str,
        old_id: str,
        compartments: [str],
        source: str,
        possible_sources: [str],
    ):
        self.id = new_id
        self.compartments = {"assembly": compartments}
        self.sources = {}
        self.in_models = {"models_amount": 1, "models_list": [source]}
        self.annotation = {}
        for ps in possible_sources:
            if ps == source:
                self.compartments.update({ps: compartments})
                self.sources.update({ps: 1})
                self.annotation.update({ps: [old_id]})
            else:
                self.compartments.update({ps: []})
                self.sources.update({ps: 0})
                self.annotation.update({ps: []})

    def updateNewObject(
        self, id_to_update: str, compart_to_update: [str], source: str,
    ):
        self.sources.update({source: self.sources.get(source) + 1})
        if source not in self.in_models["models_list"]:
            self.in_models["models_amount"] = self.in_models["models_amount"] + 1
            self.in_models["models_list"].append(source)
        self.annotation.get(source).append(id_to_update)
        self.compartments.update(
            {source: self.compartments.get(source) + compart_to_update}
        )
        for c in compart_to_update:
            if c not in self.compartments["assembly"]:
                self.compartments["assembly"].append(c)


class SetofNewObjects:
    """ Setting dictionaries for all metabolites or reactions:
    selected for supermodel - self.converted and not selected - self.notconverted. """

    def addNewObjs(self, selected: dict, where_to_add: str, model_ids: list):
        for mod_id in model_ids:
            if mod_id in list(selected.keys()):
                objects = selected.get(mod_id)
                for key in objects.keys():
                    for new_id in objects[key][1]:
                        comp = objects[key][0]
                        if new_id in getattr(self, where_to_add).keys():
                            getattr(self, where_to_add).get(new_id).updateNewObject(
                                key, comp, mod_id
                            )
                        else:
                            new = NewObject(new_id, key, comp, mod_id, model_ids)
                            getattr(self, where_to_add).update({new_id: new})

    def makeSetofNew(
        self, selected: dict, not_selected: dict, model_ids, additional,
    ):
        self.addNewObjs(selected, "assembly_conv", model_ids)
        if additional:
            self.addNewObjs(additional, "assembly_conv", model_ids)
        for new_id, new_obj in self.assembly_conv.items():
            for model_id in new_obj.in_models["models_list"]:
                self.comparison[model_id].update({new_id: new_obj})
        self.addNewObjs(
            not_selected, "notconverted", model_ids
        )  # TODO connect not_converted for really not converted only with old id

    def __init__(
        self, selected: dict, not_selected: dict, model_ids: [str], additional=None,
    ):
        self.assembly_conv = {}
        self.assembly_mix = {}
        self.comparison = defaultdict(dict)
        self.notconverted = {}
        self.makeSetofNew(selected, not_selected, model_ids, additional)

    def makeForwardBackward(
        self,
        all_models: dict,
        selected: dict,
        obj_type: "metabolites" or "reactions",
        additional=None,
    ):
        """ Creating dictionaries linking metabolites/reactions:
        NewObject in supermodel with old original ID and OldObject in original models with new ID in supermodel """
        go_old_new = defaultdict(dict)
        go_new_old = defaultdict(dict)
        model_ids = list(selected.keys())
        for model_id in model_ids:
            for key, value in selected.get(model_id).items():
                new_obj = [self.assembly_conv.get(value[1][0])]
                if additional:
                    if model_id in list(additional.keys()):
                        if key in list(additional.get(model_id).keys()):
                            new_obj.append(
                                self.assembly_conv.get(
                                    additional.get(model_id).get(key)[1][0]
                                )
                            )
                go_old_new[model_id].update({key: new_obj})
        for k, v in self.assembly_conv.items():
            for mod_id in v.in_models["models_list"]:
                old_ids = v.annotation[mod_id]
                old_obj = [
                    getattr(all_models[mod_id]["preprocess_model"], obj_type).get_by_id(
                        i
                    )
                    for i in old_ids
                ]
                go_new_old[k].update({mod_id: old_obj})
        return go_old_new, go_new_old


class SetofNewMetabolites(SetofNewObjects):
    """ Metabolites class that add name and blank reaction attribute to metabolite """

    def setMetaboliteAttributes(self, database_info: pd.core.frame.DataFrame):
        attrib = ["assembly_conv", "notconverted"]
        for a in attrib:
            for obj in getattr(self, a).values():
                if a == "assembly_conv":
                    id_noc = (
                        obj.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                    )
                    name = database_info[database_info["universal_bigg_id"] == id_noc][
                        "name"
                    ].values[0]
                else:
                    name = "Not converted"
                obj.name = name
                base_keys = list(obj.sources.keys()) + ["assembly", "comparison"]
                obj.reactions = {k: [] for k in base_keys}
                obj.formula = {k: [] for k in base_keys}


class SetofNewReactions(SetofNewObjects):
    """ Reactions class that add name, reaction equation and blank reactants/products attributes to reaction """

    def setReactionAttributes(self, database_info: pd.core.frame.DataFrame):
        attrib = ["assembly_conv", "notconverted"]
        for a in attrib:
            for obj in getattr(self, a).values():
                if a == "assembly_conv":
                    id_noc = obj.id.replace("sink_", "DM_")
                    name = database_info[database_info["bigg_id"] == id_noc]["name"]
                    if (not name.empty) and (not name.isnull().values.any()):
                        name = name.values[0]
                    else:
                        name = ""
                    equation = database_info[database_info["bigg_id"] == id_noc][
                        "reaction_string"
                    ]
                    if not equation.empty:
                        equation = equation.values[0]
                    else:
                        equation = None
                    obj.name = name
                else:
                    name = "Not converted"
                    equation = None
                obj.name = name
                obj.reaction = equation
                base_keys = list(obj.sources.keys()) + ["assembly", "comparison"]
                obj.reactants = {k: [] for k in base_keys}
                obj.products = {k: [] for k in base_keys}
                obj.metabolites = {k: {} for k in base_keys}
                obj.lower_bound = {k: [] for k in base_keys}
                obj.upper_bound = {k: [] for k in base_keys}
                obj.subsystem = {k: [] for k in base_keys}
                obj.genes = {k: [] for k in base_keys}
                obj.gene_reaction_rule = {k: [] for k in base_keys}
                obj.gene_reaction_rule.update(
                    {k + "_mixed": [] for k in obj.sources.keys()}
                )


class NewGene(object):
    """Class for one gene with new or old locus tag as ID and IDs from original models in annotation"""

    def __init__(self, new_id: str, old_id: str, source: str, possible_sources: [str]):
        self.id = new_id
        self.sources = {}
        self.in_models = {"models_amount": 1, "models_list": [source]}
        self.annotation = {}
        self.reactions = {"assembly": [], "comparison": []}
        for ps in possible_sources:
            self.reactions.update({ps: []})
            if ps == source:
                self.sources.update({ps: 1})
                self.annotation.update({ps: [old_id]})
            else:
                self.sources.update({ps: 0})
                self.annotation.update({ps: []})

    def updateNewGene(self, id_to_update: str, source: str):
        self.sources.update({source: self.sources.get(source) + 1})
        if source not in self.in_models["models_list"]:
            self.in_models["models_amount"] = self.in_models["models_amount"] + 1
            self.in_models["models_list"].append(source)
        self.annotation.get(source).append(id_to_update)


class SetofNewGenes(object):
    """ Setting dictionaries for all genes selected for supermodel - self.converted and not selected - self.notconverted. """

    def addNewGenes_conv(self, all_models_data: dict, gene_folder: PosixPath):
        for model_id in list(all_models_data.keys()):
            blast_file = gene_folder / (model_id + "_blast.tsv")
            try:
                conversion_table = pd.read_csv(str(blast_file), sep="\t", header=None)
            except:
                warnings.warn(
                    "\nWarning! File str(blast_file) can't be opened."
                    "\nOld gene will be used"
                )
                for gene in all_models_data[model_id]["preprocess_model"].genes:
                    if gene.id in self.assembly_conv.keys():
                        self.assembly_conv.get(gene.id).updateNewGene(gene.id, model_id)
                    else:
                        new_gene = NewGene(
                            gene.id, gene.id, model_id, list(all_models_data.keys())
                        )
                        self.assembly_conv.update({gene.id: new_gene})
                        self.comparison[model_id].update({gene.id: new_gene})
            else:
                conversion_table.columns = [
                    "old_id",
                    "new_id",
                    "identity",
                    "length",
                    "4",
                    "5",
                    "6",
                    "7",
                    "8",
                    "9",
                    "10",
                    "11",
                ]
                for gene in all_models_data[model_id]["preprocess_model"].genes:
                    old_gene_id = gene.id
                    attr = conversion_table[conversion_table["old_id"] == old_gene_id][
                        "new_id"
                    ]
                    if attr.empty:
                        if gene.id in self.notconverted.keys():
                            self.notconverted.get(gene.id).updateNewGene(
                                gene.id, model_id
                            )
                        else:
                            new_gene = NewGene(
                                gene.id, gene.id, model_id, list(all_models_data.keys())
                            )
                            self.notconverted.update({gene.id: new_gene})
                    elif type(attr.values[0]) != str:
                        if gene.id in self.notconverted.keys():
                            self.notconverted.get(gene.id).updateNewGene(
                                gene.id, model_id
                            )
                        else:
                            new_gene = NewGene(
                                gene.id, gene.id, model_id, list(all_models_data.keys())
                            )
                            self.notconverted.update({gene.id: new_gene})
                    else:
                        new_id = attr.values[0]
                        if new_id in self.assembly_conv.keys():
                            self.assembly_conv.get(new_id).updateNewGene(
                                gene.id, model_id
                            )
                        else:
                            new_gene = NewGene(
                                new_id, gene.id, model_id, list(all_models_data.keys())
                            )
                            self.assembly_conv.update({new_id: new_gene})
                            self.comparison[model_id].update({new_id: new_gene})

    def __init__(self, all_models_data: dict, gene_folder):
        self.assembly_conv = {}
        self.assembly_mix = {}
        self.comparison = defaultdict(dict)
        self.notconverted = {}
        if gene_folder is not None:
            self.addNewGenes_conv(all_models_data, gene_folder)
        else:
            for model_id in list(all_models_data.keys()):
                for gene in all_models_data.get(model_id).genes:
                    if gene.id in self.assembly_conv.keys():
                        self.assembly_conv.get(gene.id).updateNewGene(gene.id, model_id)
                    else:
                        new_gene = NewGene(
                            gene.id, gene.id, model_id, list(all_models_data.keys())
                        )
                        self.assembly_conv.update({gene.id: new_gene})
                        self.comparison[model_id].update({gene.id: new_gene})


class SuperModel:  # TODO REAL 30.08.23 add transport reactions for periplasmic metabolites for models without periplasmic compartments
    """ Supermodel class with metabolites and reactions. Sources - names of original models used to create supermodel.
    Creating connections between metabolites and reaction via dictionaries with sources as keys and links to
    reactants/products/reactions as values.  """

    def findReactions(
        self,
        metabolite: NewObject,
        m_go_new_old: dict,
        r_go_old_new: dict,
        model_ids: [str],
        periplasmic_r: dict,
        periplasmic_m: dict,
    ):
        for model_id in model_ids:
            old_mets = m_go_new_old.get(metabolite.id).get(model_id)
            if old_mets:
                new_r = []
                for old_met in old_mets:
                    if model_id in list(periplasmic_m.keys()):
                        if old_met.id in list(periplasmic_m.get(model_id).keys()):
                            for reaction in old_met.reactions:
                                if r_go_old_new.get(model_id).get(reaction.id):
                                    if (metabolite.id.endswith("_p")) & (
                                        reaction.id
                                        in list(periplasmic_r.get(model_id).keys())
                                    ):
                                        new_r.append(
                                            r_go_old_new.get(model_id).get(reaction.id)[
                                                0
                                            ]
                                        )
                                    elif (not metabolite.id.endswith("_p")) & (
                                        reaction.id
                                        not in list(periplasmic_r.get(model_id).keys())
                                    ):
                                        new_r.append(
                                            r_go_old_new.get(model_id).get(reaction.id)[
                                                0
                                            ]
                                        )
                        else:
                            for r in old_met.reactions:
                                if r_go_old_new.get(model_id).get(r.id):
                                    new_r.append(
                                        r_go_old_new.get(model_id).get(r.id)[0]
                                    )
                    else:
                        new_r = []
                        for r in old_met.reactions:
                            if r_go_old_new.get(model_id).get(r.id):
                                new_r.append(r_go_old_new.get(model_id).get(r.id)[0])
                if new_r:
                    metabolite.reactions[model_id] = list(set(new_r))

    def findMetabolites(
        self,
        reaction: NewObject,
        r_go_new_old: dict,
        m_go_old_new: dict,
        model_ids: [str],
        periplasmic_r: dict,
    ):
        for model_id in model_ids:
            old_react = r_go_new_old.get(reaction.id).get(model_id)
            if old_react:
                old_react_reactants = old_react[0].reactants
                old_react_products = old_react[0].products
                old_react_metabolites = old_react[0].metabolites
                if model_id in periplasmic_r.keys():
                    if old_react[0].id in periplasmic_r.get(model_id).keys():
                        for reactant in old_react_reactants:
                            new_reactants = m_go_old_new.get(model_id).get(reactant.id)
                            if new_reactants:
                                if len(new_reactants) == 1:
                                    reaction.reactants.get(model_id).append(
                                        new_reactants[0]
                                    )
                                elif len(new_reactants) > 1:
                                    for new_reactant in new_reactants:
                                        if (new_reactant.id.endswith("_p")) & (
                                            reactant.id
                                            in periplasmic_r.get(model_id)
                                            .get(old_react[0].id)
                                            .keys()
                                        ):
                                            reaction.reactants.get(model_id).append(
                                                new_reactant
                                            )
                                        if (not new_reactant.id.endswith("_p")) & (
                                            reactant.id
                                            not in periplasmic_r.get(model_id)
                                            .get(old_react[0].id)
                                            .keys()
                                        ):
                                            reaction.reactants.get(model_id).append(
                                                new_reactant
                                            )
                        for product in old_react_products:
                            new_products = m_go_old_new.get(model_id).get(product.id)
                            if new_products:
                                if len(new_products) == 1:
                                    reaction.products.get(model_id).append(
                                        new_products[0]
                                    )
                                elif len(new_products) > 1:
                                    for new_product in new_products:
                                        if (new_product.id.endswith("_p")) & (
                                            product.id
                                            in periplasmic_r.get(model_id)
                                            .get(old_react[0].id)
                                            .keys()
                                        ):
                                            reaction.products.get(model_id).append(
                                                new_product
                                            )
                                        if (not new_product.id.endswith("_p")) & (
                                            product.id
                                            not in periplasmic_r.get(model_id)
                                            .get(old_react[0].id)
                                            .keys()
                                        ):
                                            reaction.products.get(model_id).append(
                                                new_product
                                            )
                        for met, koef in old_react_metabolites.items():
                            new_mets = m_go_old_new.get(model_id).get(met.id)
                            if new_mets:
                                if len(new_mets) == 1:
                                    reaction.metabolites.get(model_id).update(
                                        {new_mets[0]: koef}
                                    )
                                elif len(new_mets) > 1:
                                    for new_met in new_mets:
                                        if (new_met.id.endswith("_p")) & (
                                            met.id
                                            in periplasmic_r.get(model_id)
                                            .get(old_react[0].id)
                                            .keys()
                                        ):
                                            reaction.metabolites.get(model_id).update(
                                                {new_met: koef}
                                            )
                                        if (not new_met.id.endswith("_p")) & (
                                            met.id
                                            not in periplasmic_r.get(model_id)
                                            .get(old_react[0].id)
                                            .keys()
                                        ):
                                            reaction.metabolites.get(model_id).update(
                                                {new_met: koef}
                                            )
                    else:
                        for reactant in old_react_reactants:
                            new_reactants = m_go_old_new.get(model_id).get(reactant.id)
                            if new_reactants:
                                if len(new_reactants) == 1:
                                    reaction.reactants.get(model_id).append(
                                        new_reactants[0]
                                    )
                                elif len(new_reactants) > 1:
                                    for new_reactant in new_reactants:
                                        if not new_reactant.id.endswith("_p"):
                                            reaction.reactants.get(model_id).append(
                                                new_reactant
                                            )
                        for product in old_react_products:
                            new_products = m_go_old_new.get(model_id).get(product.id)
                            if new_products:
                                if len(new_products) == 1:
                                    reaction.products.get(model_id).append(
                                        new_products[0]
                                    )
                                elif len(new_products) > 1:
                                    for new_product in new_products:
                                        if not new_product.id.endswith("_p"):
                                            reaction.products.get(model_id).append(
                                                new_product
                                            )
                        for met, koef in old_react_metabolites.items():
                            new_mets = m_go_old_new.get(model_id).get(met.id)
                            if new_mets:
                                if len(new_mets) == 1:
                                    reaction.metabolites.get(model_id).update(
                                        {new_mets[0]: koef}
                                    )
                                elif len(new_mets) > 1:
                                    for new_met in new_mets:
                                        if not new_met.id.endswith("_p"):
                                            reaction.metabolites.get(model_id).update(
                                                {new_met: koef}
                                            )
                else:
                    for reactant in old_react_reactants:
                        new_reactants = m_go_old_new.get(model_id).get(reactant.id)
                        if new_reactants:
                            reaction.reactants.get(model_id).append(new_reactants[0])
                    for product in old_react_products:
                        new_products = m_go_old_new.get(model_id).get(product.id)
                        if new_products:
                            reaction.products.get(model_id).append(new_products[0])
                    for met, koef in old_react_metabolites.items():
                        new_mets = m_go_old_new.get(model_id).get(met.id)
                        if new_mets:
                            reaction.metabolites.get(model_id).update(
                                {new_mets[0]: koef}
                            )

    def findGenes(
        self,
        all_models_data: dict,
        r_go_old_new: dict,
        r_go_new_old: dict,
        model_ids: [str],
        gene_folder: PosixPath,
    ):
        for model_id in model_ids:
            blast_file = gene_folder / (model_id + "_blast.tsv")
            conversion_table = pd.read_csv(str(blast_file), sep="\t", header=None)
            conversion_table.columns = [
                "old_id",
                "new_id",
                "identity",
                "length",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10",
                "11",
            ]
            for gene in self.genes.assembly_conv.values():
                if model_id in gene.in_models["models_list"]:
                    old_g_ids = gene.annotation.get(model_id)
                    for old_g_id in old_g_ids:
                        oldg_r_ids = [
                            gr.id
                            for gr in all_models_data[model_id]["preprocess_model"]
                            .genes.get_by_id(old_g_id)
                            .reactions
                        ]
                        for r_id in oldg_r_ids:
                            if r_go_old_new.get(model_id).get(r_id):
                                for new_r in r_go_old_new.get(model_id).get(r_id):
                                    if new_r not in gene.reactions.get(model_id):
                                        gene.reactions.get(model_id).append(new_r)
            for reaction in self.reactions.assembly_conv.values():
                old_rs = r_go_new_old.get(reaction.id).get(model_id)
                if old_rs:
                    new_gpr_unite_r = []
                    for oldr in old_rs:
                        if oldr.genes:
                            gene_convert = {}
                            for oldrg in oldr.genes:
                                attr_new = conversion_table[
                                    conversion_table["old_id"] == oldrg.id
                                ]["new_id"]
                                if not attr_new.empty:
                                    new_g_id = attr_new.values[0]
                                    if self.genes.assembly_conv.get(
                                        new_g_id
                                    ) not in reaction.genes.get(model_id):
                                        reaction.genes.get(model_id).append(
                                            self.genes.assembly_conv.get(new_g_id)
                                        )
                                    gene_convert.update({oldrg.id: new_g_id})
                                else:
                                    gene_convert.update({oldrg.id: "not_found"})
                            old_gpr = oldr.gene_reaction_rule
                            new_gpr, mix_gpr = makeNewGPR(old_gpr, gene_convert)
                            if new_gpr:
                                new_gpr_unite_r.append(new_gpr)
                            reaction.gene_reaction_rule.get(model_id + "_mixed").append(
                                mix_gpr
                            )
                    if len(new_gpr_unite_r) == 1:
                        reaction.gene_reaction_rule.get(model_id).append(
                            new_gpr_unite_r[0]
                        )
                    elif len(new_gpr_unite_r) >= 1:
                        united_gpr = uniteGPR(new_gpr_unite_r)
                        reaction.gene_reaction_rule.get(model_id).append(united_gpr)

    def findConnections(
        self,
        m_go_new_old: dict,
        m_go_old_new: dict,
        r_go_new_old: dict,
        r_go_old_new: dict,
        all_models_data: dict,
        periplasmic_r: dict,
        periplasmic_m: dict,
        gene_folder: PosixPath,
    ):
        model_ids = list(all_models_data.keys())
        for met in self.metabolites.assembly_conv.values():
            self.findReactions(
                met, m_go_new_old, r_go_old_new, model_ids, periplasmic_r, periplasmic_m
            )
        for r in self.reactions.assembly_conv.values():
            self.findMetabolites(
                r, r_go_new_old, m_go_old_new, model_ids, periplasmic_r
            )
        self.findGenes(
            all_models_data, r_go_old_new, r_go_new_old, model_ids, gene_folder
        )

    def getAdditionalAttributes(
        self, model_ids: [str], m_go_new_old: dict, r_go_new_old: dict
    ):
        for met in self.metabolites.assembly_conv.values():
            for model_id in model_ids:
                old_mets = m_go_new_old.get(met.id).get(model_id)
                if old_mets:
                    met.formula.get(model_id).append(old_mets[0].formula)
        for r in self.reactions.assembly_conv.values():
            for mod_id in model_ids:
                old_rs = r_go_new_old.get(r.id).get(mod_id)
                if old_rs:
                    low_b = 0
                    upp_b = 0
                    subsys = []
                    for old_r in old_rs:
                        if old_r.lower_bound < low_b:
                            low_b = old_r.lower_bound
                        if old_r.upper_bound > upp_b:
                            upp_b = old_r.upper_bound
                        subsys.append(old_r.subsystem)
                    r.lower_bound.get(mod_id).append(low_b)
                    r.upper_bound.get(mod_id).append(upp_b)
                    r.subsystem.get(mod_id).append("#or#".join(subsys))

    def __init__(
        self,
        metabolites: SetofNewMetabolites,
        reactions: SetofNewReactions,
        genes: SetofNewGenes,
        m_go_new_old: dict,
        m_go_old_new: dict,
        r_go_new_old: dict,
        r_go_old_new: dict,
        all_models_data: dict,
        periplasmic_r: dict,
        additional_periplasmic_m: dict,
        gene_folder,
    ):
        self.metabolites = metabolites
        self.reactions = reactions
        self.genes = genes
        self.sources = list(all_models_data.keys())
        self.findConnections(
            m_go_new_old,
            m_go_old_new,
            r_go_new_old,
            r_go_old_new,
            all_models_data,
            periplasmic_r,
            additional_periplasmic_m,
            gene_folder,
        )
        self.getAdditionalAttributes(self.sources, m_go_new_old, r_go_new_old)
        # supermodel.addBiomass(
        #     model_type, m_go_old_new, all_models, final_r_not_sel, final_m_not_sel
        # )

    # def addBiomass(
    #     self,
    #     types: [str],
    #     m_goOldNew: dict,
    #     all_models: dict,
    #     final_r_not_sel: dict,
    #     final_m_not_sel: dict,
    # ):
    #     new_biomass = None
    #     for typ in types:
    #         for r in all_models.get(typ).reactions:
    #             if len(r.reactants) > 24:
    #                 if not new_biomass:
    #                     new_biomass = NewObject(
    #                         "Biomass",
    #                         r.id,
    #                         final_r_not_sel.get(typ).get(r.id)[0],
    #                         typ,
    #                         types,
    #                         {},
    #                     )
    #                     new_biomass.name = "Biomass"
    #                     new_biomass.reaction = ""
    #                     new_biomass.reactants = {typ: []}
    #                     new_biomass.products = {typ: []}
    #                     new_biomass.metabolites = {typ: {}}
    #                     new_biomass.genes = {typ: []}
    #                     new_biomass.gene_reaction_rule = {typ: []}
    #                     new_biomass.lower_bound = {typ: [r.lower_bound]}
    #                     new_biomass.upper_bound = {typ: [r.upper_bound]}
    #                     new_biomass.subsystem = {typ: [r.subsystem]}
    #                     nc_biomass = NewObject(
    #                         "Biomass",
    #                         r.id,
    #                         final_r_not_sel.get(typ).get(r.id)[0],
    #                         typ,
    #                         types,
    #                         {},
    #                     )
    #                     nc_biomass.name = "Biomass"
    #                     nc_biomass.reaction = ""
    #                     nc_biomass.reactants = {typ: []}
    #                     nc_biomass.products = {typ: []}
    #                     nc_biomass.metabolites = {typ: {}}
    #                     nc_biomass.genes = {typ: []}
    #                     nc_biomass.gene_reaction_rule = {typ: []}
    #                     nc_biomass.lower_bound = {typ: [r.lower_bound]}
    #                     nc_biomass.upper_bound = {typ: [r.upper_bound]}
    #                     nc_biomass.subsystem = {typ: [r.subsystem]}
    #                 else:
    #                     new_biomass.updateNewObject(
    #                         r.id, final_r_not_sel.get(typ).get(r.id)[0], {}, typ
    #                     )
    #                     new_biomass.reactants.update({typ: []})
    #                     new_biomass.products.update({typ: []})
    #                     new_biomass.metabolites.update({typ: {}})
    #                     new_biomass.genes.update({typ: []})
    #                     new_biomass.gene_reaction_rule.update({typ: []})
    #                     new_biomass.lower_bound.update({typ: [r.lower_bound]})
    #                     new_biomass.upper_bound.update({typ: [r.upper_bound]})
    #                     new_biomass.subsystem.update({typ: [r.subsystem]})
    #                     nc_biomass.updateNewObject(
    #                         r.id, final_r_not_sel.get(typ).get(r.id)[0], {}, typ
    #                     )
    #                     nc_biomass.reactants.update({typ: []})
    #                     nc_biomass.products.update({typ: []})
    #                     nc_biomass.metabolites.update({typ: {}})
    #                     nc_biomass.genes.update({typ: []})
    #                     nc_biomass.gene_reaction_rule.update({typ: []})
    #                     nc_biomass.lower_bound.update({typ: [r.lower_bound]})
    #                     nc_biomass.upper_bound.update({typ: [r.upper_bound]})
    #                     nc_biomass.subsystem.update({typ: [r.subsystem]})
    #                 biomass_react = [mr.id for mr in r.reactants]
    #                 for reactant in biomass_react:
    #                     new_reactants = m_goOldNew.get(typ).get(reactant)
    #                     if new_reactants:
    #                         new_biomass.reactants.get(typ).append(new_reactants[0])
    #                         new_biomass.metabolites.get(typ).update(
    #                             {
    #                                 new_reactants[0]: r.metabolites.get(
    #                                     all_models.get(typ).metabolites.get_by_id(
    #                                         reactant
    #                                     )
    #                                 )
    #                             }
    #                         )
    #                         new_reactants[0].reactions.get(typ).append(new_biomass)
    #                     else:
    #                         if not final_m_not_sel.get(typ).get(reactant)[1]:
    #                             nc_biomass.reactants.get(typ).append(
    #                                 self.metabolites.notconverted.get(reactant)
    #                             )
    #                             nc_biomass.metabolites.get(typ).update(
    #                                 {
    #                                     self.metabolites.notconverted.get(
    #                                         reactant
    #                                     ): r.metabolites.get(
    #                                         all_models.get(typ).metabolites.get_by_id(
    #                                             reactant
    #                                         )
    #                                     )
    #                                 }
    #                             )
    #                             self.metabolites.notconverted.get(
    #                                 reactant
    #                             ).reactions.get(typ).append(nc_biomass)
    #                         else:
    #                             for bigg_reactant in final_m_not_sel.get(typ).get(
    #                                 reactant
    #                             )[1]:
    #                                 nc_biomass.reactants.get(typ).append(
    #                                     self.metabolites.notconverted.get(bigg_reactant)
    #                                 )
    #                                 nc_biomass.metabolites.get(typ).update(
    #                                     {
    #                                         self.metabolites.notconverted.get(
    #                                             bigg_reactant
    #                                         ): r.metabolites.get(
    #                                             all_models.get(
    #                                                 typ
    #                                             ).metabolites.get_by_id(reactant)
    #                                         )
    #                                     }
    #                                 )
    #                                 self.metabolites.notconverted.get(
    #                                     bigg_reactant
    #                                 ).reactions.get(typ).append(nc_biomass)
    #                 biomass_pro = [mp.id for mp in r.products]
    #                 for product in biomass_pro:
    #                     new_products = m_goOldNew.get(typ).get(product)
    #                     if new_products:
    #                         new_biomass.products.get(typ).append(new_products[0])
    #                         new_biomass.metabolites.get(typ).update(
    #                             {
    #                                 new_products[0]: r.metabolites.get(
    #                                     all_models.get(typ).metabolites.get_by_id(
    #                                         product
    #                                     )
    #                                 )
    #                             }
    #                         )
    #                         new_products[0].reactions.get(typ).append(new_biomass)
    #                     else:
    #                         if not final_m_not_sel.get(typ).get(product)[1]:
    #                             nc_biomass.products.get(typ).append(
    #                                 self.metabolites.notconverted.get(product)
    #                             )
    #                             nc_biomass.metabolites.get(typ).update(
    #                                 {
    #                                     self.metabolites.notconverted.get(
    #                                         product
    #                                     ): r.metabolites.get(
    #                                         all_models.get(typ).metabolites.get_by_id(
    #                                             product
    #                                         )
    #                                     )
    #                                 }
    #                             )
    #                             self.metabolites.notconverted.get(
    #                                 product
    #                             ).reactions.get(typ).append(nc_biomass)
    #                         else:
    #                             for bigg_product in final_m_not_sel.get(typ).get(
    #                                 product
    #                             )[1]:
    #                                 nc_biomass.products.get(typ).append(
    #                                     self.metabolites.notconverted.get(bigg_product)
    #                                 )
    #                                 nc_biomass.metabolites.get(typ).update(
    #                                     {
    #                                         self.metabolites.notconverted.get(
    #                                             bigg_product
    #                                         ): r.metabolites.get(
    #                                             all_models.get(
    #                                                 typ
    #                                             ).metabolites.get_by_id(product)
    #                                         )
    #                                     }
    #                                 )
    #                                 self.metabolites.notconverted.get(
    #                                     bigg_product
    #                                 ).reactions.get(typ).append(nc_biomass)
    #     self.reactions.converted.update({"Biomass": new_biomass})
    #     self.reactions.notconverted.update({"Biomass": nc_biomass})
