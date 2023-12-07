import operator
import sys
import warnings
from collections import defaultdict
from math import ceil
from os.path import exists
from pathlib import PosixPath

import dill
from future.moves import itertools
from .comparison import (
    getCoreConnections,
    getCoreGPR,
    getCoreLowerBounds,
    getCoreCoefficients,
    getCoreUpperBounds,
    getCore,
    getDifference,
)
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

    def _updateNewObject(
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
    selected for supermodel - self.assembly and not selected - self.notconverted. """

    def __addNewObjs(self, selected: dict, where_to_add: str, model_ids: list):
        for mod_id in model_ids:
            if mod_id in list(selected.keys()):
                objects = selected.get(mod_id)
                for key in objects.keys():
                    for new_id in objects[key][1]:
                        comp = objects[key][0]
                        if new_id in getattr(self, where_to_add).keys():
                            getattr(self, where_to_add).get(new_id)._updateNewObject(
                                key, comp, mod_id
                            )
                        else:
                            new = NewObject(new_id, key, comp, mod_id, model_ids)
                            getattr(self, where_to_add).update({new_id: new})

    def __makeSetofNew(
        self, selected: dict, not_selected: dict, model_ids, additional,
    ):
        self.__addNewObjs(selected, "assembly", model_ids)
        if additional:
            self.__addNewObjs(additional, "assembly", model_ids)
        for new_id, new_obj in self.assembly.items():
            for model_id in new_obj.in_models["models_list"]:
                getattr(self, model_id).update({new_id: new_obj})
        self.__addNewObjs(
            not_selected, "notconverted", model_ids
        )  # TODO connect not_converted for really not converted only with old id

    def __init__(
        self, selected: dict, not_selected: dict, model_ids: [str], additional=None,
    ):
        self.assembly = {}
        for source in selected.keys():
            setattr(self, source, {})
        self.assembly_mix = {}
        self.comparison = defaultdict(dict)
        self.notconverted = {}
        self.__makeSetofNew(selected, not_selected, model_ids, additional)

    def _makeForwardBackward(
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
                new_obj = [self.assembly.get(value[1][0])]
                if additional:
                    if model_id in list(additional.keys()):
                        if key in list(additional.get(model_id).keys()):
                            new_obj.append(
                                self.assembly.get(
                                    additional.get(model_id).get(key)[1][0]
                                )
                            )
                go_old_new[model_id].update({key: new_obj})
        for k, v in self.assembly.items():
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

    def _setMetaboliteAttributes(self, database_info: pd.core.frame.DataFrame):
        attrib = ["assembly", "notconverted"]
        for a in attrib:
            for obj in getattr(self, a).values():
                if a == "assembly":
                    id_noc = (
                        obj.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                    )
                    name = database_info[database_info["universal_bigg_id"] == id_noc][
                        "name"
                    ].values[0]
                else:
                    name = "Not converted"
                obj.name = name
                obj.reactions = {k: [] for k in obj.sources.keys()}
                obj.reactions.update({"assembly": [], "comparison": {}})
                obj.formula = {k: [] for k in obj.sources.keys()}


class SetofNewReactions(SetofNewObjects):
    """ Reactions class that add name, reaction equation and blank reactants/products attributes to reaction """

    def _setReactionAttributes(self, database_info: pd.core.frame.DataFrame):
        attrib = ["assembly", "notconverted"]
        for a in attrib:
            for obj in getattr(self, a).values():
                if a == "assembly":
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
                base_keys = list(obj.sources.keys()) + ["assembly"]
                obj.reactants = {k: [] for k in base_keys}
                obj.reactants.update({"comparison": {}})
                obj.products = {k: [] for k in base_keys}
                obj.products.update({"comparison": {}})
                obj.metabolites = {k: {} for k in base_keys + ["comparison"]}
                obj.lower_bound = {k: [] for k in base_keys}
                obj.lower_bound.update({"comparison": {}})
                obj.upper_bound = {k: [] for k in base_keys}
                obj.upper_bound.update({"comparison": {}})
                obj.subsystem = {k: [] for k in list(obj.sources.keys())}
                obj.genes = {k: [] for k in base_keys}
                obj.genes.update({"comparison": {}})
                obj.gene_reaction_rule = {k: [] for k in base_keys}
                obj.gene_reaction_rule.update(
                    {k + "_mixed": [] for k in obj.sources.keys()}
                )
                obj.gene_reaction_rule.update({"comparison": {}})


class NewGene(object):
    """Class for one gene with new or old locus tag as ID and IDs from original models in annotation"""

    def __init__(self, new_id: str, old_id: str, source: str, possible_sources: [str]):
        self.id = new_id
        self.sources = {}
        self.in_models = {"models_amount": 1, "models_list": [source]}
        self.annotation = {}
        self.reactions = {"assembly": [], "comparison": {}}
        for ps in possible_sources:
            self.reactions.update({ps: []})
            if ps == source:
                self.sources.update({ps: 1})
                self.annotation.update({ps: [old_id]})
            else:
                self.sources.update({ps: 0})
                self.annotation.update({ps: []})

    def _updateNewGene(self, id_to_update: str, source: str):
        self.sources.update({source: self.sources.get(source) + 1})
        if source not in self.in_models["models_list"]:
            self.in_models["models_amount"] = self.in_models["models_amount"] + 1
            self.in_models["models_list"].append(source)
        self.annotation.get(source).append(id_to_update)


class SetofNewGenes(object):
    """ Setting dictionaries for all genes selected for supermodel - self.converted and not selected - self.notconverted. """

    def __addNewGenes_conv(self, all_models_data: dict, gene_folder: PosixPath):
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
                    if gene.id in self.assembly.keys():
                        self.assembly.get(gene.id)._updateNewGene(gene.id, model_id)
                    else:
                        new_gene = NewGene(
                            gene.id, gene.id, model_id, list(all_models_data.keys())
                        )
                        self.assembly.update({gene.id: new_gene})
                        getattr(self, model_id).update({gene.id: new_gene})
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
                        if new_id in self.assembly.keys():
                            self.assembly.get(new_id)._updateNewGene(gene.id, model_id)
                        else:
                            new_gene = NewGene(
                                new_id, gene.id, model_id, list(all_models_data.keys())
                            )
                            self.assembly.update({new_id: new_gene})
                            getattr(self, model_id).update({new_id: new_gene})

    def __init__(self, all_models_data: dict, gene_folder):
        self.assembly = {}
        for source in list(all_models_data.keys()):
            setattr(self, source, {})
        self.assembly_mix = {}
        self.comparison = defaultdict(dict)
        self.notconverted = {}
        if gene_folder is not None:
            self.__addNewGenes_conv(all_models_data, gene_folder)
        else:
            for model_id in list(all_models_data.keys()):
                for gene in all_models_data.get(model_id).genes:
                    if gene.id in self.assembly.keys():
                        self.assembly.get(gene.id)._updateNewGene(gene.id, model_id)
                    else:
                        new_gene = NewGene(
                            gene.id, gene.id, model_id, list(all_models_data.keys())
                        )
                        self.assembly.update({gene.id: new_gene})
                        getattr(self, model_id).update({gene.id: new_gene})


class SuperModel:  # TODO REAL 30.08.23 add transport reactions for periplasmic metabolites for models without periplasmic compartments
    """ Supermodel class with metabolites and reactions. Sources - names of original models used to create supermodel.
    Creating connections between metabolites and reaction via dictionaries with sources as keys and links to
    reactants/products/reactions as values.  """

    def __find_reactions(
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

    def __find_metabolites(
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

    def __find_genes(
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
            for gene in self.genes.assembly.values():
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
            for reaction in self.reactions.assembly.values():
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
                                    if self.genes.assembly.get(
                                        new_g_id
                                    ) not in reaction.genes.get(model_id):
                                        reaction.genes.get(model_id).append(
                                            self.genes.assembly.get(new_g_id)
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

    def __find_connections(
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
        for met in self.metabolites.assembly.values():
            self.__find_reactions(
                met, m_go_new_old, r_go_old_new, model_ids, periplasmic_r, periplasmic_m
            )
        for r in self.reactions.assembly.values():
            self.__find_metabolites(
                r, r_go_new_old, m_go_old_new, model_ids, periplasmic_r
            )
        self.__find_genes(
            all_models_data, r_go_old_new, r_go_new_old, model_ids, gene_folder
        )

    def __get_additional_attributes(
        self, model_ids: [str], m_go_new_old: dict, r_go_new_old: dict
    ):
        for met in self.metabolites.assembly.values():
            for model_id in model_ids:
                old_mets = m_go_new_old.get(met.id).get(model_id)
                if old_mets:
                    met.formula.get(model_id).append(old_mets[0].formula)
        for r in self.reactions.assembly.values():
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

    def __swapReactantsAndProducts(self, r: NewObject, sources_to_swap: list):
        for s in sources_to_swap:
            a = r.reactants.get(s)
            b = r.products.get(s)
            r.reactants[s] = b
            r.products[s] = a
            aa = r.lower_bound.get(s)[0] * -1
            bb = r.upper_bound.get(s)[0] * -1
            r.lower_bound[s] = [bb]
            r.upper_bound[s] = [aa]
            for met, koef in r.metabolites.get(s).items():
                r.metabolites.get(s)[met] = koef * -1

    def __runSwitchedMetabolites(self):
        for r in self.reactions.assembly.values():
            ex = False
            react_in = r.reactants[r.in_models["models_list"][0]]
            pro_in = r.products[r.in_models["models_list"][0]]
            for tmp in r.in_models["models_list"]:
                react_in = list(set(react_in) & set(r.reactants.get(tmp)))
                pro_in = list(set(pro_in) & set(r.products.get(tmp)))
                if (not r.reactants.get(tmp)) | (not r.products.get(tmp)):
                    ex = True
            if not ex:
                if (not react_in) | (not pro_in):
                    up = r.in_models["models_amount"] - 1
                    down = ceil(r.in_models["models_amount"] / 2) - 1
                    consist = []
                    for i in range(up, down, -1):
                        combinations = list(
                            itertools.combinations(r.in_models["models_list"], i)
                        )
                        for comb in combinations:
                            react_in_comb = r.reactants[comb[0]]
                            for c in comb:
                                react_in_comb = list(
                                    set(react_in_comb) & set(r.reactants.get(c))
                                )
                            if react_in_comb:
                                consist.append(comb)
                        if consist != []:
                            break
                    if len(consist) == 1:
                        # "Case 1: majority"
                        source_to_swap = list(
                            set(r.in_models["models_list"]) - set(consist[0])
                        )
                        self.__swapReactantsAndProducts(r, source_to_swap)
                    elif len(consist) == 2:
                        lb1 = 0
                        lb2 = 0
                        for tmp in r.in_models["models_list"]:
                            if tmp in consist[0]:
                                if r.lower_bound.get(tmp)[0] < lb1:
                                    lb1 = r.lower_bound.get(tmp)[0]
                            if tmp in consist[1]:
                                if r.lower_bound.get(tmp)[0] < lb2:
                                    lb2 = r.lower_bound.get(tmp)[0]
                        swap = None
                        if (lb1 >= 0) & (lb2 < 0):
                            swap = consist[1]
                        if (lb1 < 0) & (lb2 >= 0):
                            swap = consist[0]
                        if swap:
                            # "Case 2: boundary"
                            self.__swapReactantsAndProducts(r, swap)
                        else:
                            # "Case 3: Nothing sort"
                            sel = sorted(r.in_models["models_list"])[0]
                            not_sel = []
                            for tmp in sorted(r.in_models["models_list"])[1:]:
                                if not (
                                    set(r.reactants.get(tmp))
                                    & set(r.reactants.get(sel))
                                ):
                                    not_sel.append(tmp)
                            self.__swapReactantsAndProducts(r, not_sel)
                    # Maybe remove, since len(consist) is expected to be only 1 or 2
                    else:
                        print("Can enter consist more 2")
                        # "Case 3: Nothing sort"
                        sel = sorted(r.in_models["models_list"])[0]
                        not_sel = []
                        for tmp in sorted(r.in_models["models_list"])[1:]:
                            if not (
                                set(r.reactants.get(tmp)) & set(r.reactants.get(sel))
                            ):
                                not_sel.append(tmp)
                        self.__swapReactantsAndProducts(r, not_sel)

    def __assemble_attributes(self, and_as_solid: bool):
        for met in self.metabolites.assembly.values():
            ass_r = getCoreConnections(met.reactions, 1, operator.ge, self.sources)
            met.reactions.update({"assembly": ass_r})
        for gene in self.genes.assembly.values():
            ass_rg = getCoreConnections(gene.reactions, 1, operator.ge, self.sources)
            gene.reactions.update({"assembly": ass_rg})
        for react in self.reactions.assembly.values():
            ass_reactants = getCoreConnections(
                react.reactants, 1, operator.ge, self.sources
            )
            ass_products = getCoreConnections(
                react.products, 1, operator.ge, self.sources
            )
            ass_genes = getCoreConnections(react.genes, 1, operator.ge, self.sources)
            ass_gpr = getCoreGPR(
                react.gene_reaction_rule, 1, operator.ge, self.sources, and_as_solid,
            )
            ass_lower_bound = getCoreLowerBounds(
                react.lower_bound, 1, react.in_models["models_list"]
            )
            ass_upper_bound = getCoreUpperBounds(
                react.upper_bound, 1, react.in_models["models_list"]
            )
            react.reactants.update({"assembly": ass_reactants})
            react.products.update({"assembly": ass_products})
            react.genes.update({"assembly": ass_genes})
            react.gene_reaction_rule.update({"assembly": ass_gpr})
            react.lower_bound.update({"assembly": ass_lower_bound})
            react.upper_bound.update({"assembly": ass_upper_bound})
            core_metabolites = getCoreCoefficients(
                react.metabolites,
                react.reactants,
                react.products,
                "assembly",
                1,
                react.in_models["models_list"],
            )
            react.metabolites.update({"assembly": core_metabolites})

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
        and_as_solid: bool,
    ):
        self.metabolites = metabolites
        self.reactions = reactions
        self.genes = genes
        self.sources = list(all_models_data.keys())
        self.__find_connections(
            m_go_new_old,
            m_go_old_new,
            r_go_new_old,
            r_go_old_new,
            all_models_data,
            periplasmic_r,
            additional_periplasmic_m,
            gene_folder,
        )
        self.__get_additional_attributes(self.sources, m_go_new_old, r_go_new_old)
        self.__runSwitchedMetabolites()
        self.__assemble_attributes(and_as_solid)

    def get_short_name_len(self) -> int:
        for i in range(len(max(self.sources, key=len)) + 1):
            short = []
            for source in self.sources:
                short.append(source[:i])
            if len(set(short)) == len(self.sources):
                return i

    def at_least_in(self, number_of_model: int, and_as_solid=False):
        if (
            type(number_of_model) != int
            or number_of_model < 1
            or number_of_model > len(self.sources)
        ):
            raise ValueError("Number to check does not fit the number of models")
        elif number_of_model == 1:
            raise ValueError(
                "Features in at least 1 model are already found in assembly. "
                "You do not need to run this comparison separately"
            )
        else:
            getCore(self, number_of_model, operator.ge, and_as_solid)

    def exactly_in(self, number_of_model: int, and_as_solid=False):
        if (
            type(number_of_model) != int
            or number_of_model < 1
            or number_of_model > len(self.sources)
        ):
            raise ValueError("Number to check does not fit the number of models")
        else:
            getCore(self, number_of_model, operator.eq, and_as_solid)

    def present(self, yes=None, no=None, short_name_len=None, and_as_solid=False):
        if yes is None and no is None:
            raise ValueError(
                "Both models present and models not present are not provided. "
                "Please provide at least one of the list"
            )
        elif (yes is not None and type(yes) != list) or (
            no is not None and type(no) != list
        ):
            raise ValueError(
                "Present or not present models are in wrong type. "
                "Please provide lists"
            )
        else:
            if yes is None:
                yes = []
            if no is None:
                no = []
            wrong_yes = set(yes) - set(self.sources)
            wrong_no = set(no) - set(self.sources)
            if wrong_yes or wrong_no:
                raise ValueError(
                    "Some of input models are not in supermodel. "
                    "Please check the input ids"
                )
            else:
                if short_name_len is None:
                    short_name_len = self.get_short_name_len()
                getDifference(self, yes, no, and_as_solid, short_name_len)

    def get_venn_segments(self, short_name_len=None, and_as_solid=False):
        """ Getting metabolites and reactions networks for each Venn segment in Venn
        diagram."""
        if short_name_len is None:
            short_name_len = self.get_short_name_len()
        combinations = []
        for i in range(1, len(self.sources)):
            combinations.extend(itertools.combinations(self.sources, i))
        for combo in combinations:
            yes = sorted(list(combo))
            no = sorted((list(set(self.sources) - set(combo))))
            getDifference(self, yes, no, and_as_solid, short_name_len)

    def get_intersection(self, and_as_solid=False):
        getCore(self, len(self.sources), operator.ge, and_as_solid)

    def write_supermodel_to_pkl(self, output_name: str, recursion_limit=None):
        if not output_name.endswith(".pkl"):
            raise ValueError("Wrong extension of the file")
        if exists(output_name):
            raise ValueError("File already exist, change the name")
        else:
            if recursion_limit is None:
                recursion_limit = 50000
            sys.setrecursionlimit(recursion_limit)
            dill.dump(self, open(output_name, "wb"))


def read_supermodel_from_pkl(input_name: str):
    supermodel = dill.load(open(input_name, "rb"))
    return supermodel
