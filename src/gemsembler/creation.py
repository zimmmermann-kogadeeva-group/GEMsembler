import itertools
import operator
import re
import resource
import sys
import warnings
from collections import defaultdict
from math import ceil
from os.path import exists
from pathlib import PosixPath

import dill
import pandas as pd

from .comparison import (
    getCore,
    getCoreCoefficients,
    getCoreConnections,
    getCoreGPR,
    getCoreLowerBounds,
    getCoreUpperBounds,
    getDifference,
)
from .genes import makeNewGPR, uniteGPR


class NewElement:
    """ New object class - one metabolite or reaction for supermodel. """

    def __init__(
        self,
        new_id: str,
        old_id: str,
        compartments: [str],
        source: str,
        possible_sources: [str],
        converted: bool,
    ):
        self.id = new_id
        self.compartments = {"assembly": compartments}
        self.sources = {}
        self.in_models = {"models_amount": 1, "models_list": [source]}
        self.annotation = {}
        self.converted = converted
        for ps in possible_sources:
            if ps == source:
                self.compartments.update({ps: compartments})
                self.sources.update({ps: 1})
                self.annotation.update({ps: [old_id]})
            else:
                self.compartments.update({ps: []})
                self.sources.update({ps: 0})
                self.annotation.update({ps: []})

    def _update_new_element(
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


class NewMetabolite(NewElement):
    def __init__(
        self,
        new_id: str,
        old_id: str,
        compartments: [str],
        source: str,
        possible_sources: [str],
        converted: bool,
        m_database_info: pd.core.frame.DataFrame,
    ):
        super().__init__(
            new_id, old_id, compartments, source, possible_sources, converted
        )
        if converted:
            id_noc = re.sub("_([cep])$", "", new_id)
            name = m_database_info[m_database_info["universal_bigg_id"] == id_noc][
                "name"
            ].values[0]
        else:
            name = "Not converted"
        self.name = name
        self.reactions = {k: [] for k in possible_sources}
        self.reactions.update({"assembly": [], "comparison": {}})
        self.formula = {k: [] for k in possible_sources}

    def _find_reactions(
        self,
        m_go_new_old: dict,
        r_go_old_new: dict,
        periplasmic_r: dict,
        periplasmic_m: dict,
    ):
        for model_id in self.sources.keys():
            old_mets = m_go_new_old.get(self.id).get(model_id)
            if not old_mets:
                continue
            new_r = []
            for old_met in old_mets:
                if old_met.id in list(periplasmic_m.get(model_id, {}).keys()):
                    # Old metabolite has additional periplasmic version in supermodel.
                    # So we need to split its reactions between them
                    for reaction in old_met.reactions:
                        if not r_go_old_new.get(model_id).get(reaction.id):
                            continue
                        met_is_periplasmic = self.id.endswith("_p")
                        met_was_converted_to_periplasmic = (
                            old_met.id
                            in periplasmic_r.get(model_id).get(reaction.id, {}).keys()
                        )
                        if (met_is_periplasmic & met_was_converted_to_periplasmic) or (
                            (not met_is_periplasmic)
                            & (not met_was_converted_to_periplasmic)
                        ):
                            new_r.append(r_go_old_new.get(model_id).get(reaction.id)[0])
                else:
                    # No need to split reactions
                    for r in old_met.reactions:
                        if r_go_old_new.get(model_id).get(r.id):
                            new_r.append(r_go_old_new.get(model_id).get(r.id)[0])
            if new_r:
                self.reactions[model_id] = list(set(new_r))


class NewReaction(NewElement):
    def __init__(
        self,
        new_id: str,
        old_id: str,
        compartments: [str],
        source: str,
        possible_sources: [str],
        converted: bool,
        r_database_info: pd.core.frame.DataFrame,
    ):
        super().__init__(
            new_id, old_id, compartments, source, possible_sources, converted
        )
        if converted:
            id_noc = new_id.replace("sink_", "DM_")
            name = r_database_info[r_database_info["bigg_id"] == id_noc]["name"]
            if (not name.empty) and (not name.isnull().values.any()):
                name = name.values[0]
            else:
                name = ""
            equation = r_database_info[r_database_info["bigg_id"] == id_noc][
                "reaction_string"
            ]
            if not equation.empty:
                equation = equation.values[0]
            else:
                equation = None
        else:
            name = "Not converted"
            equation = None
        self.name = name
        self.reaction = equation
        base_keys = possible_sources + ["assembly"]
        self.reactants = {k: [] for k in base_keys}
        self.reactants.update({"comparison": {}})
        self.products = {k: [] for k in base_keys}
        self.products.update({"comparison": {}})
        self.metabolites = {k: {} for k in base_keys + ["comparison"]}
        self.lower_bound = {k: [] for k in base_keys}
        self.lower_bound.update({"comparison": {}})
        self.upper_bound = {k: [] for k in base_keys}
        self.upper_bound.update({"comparison": {}})
        self.subsystem = {k: [] for k in possible_sources}
        self.genes = {k: [] for k in base_keys}
        self.genes.update({"comparison": {}})
        self.gene_reaction_rule = {k: [] for k in base_keys}
        self.gene_reaction_rule.update({k + "_mixed": [] for k in possible_sources})
        self.gene_reaction_rule.update({"comparison": {}})

    def __sel_met_from_p_model_for_p_r(
        self,
        old_react_metabolites: list,
        old_react_id: str,
        model_id: str,
        m_go_old_new: dict,
        periplasmic_r: dict,
    ):
        out_met = {}
        for met in old_react_metabolites:
            new_mets = m_go_old_new.get(model_id).get(met.id)
            if not new_mets:
                continue
            if len(new_mets) == 1:
                out_met.update({new_mets[0]: met})
            else:
                test_of_single_entry = True
                for new_met in new_mets:
                    new_met_is_periplasmic = new_met.id.endswith("_p")
                    new_met_was_converted_to_periplasmic = (
                        met.id in periplasmic_r.get(model_id).get(old_react_id).keys()
                    )
                    if (
                        new_met_is_periplasmic & new_met_was_converted_to_periplasmic
                    ) or (
                        (not new_met_is_periplasmic)
                        & (not new_met_was_converted_to_periplasmic)
                    ):
                        out_met.update({new_met: met})
                        if test_of_single_entry:
                            test_of_single_entry = False
                            continue
                        if not test_of_single_entry:
                            problem_m = [n.id for n, m in out_met.items() if m == met]
                            warnings.warn(
                                f"Something went wrong with periplasmic"
                                f"connections of {self.id}. Problematic"
                                f"metabolites are {' '.join(problem_m)}"
                            )
        return out_met

    def __sel_met_from_p_model_for_not_p_r(
        self, old_react_metabolites: list, model_id: str, m_go_old_new: dict
    ):
        out_met = {}
        for met in old_react_metabolites:
            new_mets = m_go_old_new.get(model_id).get(met.id)
            if not new_mets:
                continue
            if len(new_mets) == 1:
                out_met.update({new_mets[0]: met})
            else:
                test_of_single_entry = True
                for new_met in new_mets:
                    if not new_met.id.endswith("_p"):
                        out_met.update({new_met: met})
                        if test_of_single_entry:
                            test_of_single_entry = False
                            continue
                        if not test_of_single_entry:
                            problem_m = [n.id for n, m in out_met.items() if m == met]
                            warnings.warn(
                                f"Something went wrong with periplasmic"
                                f"connections of {self.id}. Problematic"
                                f"metabolites are {' '.join(problem_m)}"
                            )
        return out_met

    def _find_reactants_products(
        self, r_go_new_old: dict, m_go_old_new: dict, periplasmic_r: dict, m_type: str
    ):
        for model_id in self.sources.keys():
            old_react = r_go_new_old.get(self.id).get(model_id)
            # old_react is list, usually with 1 element,
            # but even if not 1, these reactions have the same r equation,
            # so we can take any and I tool the 1st
            if not old_react:
                continue
            old_react_react_prod = getattr(old_react[0], m_type)
            model_has_periplasmic_changes = model_id in periplasmic_r.keys()
            reaction_has_periplasmic_changes = (
                old_react[0].id in periplasmic_r.get(model_id, {}).keys()
            )
            if (model_has_periplasmic_changes) & (reaction_has_periplasmic_changes):
                met_for_p_r = self.__sel_met_from_p_model_for_p_r(
                    old_react_react_prod,
                    old_react[0].id,
                    model_id,
                    m_go_old_new,
                    periplasmic_r,
                )
                for m in met_for_p_r.keys():
                    getattr(self, m_type).get(model_id).append(m)
            elif model_has_periplasmic_changes & (not reaction_has_periplasmic_changes):
                met_for_not_p_r = self.__sel_met_from_p_model_for_not_p_r(
                    old_react_react_prod, model_id, m_go_old_new,
                )
                for m in met_for_not_p_r.keys():
                    getattr(self, m_type).get(model_id).append(m)
            else:
                # There was no periplasmic perturbation in the model
                # Only 1 element in new_reacts_prods is expected
                for react_prod in old_react_react_prod:
                    new_reacts_prods = m_go_old_new.get(model_id).get(react_prod.id)
                    if new_reacts_prods:
                        getattr(self, m_type).get(model_id).append(new_reacts_prods[0])
                        if len(new_reacts_prods) > 1:
                            warnings.warn(
                                f"Unexpected not unique connections between new "
                                f"and old metabolite without periplasmic story."
                                f"Model: {model_id}. Old metabolite: {react_prod.id}."
                                f"New metabolites: "
                                f"{' '.join([n.id for n in new_reacts_prods])}"
                            )

    def _find_metabolites(
        self, r_go_new_old: dict, m_go_old_new: dict, periplasmic_r: dict,
    ):
        for model_id in self.sources.keys():
            old_react = r_go_new_old.get(self.id).get(model_id)
            # old_react is list, usually with 1 element,
            # but even if not 1, these reactions have the same r equation,
            # so we can take any and I tool the 1st
            if not old_react:
                continue
            old_react_metabolites = old_react[0].metabolites
            model_has_periplasmic_changes = model_id in periplasmic_r.keys()
            reaction_has_periplasmic_changes = (
                old_react[0].id in periplasmic_r.get(model_id, {}).keys()
            )
            if model_has_periplasmic_changes & reaction_has_periplasmic_changes:
                met_for_p_r = self.__sel_met_from_p_model_for_p_r(
                    list(old_react_metabolites.keys()),
                    old_react[0].id,
                    model_id,
                    m_go_old_new,
                    periplasmic_r,
                )
                for m, v in met_for_p_r.items():
                    self.metabolites.get(model_id).update({m: old_react_metabolites[v]})
            elif model_has_periplasmic_changes & (not reaction_has_periplasmic_changes):
                met_for_not_p_r = self.__sel_met_from_p_model_for_not_p_r(
                    list(old_react_metabolites.keys()), model_id, m_go_old_new,
                )
                for m, v in met_for_not_p_r.items():
                    self.metabolites.get(model_id).update({m: old_react_metabolites[v]})
            else:
                # There was no periplasmic perturbation in the model
                # Only 1 element in new_mets is expected
                for met, koef in old_react_metabolites.items():
                    new_mets = m_go_old_new.get(model_id).get(met.id)
                    if new_mets:
                        self.metabolites.get(model_id).update({new_mets[0]: koef})
                        if len(new_mets) > 1:
                            warnings.warn(
                                f"Unexpected not unique connections between new "
                                f"and old metabolite without periplasmic story."
                                f"Model: {model_id}. Old metabolite: {met.id}."
                                f"New metabolites: {' '.join([n.id for n in new_mets])}"
                            )


class SetofNewElements:
    """ Setting dictionaries for all metabolites or reactions:
    selected for supermodel - self.assembly and not selected - self.notconverted. """

    def __add_new_elements(
        self,
        element_type: str,
        selected: dict,
        where_to_add: str,
        model_ids: list,
        convered: bool,
        db_info: pd.core.frame.DataFrame,
    ):
        new_elements = {"metabolites": NewMetabolite, "reactions": NewReaction}
        for mod_id in model_ids:
            if mod_id not in list(selected.keys()):
                continue
            objects = selected.get(mod_id)
            for key in objects.keys():
                for new_id in objects[key][1]:
                    comp = objects[key][0]
                    if new_id in getattr(self, where_to_add).keys():
                        getattr(self, where_to_add).get(new_id)._update_new_element(
                            key, comp, mod_id
                        )
                    else:
                        new = new_elements[element_type](
                            new_id, key, comp, mod_id, model_ids, convered, db_info
                        )
                        getattr(self, where_to_add).update({new_id: new})

    def __init__(
        self,
        element_type: str,
        selected: dict,
        not_selected: dict,
        model_ids: [str],
        db_info: pd.core.frame.DataFrame,
        additional=None,
    ):
        self.assembly = {}
        for source in selected.keys():
            setattr(self, source, {})
        self.assembly_mix = {}
        self.comparison = defaultdict(dict)
        self.notconverted = {}
        self.__add_new_elements(
            element_type, selected, "assembly", model_ids, True, db_info
        )
        if additional:
            self.__add_new_elements(
                element_type, additional, "assembly", model_ids, True, db_info
            )
        for new_id, new_obj in self.assembly.items():
            for model_id in new_obj.in_models["models_list"]:
                getattr(self, model_id).update({new_id: new_obj})
        self.__add_new_elements(
            element_type, not_selected, "notconverted", model_ids, False, db_info
        )

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
                    if key in list(additional.get(model_id, {}).keys()):
                        new_obj.append(
                            self.assembly.get(additional.get(model_id).get(key)[1][0])
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
                    f"\nWarning! File {str(blast_file)} can't be opened."
                    f"\nOld gene will be used"
                )
                for gene in all_models_data[model_id]["preprocess_model"].genes:
                    if gene.id in self.assembly.keys():
                        self.assembly.get(gene.id)._updateNewGene(gene.id, model_id)
                        getattr(self, model_id).update(
                            {gene.id: self.assembly.get(gene.id)}
                        )
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
                            self.notconverted.get(gene.id)._updateNewGene(
                                gene.id, model_id
                            )
                        else:
                            new_gene = NewGene(
                                gene.id, gene.id, model_id, list(all_models_data.keys())
                            )
                            self.notconverted.update({gene.id: new_gene})
                    elif type(attr.values[0]) != str:
                        if gene.id in self.notconverted.keys():
                            self.notconverted.get(gene.id)._updateNewGene(
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
                            getattr(self, model_id).update(
                                {new_id: self.assembly.get(new_id)}
                            )
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
                for gene in all_models_data[model_id]["preprocess_model"].genes:
                    if gene.id in self.assembly.keys():
                        self.assembly.get(gene.id)._updateNewGene(gene.id, model_id)
                    else:
                        new_gene = NewGene(
                            gene.id, gene.id, model_id, list(all_models_data.keys())
                        )
                        self.assembly.update({gene.id: new_gene})
            for gene in self.assembly.values():
                for model_id in gene.in_models["models_list"]:
                    getattr(self, model_id).update({gene.id: gene})


class SuperModel:  # TODO REAL 30.08.23 add transport reactions for periplasmic metabolites for models without periplasmic compartments
    """ Supermodel class with metabolites and reactions. Sources - names of original models used to create supermodel.
    Creating connections between metabolites and reaction via dictionaries with sources as keys and links to
    reactants/products/reactions as values.  """

    def __find_genes(
        self,
        all_models_data: dict,
        r_go_old_new: dict,
        r_go_new_old: dict,
        model_ids: [str],
        gene_folder: PosixPath,
    ):
        for model_id in model_ids:
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
            if gene_folder is not None:
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
            for reaction in self.reactions.assembly.values():
                old_rs = r_go_new_old.get(reaction.id).get(model_id)
                if old_rs:
                    new_gpr_unite_r = []
                    for oldr in old_rs:
                        if oldr.genes:
                            gene_convert = {}
                            for oldrg in oldr.genes:
                                if gene_folder is not None:
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
                                else:
                                    gene_convert.update({oldrg.id: oldrg.id})
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
            met._find_reactions(
                m_go_new_old, r_go_old_new, periplasmic_r, periplasmic_m
            )
        for r in self.reactions.assembly.values():
            r._find_reactants_products(
                r_go_new_old, m_go_old_new, periplasmic_r, "reactants"
            )
            r._find_reactants_products(
                r_go_new_old, m_go_old_new, periplasmic_r, "products"
            )
            r._find_metabolites(r_go_new_old, m_go_old_new, periplasmic_r)
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

    def __swapReactantsAndProducts(self, r: NewElement, sources_to_swap: list):
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
                    # len(consist) is expected to be only 1 or 2
                    else:
                        warnings.warn(
                            f"Warning! Something went wrong with swaping "
                            f"metabolites for {r.id}."
                            f"Can enter consist more 2."
                            f"Reactants: {r.reactants}. "
                            f"Products: {r.products}. Consist: {consist}"
                        )
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
        final_m_sel: dict,
        final_m_not_sel: dict,
        final_r_sel: dict,
        final_r_not_sel: dict,
        all_models_data: dict,
        additional_periplasmic_m: dict,
        periplasmic_r: dict,
        m_db_info: pd.core.frame.DataFrame,
        r_db_info: pd.core.frame.DataFrame,
        gene_folder,
        and_as_solid: bool,
    ):
        self.sources = list(all_models_data.keys())
        self.metabolites = SetofNewElements(
            "metabolites",
            final_m_sel,
            final_m_not_sel,
            self.sources,
            m_db_info,
            additional_periplasmic_m,
        )
        self.reactions = SetofNewElements(
            "reactions", final_r_sel, final_r_not_sel, self.sources, r_db_info,
        )
        self.genes = SetofNewGenes(all_models_data, gene_folder)

        m_go_old_new, m_go_new_old = self.metabolites._makeForwardBackward(
            all_models_data, final_m_sel, "metabolites", additional_periplasmic_m,
        )
        r_go_old_new, r_go_new_old = self.reactions._makeForwardBackward(
            all_models_data, final_r_sel, "reactions"
        )
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
                    f"Some of input models are not in supermodel. "
                    f"Maybe {wrong_yes} or {wrong_no}"
                    f"Please check the input ids"
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

    def get_all_confident_levels(self, and_as_solid=False):
        for i in range(len(self.sources), 1, -1):
            self.at_least_in(i, and_as_solid=and_as_solid)

    def write_supermodel_to_pkl(self, output_name: str, recursion_limit=None):
        if not output_name.endswith(".pkl"):
            raise ValueError("Wrong extension of the file")
        if exists(output_name):
            raise ValueError("File already exist, change the name")
        else:
            max_rec = 0x100000
            resource.setrlimit(
                resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY]
            )
            sys.setrecursionlimit(max_rec)
            with open(output_name, "wb") as fh:
                dill.dump(self, fh)


def read_supermodel_from_pkl(input_name: str):
    supermodel = dill.load(open(input_name, "rb"))
    return supermodel
