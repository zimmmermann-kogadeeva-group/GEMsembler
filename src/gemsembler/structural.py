import itertools
from collections import defaultdict
from copy import deepcopy
from itertools import combinations

import cobra

from .selection import Selected


class StructuralR(object):
    def __init__(self, bigg_structural: dict, selected: Selected):
        comment = list(bigg_structural.values())[0]
        if list(bigg_structural.keys())[0] == "NOT_found" or comment.startswith(
            "Potentially_found"
        ):
            self.structural = []
        else:
            self.structural = list(bigg_structural.keys())
        if len(bigg_structural) == 1 and comment.startswith(
            "Found_via_adding_periplasmic_compartment"
        ):
            self.comment = "Found_via_adding_periplasmic_compartment"
            self.suggestions = {"orig_m": [], "b_m": [], "p_m": []}
            if comment.split("-")[1].split(" ")[0] != "":
                self.suggestions["orig_m"] = self.suggestions["orig_m"] + comment.split(
                    "-"
                )[1].split(" ")
                self.suggestions["b_m"] = self.suggestions["b_m"] + comment.split("-")[
                    2
                ].split(" ")
                self.suggestions["p_m"] = self.suggestions["p_m"] + comment.split("-")[
                    3
                ].split(" ")
            if comment.split("-")[4].split(" ")[0] != "":
                self.suggestions["orig_m"] = self.suggestions["orig_m"] + comment.split(
                    "-"
                )[4].split(" ")
                self.suggestions["b_m"] = self.suggestions["b_m"] + comment.split("-")[
                    5
                ].split(" ")
                self.suggestions["p_m"] = self.suggestions["p_m"] + comment.split("-")[
                    6
                ].split(" ")
            self.compartments = list(bigg_structural.values())[0].split("-")[7].split()
        else:
            self.compartments = selected.compartments
            if len(bigg_structural) == 1 and comment.startswith(
                "Found_via_one_to_many_metabolites"
            ):
                self.comment = "Found_via_one_to_many_metabolites"
                self.suggestions = {
                    "orig_m": comment.split("-")[1].split(" "),
                    "b_m": comment.split("-")[2].split(" "),
                }
            elif len(bigg_structural) == 1 and comment.startswith(
                "Potentially_found_via_many_to_one_metabolites"
            ):
                self.comment = "Potentially_found_via_many_to_one_metabolites"
                self.suggestions = {
                    "reaction": list(bigg_structural.keys())[0],
                    "orig_m": [],
                    "b_m": [],
                }
                if comment.split("-")[1].split(" ")[0] != "":
                    self.suggestions["orig_m"] = self.suggestions[
                        "orig_m"
                    ] + comment.split("-")[1].split(" ")
                    self.suggestions["b_m"] = self.suggestions["b_m"] + comment.split(
                        "-"
                    )[2].split(" ")
                if comment.split("-")[3].split(" ")[0] != "":
                    self.suggestions["orig_m"] = self.suggestions[
                        "orig_m"
                    ] + comment.split("-")[3].split(" ")
                    self.suggestions["b_m"] = self.suggestions["b_m"] + comment.split(
                        "-"
                    )[4].split(" ")
            elif len(bigg_structural) == 1 and comment.startswith(
                (
                    "Not_found_via_many_to_one_metabolites",
                    "Potentially_found_but_no_confidence",
                    "Not_all_metabolites_are_converted_but_important_for_many_original",
                )
            ):
                if comment.startswith("Not_found_via_many_to_one_metabolites"):
                    list_name = "not_fit"
                else:
                    list_name = "no_data"
                self.comment = comment.split("-")[0]
                self.suggestions = {
                    "reaction": list(bigg_structural.keys())[0],
                    "orig_m": [],
                    "b_m": [],
                }
                if comment.split("-")[1].split(" ")[0] != "":
                    self.suggestions["orig_m"] = self.suggestions[
                        "orig_m"
                    ] + comment.split("-")[1].split(" ")
                    self.suggestions["b_m"] = self.suggestions["b_m"] + [
                        list_name
                    ] * len(comment.split("-")[1].split(" "))
                if comment.split("-")[2].split(" ")[0] != "":
                    self.suggestions["orig_m"] = self.suggestions[
                        "orig_m"
                    ] + comment.split("-")[2].split(" ")
                    self.suggestions["b_m"] = self.suggestions["b_m"] + [
                        list_name
                    ] * len(comment.split("-")[2].split(" "))
            else:
                self.suggestions = None
                if len(bigg_structural) == 1:
                    self.comment = comment
                else:
                    self.comment = f"Several result ids and comments: {' -NEW_COMMENT- '.join(list(bigg_structural.values()))}"


class StructuralM(object):
    def __init__(self, ids: list, comment: str, compartments: list):
        self.structural = ids
        self.compartments = compartments
        self.comment = comment


def getReaction(bigg_met1, bigg_met2, bigg_network_r, comment):
    """ Find reaction id from reaction's metabolites """

    bigg_met1_str = " ".join(sorted(bigg_met1))
    bigg_met2_str = " ".join(sorted(bigg_met2))
    bigg_equation = "<->".join(sorted([bigg_met1_str, bigg_met2_str]))
    if bigg_equation in list(bigg_network_r.keys()):
        return {bigg_network_r[bigg_equation]: comment}
    else:
        return {}


def Hydrogens(
    bigg_met1, bigg_met2, compart1, compart2, bigg_network_r, addit_comment=""
):
    """ Trying to convert reaction after add/removing H from/to reactants/products """
    bigg_r = {}
    if addit_comment != "":
        addit_comment = addit_comment + "-"
    for c1 in list(set(compart1)):
        bigg_met1_mod = deepcopy(bigg_met1)
        if "h_" + c1 in bigg_met1_mod:
            bigg_met1_mod.remove("h_" + c1)
            tmp_bigg_r = getReaction(
                bigg_met1_mod,
                bigg_met2,
                bigg_network_r,
                addit_comment + "Found_via_removing_H_from_reactants",
            )
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
            else:
                for cc2 in list(set(compart2)):
                    bigg_met2_h = deepcopy(bigg_met2)
                    bigg_met2_h.append("h_" + cc2)
                    tmp_bigg_r = getReaction(
                        bigg_met1_mod,
                        bigg_met2_h,
                        bigg_network_r,
                        addit_comment
                        + "Found_via_removing_H_from_reactants_adding_H_to_products",
                    )
                    if tmp_bigg_r:
                        bigg_r.update(tmp_bigg_r)
        else:
            bigg_met1_mod.append("h_" + c1)
            tmp_bigg_r = getReaction(
                bigg_met1_mod,
                bigg_met2,
                bigg_network_r,
                addit_comment + "Found_via_adding_H_to_reactants",
            )
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
    for c2 in list(set(compart2)):
        bigg_met2_mod = deepcopy(bigg_met2)
        if "h_" + c2 in bigg_met2_mod:
            bigg_met2_mod.remove("h_" + c2)
            tmp_bigg_r = getReaction(
                bigg_met1,
                bigg_met2_mod,
                bigg_network_r,
                addit_comment + "Found_via_removing_H_from_products",
            )
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
            else:
                for cc1 in list(set(compart1)):
                    bigg_met1_h = deepcopy(bigg_met1)
                    bigg_met1_h.append("h_" + cc1)
                    tmp_bigg_r = getReaction(
                        bigg_met1_h,
                        bigg_met2_mod,
                        bigg_network_r,
                        addit_comment
                        + "Found_via_adding_H_from_reactants_removing_H_to_products",
                    )
                    if tmp_bigg_r:
                        bigg_r.update(tmp_bigg_r)
        else:
            bigg_met2_mod.append("h_" + c2)
            tmp_bigg_r = getReaction(
                bigg_met1,
                bigg_met2_mod,
                bigg_network_r,
                addit_comment + "Found_via_adding_H_to_products",
            )
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
    return bigg_r


def Periplasmic(bigg_met1, bigg_met2, bigg_network_r):
    """ Trying to convert reaction after replacing compartments to periplasmic for all possible set of metabolites """
    bigg_r = {}
    bigg1_comb = []
    for r1 in range(len(bigg_met1) + 1):
        bigg1_comb = bigg1_comb + list(combinations(list(bigg_met1.keys()), r1))
    bigg2_comb = []
    for r2 in range(len(bigg_met2) + 1):
        bigg2_comb = bigg2_comb + list(combinations(list(bigg_met2.keys()), r2))
    for comb1 in bigg1_comb:
        for comb2 in bigg2_comb:
            bigg1 = list(set(bigg_met1) - set(comb1))
            c1 = list(set([b1[-1] for b1 in bigg1]))
            orig1 = [v for k, v in bigg_met1.items() if k in comb1]
            if comb1:
                bigg1_p = [m1[:-1] + "p" for m1 in comb1]
                c1.append("p")
            else:
                bigg1_p = []
            bigg2 = list(set(bigg_met2) - set(comb2))
            c2 = list(set([b2[-1] for b2 in bigg2]))
            orig2 = [vv for kk, vv in bigg_met2.items() if kk in comb2]
            if comb2:
                bigg2_p = [m2[:-1] + "p" for m2 in comb2]
                c2.append("p")
            else:
                bigg2_p = []
            tmp_bigg_r = getReaction(
                bigg1 + bigg1_p,
                bigg2 + bigg2_p,
                bigg_network_r,
                f"Found_via_adding_periplasmic_compartment-{' '.join(orig1)}-{' '.join(comb1)}-{' '.join(bigg1_p)}"
                f"-{' '.join(orig2)}-{' '.join(comb2)}-{' '.join(bigg2_p)}-{' '.join(list(set(c1 + c2)))}",
            )
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
    return bigg_r


def SeveralMetabolites(
    bigg_met1,
    bigg_met2,
    met1_to_many,
    met2_to_many,
    compart1,
    compart2,
    bigg_network_r,
):
    """ Trying to convert reaction with several variants coming from 1-n metabolites by variants from 1-n metabolites to 1-1 converted metabolites """
    bigg_r = {}
    bigg_r_h = {}
    bigg_to_many_variants1 = []
    orig_to_many_variants1 = []
    if met1_to_many != {}:
        for key1, val1 in met1_to_many.items():
            orig_to_many_variants1.append(key1)
            bigg_to_many_variants1.append([met_variant1 for met_variant1 in val1])
        bigg_to_many_variants1 = list(itertools.product(*bigg_to_many_variants1))
    bigg_to_many_variants2 = []
    orig_to_many_variants2 = []
    if met2_to_many != {}:
        for key2, val2 in met2_to_many.items():
            orig_to_many_variants2.append(key2)
            bigg_to_many_variants2.append([met_variant2 for met_variant2 in val2])
        bigg_to_many_variants2 = list(itertools.product(*bigg_to_many_variants2))
    if bigg_to_many_variants1 and bigg_to_many_variants2:
        for variant1 in bigg_to_many_variants1:
            for variant2 in bigg_to_many_variants2:
                bigg_met1_mod = list(bigg_met1.keys()) + list(variant1)
                bigg_met2_mod = list(bigg_met2.keys()) + list(variant2)
                tmp_bigg_r = getReaction(
                    bigg_met1_mod,
                    bigg_met2_mod,
                    bigg_network_r,
                    f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}",
                )
                tmp_bigg_r_h = Hydrogens(
                    bigg_met1_mod,
                    bigg_met2_mod,
                    compart1,
                    compart2,
                    bigg_network_r,
                    f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}",
                )

                if tmp_bigg_r:
                    bigg_r.update(tmp_bigg_r)
                if tmp_bigg_r_h:
                    bigg_r_h.update(tmp_bigg_r_h)
    elif bigg_to_many_variants1 and (not bigg_to_many_variants2):
        for variant1 in bigg_to_many_variants1:
            bigg_met1_mod = list(bigg_met1.keys()) + list(variant1)
            tmp_bigg_r = getReaction(
                bigg_met1_mod,
                list(bigg_met2.keys()),
                bigg_network_r,
                f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}",
            )
            tmp_bigg_r_h = Hydrogens(
                bigg_met1_mod,
                list(bigg_met2.keys()),
                compart1,
                compart2,
                bigg_network_r,
                f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}",
            )
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
            if tmp_bigg_r_h:
                bigg_r_h.update(tmp_bigg_r_h)
    elif bigg_to_many_variants2 and (not bigg_to_many_variants1):
        for variant2 in list(bigg_to_many_variants2):
            bigg_met2_mod = list(bigg_met2.keys()) + list(variant2)
            tmp_bigg_r = getReaction(
                list(bigg_met1.keys()),
                bigg_met2_mod,
                bigg_network_r,
                f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}",
            )
            tmp_bigg_r_h = Hydrogens(
                list(bigg_met1.keys()),
                bigg_met2_mod,
                compart1,
                compart2,
                bigg_network_r,
                f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}",
            )
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
            if tmp_bigg_r_h:
                bigg_r_h.update(tmp_bigg_r_h)
    if bigg_r:
        return bigg_r
    else:
        return bigg_r_h


def convertReactionViaNetworkStructure(
    orig_met1: list,
    orig_met2: list,
    selected_met: dict,
    bigg_network_r: dict,
    do_priplasmic: bool,
):
    """ Converting one reaction with several strategies. 1) Try just via it's metabolites if all them are converted 1-1.
     2) If not successful try to add/remove H 3) If not successful try to consider periplasmic compartment 4) If all
     metabolites in the reaction are converted 1-1 or 1-n try different variants from this n options. Writing comments
     with information how it was converted and what happened to metabolites for further suggestions"""
    compart1 = []
    bigg_met1 = {}
    bigg_met1_to_many = {}
    bigg_met1_from_many = {}
    for met1 in orig_met1:
        compart1 = compart1 + selected_met[met1].compartments
        if (
            selected_met[met1].to_one_id == True
            and selected_met[met1].from_one_id == True
        ):
            bigg_met1.update({selected_met[met1].highest_consistent[0]: met1})
        elif (
            selected_met[met1].to_one_id == False
            and selected_met[met1].from_one_id == True
        ):
            bigg_met1_to_many.update({met1: selected_met[met1].highest_consistent})
        elif (
            selected_met[met1].to_one_id == True
            and selected_met[met1].from_one_id == False
        ):
            bigg_met1_from_many.update({met1: selected_met[met1].highest_consistent[0]})
    compart2 = []
    bigg_met2 = {}
    bigg_met2_to_many = {}
    bigg_met2_from_many = {}
    for met2 in orig_met2:
        compart2 = compart2 + selected_met[met2].compartments
        if (
            selected_met[met2].to_one_id == True
            and selected_met[met2].from_one_id == True
        ):
            bigg_met2.update({selected_met[met2].highest_consistent[0]: met2})
        elif (
            selected_met[met2].to_one_id == False
            and selected_met[met2].from_one_id == True
        ):
            bigg_met2_to_many.update({met2: selected_met[met2].highest_consistent})
        elif (
            selected_met[met2].to_one_id == True
            and selected_met[met2].from_one_id == False
        ):
            bigg_met2_from_many.update({met2: selected_met[met2].highest_consistent[0]})
    if (len(bigg_met1) == len(orig_met1)) & (len(bigg_met2) == len(orig_met2)):
        bigg_r = getReaction(
            list(bigg_met1.keys()),
            list(bigg_met2.keys()),
            bigg_network_r,
            "Found_via_pure_reaction_equation",
        )
        if not bigg_r:
            bigg_r = Hydrogens(
                list(bigg_met1.keys()),
                list(bigg_met2.keys()),
                compart1,
                compart2,
                bigg_network_r,
            )
            if not bigg_r:
                if do_priplasmic:
                    bigg_r = Periplasmic(bigg_met1, bigg_met2, bigg_network_r,)
                if not bigg_r:
                    bigg_r = {"NOT_found": "Not_found_in_BiGG_network"}
    elif (len(bigg_met1) + len(bigg_met1_to_many) == len(orig_met1)) & (
        len(bigg_met2) + len(bigg_met2_to_many) == len(orig_met2)
    ):
        bigg_r = SeveralMetabolites(
            bigg_met1,
            bigg_met2,
            bigg_met1_to_many,
            bigg_met2_to_many,
            compart1,
            compart2,
            bigg_network_r,
        )
        if not bigg_r:
            bigg_r = {"NOT_found": "Not_found_via_one_to_many_metabolites"}
    elif (len(bigg_met1) + len(bigg_met1_from_many) == len(orig_met1)) & (
        len(bigg_met2) + len(bigg_met2_from_many) == len(orig_met2)
    ):
        if len(bigg_met1) + len(bigg_met2) == 0:
            bigg_r = getReaction(
                list(bigg_met1.keys()) + list(bigg_met1_from_many.values()),
                list(bigg_met2.keys()) + list(bigg_met2_from_many.values()),
                bigg_network_r,
                f"Potentially_found_but_no_confidence-"
                f"{' '.join(list(bigg_met1_from_many.keys()))}-{' '.join(list(bigg_met2_from_many.keys()))}",
            )
        else:
            bigg_r = getReaction(
                list(bigg_met1.keys()) + list(bigg_met1_from_many.values()),
                list(bigg_met2.keys()) + list(bigg_met2_from_many.values()),
                bigg_network_r,
                f"Potentially_found_via_many_to_one_metabolites-"
                f"{' '.join(list(bigg_met1_from_many.keys()))}-{' '.join(list(bigg_met1_from_many.values()))}-"
                f"{' '.join(list(bigg_met2_from_many.keys()))}-{' '.join(list(bigg_met2_from_many.values()))}",
            )
        if not bigg_r:
            bigg_r = {
                "NOT_found": f"Not_found_via_many_to_one_metabolites-"
                f"{' '.join(list(bigg_met1_from_many.keys()))}-{' '.join(list(bigg_met2_from_many.keys()))}"
            }
    elif len(bigg_met1_from_many) + len(bigg_met2_from_many) > 0:
        bigg_r = {
            "NOT_found": f"Not_all_metabolites_are_converted_but_important_for_many_original-"
            f"{' '.join(list(bigg_met1_from_many.keys()))}-{' '.join(list(bigg_met2_from_many.keys()))}"
        }
    else:
        bigg_r = {"NOT_found": "Not_all_metabolites_for_reaction_are_converted"}
    return bigg_r


def runStructuralCheck(
    react_checked: dict,
    met_checked: dict,
    model: cobra.core.model.Model,
    bigg_network: dict,
):
    """ Checking reactions equations for models with no conversion need (with BiGG ids originally). Should be in BiGG
    database if not exchange or biomass reaction. """
    react_struct_checked = {}
    for id_r, sel in react_checked.items():
        bigg_met1 = [met1.id for met1 in model.reactions.get_by_id(id_r).reactants]
        bigg_met2 = [met2.id for met2 in model.reactions.get_by_id(id_r).products]
        bigg_met1_checked = []
        bigg_met2_checked = []
        for m1 in bigg_met1:
            if (
                met_checked[m1].to_one_id == True
                and met_checked[m1].from_one_id == True
            ):
                bigg_met1_checked.append(met_checked[m1].highest_consistent[0])
        for m2 in bigg_met2:
            if (
                met_checked[m2].to_one_id == True
                and met_checked[m2].from_one_id == True
            ):
                bigg_met2_checked.append(met_checked[m2].highest_consistent[0])
        if (len(bigg_met1) > 20) | (len(bigg_met2) > 20):
            react_struct_checked.update(
                {id_r: StructuralR({"Biomass": "growth_reaction"}, sel)}
            )
            continue
        bigg_r = getReaction(
            bigg_met1_checked,
            bigg_met2_checked,
            bigg_network,
            "checked_via_structural_reaction_equation",
        )
        if bigg_r:
            react_struct_checked.update({id_r: StructuralR(bigg_r, sel)})
        elif sel.highest_consistent:
            react_struct_checked.update(
                {
                    id_r: StructuralR(
                        {sel.highest_consistent[0]: "id_in_bigg_originally"}, sel
                    )
                }
            )
        elif (
            (len(bigg_met1) == 1) & (len(bigg_met2) == 0)
            | (len(bigg_met1) == 0) & (len(bigg_met2) == 1)
        ) & (id_r.startswith("EX_")):
            react_struct_checked.update(
                {id_r: StructuralR({id_r: "not_found_but_exchange_reaction"}, sel)}
            )
        else:
            react_struct_checked.update(
                {
                    id_r: StructuralR(
                        {"NOT_found": f"{id_r} not_pass_id_and_structural_checking"},
                        sel,
                    )
                }
            )
    return react_struct_checked


def runStructuralConversion(
    model_db: str,
    first_stage_selected_r: dict,
    first_stage_selected_m: dict,
    model: cobra.core.model.Model,
    bigg_network: dict,
    models_periplasmic: bool,
):
    """ Running structural conversion for all reactions. Selection reactions that have only 1 id as result """
    if model_db == "bigg":
        structural_conversion_r = runStructuralCheck(
            first_stage_selected_r, first_stage_selected_m, model, bigg_network
        )
    else:
        structural_conversion_r = {}
        for orig_id, selected in first_stage_selected_r.items():
            orig_met1 = [
                react.id for react in model.reactions.get_by_id(orig_id).reactants
            ]
            orig_met2 = [pro.id for pro in model.reactions.get_by_id(orig_id).products]
            if (len(orig_met1) > 20) | (len(orig_met2) > 20):
                structural_bigg_id = {
                    "Biomass": "No structural conversion since growth_reaction"
                }
            else:
                structural_bigg_id = convertReactionViaNetworkStructure(
                    orig_met1,
                    orig_met2,
                    first_stage_selected_m,
                    bigg_network,
                    models_periplasmic,
                )
            structural_conversion_r.update(
                {orig_id: StructuralR(structural_bigg_id, selected)}
            )
    return structural_conversion_r


def runSuggestionsMet(
    model_db: str, structural_rs: dict, struct_rs_sel: dict, met_selected: dict,
):
    # model uses bigg originally so there was no structural conversion, only checking
    if model_db == "bigg":
        structural_met = {
            orig_id: StructuralM(
                sel.highest_consistent,
                "No need for structural, because of bigg database",
                sel.compartments,
            )
            for orig_id, sel in met_selected.items()
        }
        sug_from_many = {}
    # actually getting suggestions based on reactions that were converted 1-1 with reaction equation
    else:
        structural_met = {}
        sug_to_many = defaultdict(list)
        sug_from_many = defaultdict(list)
        for r_old_id, r_sel in struct_rs_sel.items():
            if (
                r_sel.to_one_id == True
                and r_sel.from_one_id == True
                and structural_rs[r_old_id].comment.startswith(
                    "Found_via_one_to_many_metabolites"
                )
            ):
                for i in range(len(structural_rs[r_old_id].suggestions["orig_m"])):
                    sug_to_many[
                        structural_rs[r_old_id].suggestions["orig_m"][i]
                    ].append(structural_rs[r_old_id].suggestions["b_m"][i])
            elif structural_rs[r_old_id].comment.startswith(
                (
                    "Potentially_found_via_many_to_one_metabolites",
                    "Not_found_via_many_to_one_metabolites",
                    "Potentially_found_but_no_confidence",
                    "Not_all_metabolites_are_converted_but_important_for_many_original",
                )
            ):
                for ii in range(len(structural_rs[r_old_id].suggestions["orig_m"])):
                    sug_from_many[
                        structural_rs[r_old_id].suggestions["orig_m"][ii]
                    ].append(structural_rs[r_old_id].suggestions["b_m"][ii])
        sug_to_many = {k: v[0] for k, v in sug_to_many.items() if len(set(v)) == 1}
        for m_orig_id, met_sel in met_selected.items():
            if met_sel.to_one_id == True and met_sel.from_one_id == True:
                structural_met.update(
                    {
                        m_orig_id: StructuralM(
                            met_sel.highest_consistent,
                            "Already 1-1 conversion",
                            met_sel.compartments,
                        )
                    }
                )
            elif met_sel.to_one_id == False and met_sel.from_one_id == True:
                if m_orig_id in sug_to_many.keys():
                    structural_met.update(
                        {
                            m_orig_id: StructuralM(
                                [sug_to_many[m_orig_id]],
                                "Suggestion from one to many options",
                                met_sel.compartments,
                            )
                        }
                    )
                else:
                    structural_met.update(
                        {
                            m_orig_id: StructuralM(
                                [],
                                "No suggestion from one to many options",
                                met_sel.compartments,
                            )
                        }
                    )
            elif met_sel.to_one_id == True and met_sel.from_one_id == False:
                # checking how the same original id was converted in other models with different db
                # no information about the id in models with other db
                if not met_sel.in_other_models:
                    in_others = False
                else:
                    # checking whether in models with other db the id is converted as 1-1 ([True, True])
                    in_others = True
                    for val in met_sel.in_other_models.values():
                        if val != [True, True]:
                            in_others = False
                    # if the id is 1-1 in other models, checking whether other ids from the model, group that lead to n-1, are not 1-1 in models with other db
                    if in_others:
                        for oth_id in met_sel.from_many_other_ids:
                            for othval in met_selected[oth_id].in_other_models.values():
                                if othval == [True, True]:
                                    in_others = False
                if in_others:
                    structural_met.update(
                        {
                            m_orig_id: StructuralM(
                                met_sel.highest_consistent,
                                "Suggestion from many original ids to one bigg based on unique conversion from other models",
                                met_sel.compartments,
                            )
                        }
                    )
                # checking whether usage of conversion for one original ids from the group leads to better reaction equation conversion
                else:
                    in_r_eq = False
                    if m_orig_id in sug_from_many.keys():
                        if ("not_fit" not in set(sug_from_many[m_orig_id])) and (
                            {"no_data"} != set(sug_from_many[m_orig_id])
                        ):
                            in_r_eq = True
                            for ot_id in met_sel.from_many_other_ids:
                                if not {"not_fit", "no_data"} >= set(
                                    sug_from_many[ot_id]
                                ):
                                    in_r_eq = False
                    if in_r_eq:
                        structural_met.update(
                            {
                                m_orig_id: StructuralM(
                                    met_sel.highest_consistent,
                                    "Suggestion from many original ids to one bigg based on better reaction equation conversions for that original id in a group",
                                    met_sel.compartments,
                                )
                            }
                        )
                    else:
                        structural_met.update(
                            {
                                m_orig_id: StructuralM(
                                    [],
                                    "No suggestion from many original ids to one bigg",
                                    met_sel.compartments,
                                )
                            }
                        )
            elif met_sel.to_one_id == False and met_sel.from_one_id == False:
                structural_met.update(
                    {
                        m_orig_id: StructuralM(
                            [],
                            "No suggestion from many original ids to many bigg ids",
                            met_sel.compartments,
                        )
                    }
                )
            else:
                structural_met.update(
                    {
                        m_orig_id: StructuralM(
                            [],
                            "No suggestion from not converted",
                            met_sel.compartments,
                        )
                    }
                )
    return structural_met, sug_from_many
