import itertools
from itertools import groupby, combinations
import operator
from copy import deepcopy
import pandas as pd
import cobra
from collections import Counter
import general


def findReactionViaPureEquation(bigg_met1, bigg_met2, BiGG_network_r, comment):
    bigg_met1_str = " ".join(sorted(bigg_met1))
    bigg_met2_str = " ".join(sorted(bigg_met2))
    bigg_r = BiGG_network_r[
        ((BiGG_network_r["1metabolites"] == bigg_met1_str) & (BiGG_network_r["2metabolites"] == bigg_met2_str)) | (
                (BiGG_network_r["1metabolites"] == bigg_met2_str) & (BiGG_network_r["2metabolites"] == bigg_met1_str))][
        "reaction"]
    if not bigg_r.empty:
        return {bigg_r.values[0]: comment}
    else:
        return {}


def workWithH(bigg_met1, bigg_met2, compart1, compart2, BiGG_network_r, addit_comment=""):
    bigg_r = {}
    if addit_comment != "": addit_comment = addit_comment + "-"
    for c1 in list(set(compart1)):
        bigg_met1_mod = deepcopy(bigg_met1)
        if "h_" + c1 in bigg_met1_mod:
            bigg_met1_mod.remove("h_" + c1)
            tmp_bigg_r = findReactionViaPureEquation(bigg_met1_mod, bigg_met2, BiGG_network_r,
                                                     addit_comment + "Found_via_removing_H_from_reactants")
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
            else:
                for cc2 in list(set(compart2)):
                    bigg_met2_h = deepcopy(bigg_met2)
                    bigg_met2_h.append("h_" + cc2)
                    tmp_bigg_r = findReactionViaPureEquation(bigg_met1_mod, bigg_met2_h, BiGG_network_r,
                                                            addit_comment + "Found_via_removing_H_from_reactants_adding_H_to_products")
                    if tmp_bigg_r: bigg_r.update(tmp_bigg_r)
        else:
            bigg_met1_mod.append("h_" + c1)
            tmp_bigg_r = findReactionViaPureEquation(bigg_met1_mod, bigg_met2, BiGG_network_r,
                                                     addit_comment + "Found_via_adding_H_to_reactants")
            if tmp_bigg_r: bigg_r.update(tmp_bigg_r)
    for c2 in list(set(compart2)):
        bigg_met2_mod = deepcopy(bigg_met2)
        if "h_" + c2 in bigg_met2_mod:
            bigg_met2_mod.remove("h_" + c2)
            tmp_bigg_r = findReactionViaPureEquation(bigg_met1, bigg_met2_mod, BiGG_network_r,
                                                     addit_comment + "Found_via_removing_H_from_products")
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
            else:
                for cc1 in list(set(compart1)):
                    bigg_met1_h = deepcopy(bigg_met1)
                    bigg_met1_h.append("h_" + cc1)
                    tmp_bigg_r = findReactionViaPureEquation(bigg_met1_h, bigg_met2_mod, BiGG_network_r,
                                                             addit_comment + "Found_via_adding_H_from_reactants_removing_H_to_products")
                    if tmp_bigg_r: bigg_r.update(tmp_bigg_r)
        else:
            bigg_met2_mod.append("h_" + c2)
            tmp_bigg_r = findReactionViaPureEquation(bigg_met1, bigg_met2_mod, BiGG_network_r,
                                                     addit_comment + "Found_via_adding_H_to_products")
            if tmp_bigg_r: bigg_r.update(tmp_bigg_r)
    return bigg_r


def addPeriplasmic(bigg_met1, bigg_met2, compart1, compart2, BiGG_network_r, addit_comment=""):

    bigg_r = {}
    # bigg_r_H = {}
    bigg1_comb = []
    for r1 in range(len(bigg_met1)+1):
        bigg1_comb = bigg1_comb + list(combinations(bigg_met1, r1))
    bigg2_comb = []
    for r2 in range(len(bigg_met2)+1):
        bigg2_comb = bigg2_comb + list(combinations(bigg_met2, r2))
    for comb1 in bigg1_comb:
        for comb2 in bigg2_comb:
            bigg1 = list(set(bigg_met1) - set(comb1))
            if comb1:
                bigg1_p = [m1[:-1]+"p" for m1 in comb1]
            else:
                bigg1_p = []
            bigg2 = list(set(bigg_met2) - set(comb2))
            if comb2:
                bigg2_p = [m2[:-1]+"p" for m2 in comb2]
            else:
                bigg2_p = []
            tmp_bigg_r = findReactionViaPureEquation(bigg1+bigg1_p, bigg2+bigg2_p, BiGG_network_r,
                                                     "Found_via_adding_periplasmic_compartment")
            # tmp_bigg_r_H = workWithH(bigg1+bigg1_p, bigg2+bigg2_p, compart1+["p"], compart2+["p"], BiGG_network_r,
            #                          "Found_via_adding_periplasmic_compartment")
            if tmp_bigg_r: bigg_r.update(tmp_bigg_r)
            # if tmp_bigg_r_H: bigg_r_H.update(tmp_bigg_r_H)
    if bigg_r:
        return bigg_r
    # else:
    #     return bigg_r_H


def trySeveralMetabolites(bigg_met1, bigg_met2, met1_to_many, met2_to_many, compart1, compart2, BiGG_network_r):
    bigg_r = {}
    bigg_r_H = {}
    bigg_to_many_variants1 = []
    orig_to_many_variants1 = []
    if met1_to_many != {}:
        for key1, val1 in met1_to_many.items():
            orig_to_many_variants1.append(key1)
            bigg_to_many_variants1.append([met_variant1 for met_variant1 in val1[1]])
        bigg_to_many_variants1 = list(itertools.product(*bigg_to_many_variants1))
    bigg_to_many_variants2 = []
    orig_to_many_variants2 = []
    if met2_to_many != {}:
        for key2, val2 in met2_to_many.items():
            orig_to_many_variants2.append(key2)
            bigg_to_many_variants2.append([met_variant2 for met_variant2 in val2[1]])
        bigg_to_many_variants2 = list(itertools.product(*bigg_to_many_variants2))
    if bigg_to_many_variants1 and bigg_to_many_variants2:
        for variant1 in bigg_to_many_variants1:
            for variant2 in bigg_to_many_variants2:
                bigg_met1_mod = bigg_met1 + list(variant1)
                bigg_met2_mod = bigg_met2 + list(variant2)
                tmp_bigg_r = findReactionViaPureEquation(bigg_met1_mod, bigg_met2_mod, BiGG_network_r,
                                                         f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}")
                tmp_bigg_r_H = workWithH(bigg_met1_mod, bigg_met2_mod, compart1, compart2, BiGG_network_r,
                                         f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}")

                if tmp_bigg_r: bigg_r.update(tmp_bigg_r)
                if tmp_bigg_r_H: bigg_r_H.update(tmp_bigg_r_H)
    if bigg_to_many_variants1 and (not bigg_to_many_variants2):
        for variant1 in bigg_to_many_variants1:
            bigg_met1_mod = bigg_met1 + list(variant1)
            tmp_bigg_r = findReactionViaPureEquation(bigg_met1_mod, bigg_met2, BiGG_network_r,
                                                     f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}")
            tmp_bigg_r_H = workWithH(bigg_met1_mod, bigg_met2, compart1, compart2, BiGG_network_r,
                                   f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}")
            if tmp_bigg_r: bigg_r.update(tmp_bigg_r)
            if tmp_bigg_r_H: bigg_r_H.update(tmp_bigg_r_H)
    if bigg_to_many_variants2 and (not bigg_to_many_variants1):
        for variant2 in list(bigg_to_many_variants2):
            bigg_met2_mod = bigg_met2 + list(variant2)
            tmp_bigg_r = findReactionViaPureEquation(bigg_met1, bigg_met2_mod, BiGG_network_r,
                                                     f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}")
            tmp_bigg_r_H = workWithH(bigg_met1, bigg_met2_mod, compart1, compart2, BiGG_network_r,
                                   f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}")
            if tmp_bigg_r: bigg_r.update(tmp_bigg_r)
            if tmp_bigg_r_H: bigg_r_H.update(tmp_bigg_r_H)
    if bigg_r:
        return bigg_r
    else:
        return bigg_r_H


def convertReactionViaNetworkStructure(reaction_id: str, model: cobra.core.model.Model,
                                       first_confidence_conversion: dict, highest: dict,
                                       BiGG_network_r: pd.core.frame.DataFrame, one_to_many_conversion=None):
    orig_met1 = [react.id for react in model.reactions.get_by_id(reaction_id).reactants]
    bigg_met1 = []
    compart1 = []
    met1_to_many = {}
    for met1 in orig_met1:
        compart1.append(highest[met1][0][0])
        if met1 in first_confidence_conversion:
            bigg_met1.append(first_confidence_conversion[met1][1][0])
        if one_to_many_conversion:
            if met1 in one_to_many_conversion:
                met1_to_many.update({met1: one_to_many_conversion[met1]})
    orig_met2 = [pro.id for pro in model.reactions.get_by_id(reaction_id).products]
    bigg_met2 = []
    compart2 = []
    met2_to_many = {}
    for met2 in orig_met2:
        compart2.append(highest[met2][0][0])
        if met2 in first_confidence_conversion:
            bigg_met2.append(first_confidence_conversion[met2][1][0])
        if one_to_many_conversion:
            if met2 in one_to_many_conversion:
                met2_to_many.update({met2: one_to_many_conversion[met2]})
    if (len(bigg_met1) == len(orig_met1)) & (len(bigg_met2) == len(orig_met2)):
        bigg_r = findReactionViaPureEquation(bigg_met1, bigg_met2, BiGG_network_r,
                                                 "Found_via_pure_reaction_equation")
        if not bigg_r:
            bigg_r = workWithH(bigg_met1, bigg_met2, compart1, compart2, BiGG_network_r)
            if not bigg_r:
                bigg_r = addPeriplasmic(bigg_met1, bigg_met2, compart1, compart2, BiGG_network_r)
                if not bigg_r:
                    bigg_r = {"NOT_found": "Not_found_in_BiGG_network"}
    elif (len(bigg_met1) + len(met1_to_many.keys()) == len(orig_met1)) & (
            len(bigg_met2) + len(met2_to_many.keys()) == len(orig_met2)):
        bigg_r = trySeveralMetabolites(bigg_met1, bigg_met2, met1_to_many, met2_to_many, compart1, compart2, BiGG_network_r)
        if not bigg_r:
            bigg_r = {"NOT_found": "Not_found_via_one_to_many_metabolites"}
    else:
        bigg_r = {"NOT_found": "Not_all_metabolites_for_reaction_are_converted"}
    return bigg_r


def runAdditionalConversion(model_types: [str], first_stage_selected_m: dict, first_stage_selected_r: dict,
                            all_models: dict, bigg_network: dict,
                            obj_type: "metabolites" or "reactions"):
    structural_conversion = {}
    structural_conversion_to_one = {}
    for typ in model_types:
        structural_conversion.update({typ: {}})
        structural_conversion_to_one.update({typ: {}})
        if obj_type == "reactions":
            one_one_m = first_stage_selected_m.get("one_to_one").get(typ)
            one_to_many_m = first_stage_selected_m.get("one_to_many").get(typ)
            highest_m = first_stage_selected_m.get("intermediate_data").get("highest").get(typ)
            for key, additional in first_stage_selected_r.items():
                if typ in additional.keys():
                    for orig_id, select_bigg_id in additional.get(typ).items():
                        structural_bigg_id = convertReactionViaNetworkStructure(orig_id, all_models.get(typ),
                                                                                one_one_m, highest_m,
                                                                                bigg_network.get(obj_type), one_to_many_m)
                        structural_conversion.get(typ).update({orig_id: [select_bigg_id, structural_bigg_id, {"first_selection_type": key}]})
                        if (len(structural_bigg_id.keys()) == 1) & (list(structural_bigg_id.keys())[0] != "NOT_found"):
                            structural_conversion_to_one.get(typ).update({orig_id: [select_bigg_id[0], list(structural_bigg_id.keys())]})

    return structural_conversion, structural_conversion_to_one


def getSuggestionForMetabolites(model_types: [str], structural_r_one_one: dict, structural_r_all: dict):
    suggestions_m = {}
    suggestions_m_sel = {}
    for typ in model_types:
        suggestions_m.update({typ: {}})
        suggestions_m_sel.update({typ: {}})
        for orig_id, struct_id in structural_r_one_one.get(typ).items():
            struct_id = struct_id[1]
            comment = structural_r_all.get(typ).get(orig_id)[1].get(struct_id)
            if comment.startswith("Found_via_one_to_many_metabolites"):
                orig_m = comment.split("-")[1].split(" ")
                suggest_m = comment.split("-")[2].split(" ")
                for i in range(len(orig_m)):
                    if orig_m[i] in suggestions_m.get(typ).keys():
                        suggestions_m.get(typ).get(orig_m[i]).append(suggest_m[i])
                    else:
                        suggestions_m.get(typ).update({orig_m[i]: [suggest_m[i]]})
        for key, value in suggestions_m.get(typ).items():
            occurrence = Counter(value)
            suggestions_m.get(typ)[key] = {variant: number for variant, number in occurrence.items()}
            if len(occurrence.keys()) == 1:
                suggestions_m_sel.get(typ).update({key: list(occurrence.keys())})
    return suggestions_m, suggestions_m_sel


def ManyOneMetFromStructuralMet(many_to_one: dict, model: cobra.core.model.Model,
                                       first_structural_met: dict,
                                       BiGG_network_r: pd.core.frame.DataFrame):
    by_value = operator.itemgetter(1)
    grouped_many_one = [dict(g) for k, g in groupby(sorted(many_to_one.items(), key=by_value), by_value)]
    many_to_one_suggestions = {}
    for group in grouped_many_one:
        found_group = {}
        for orig_id, value in group.items():
            reactions = model.metabolites.get_by_id(orig_id).reactions
            found = 0
            n = 0
            for reaction in reactions:
                orig_met1 = [m1.id for m1 in reaction.reactants]
                orig_met2 = [m2.id for m2 in reaction.products]
                orig_met = [m.id for m in reaction.metabolites if m.id != orig_id]
                if (orig_met != []) & (set(orig_met) <= set(first_structural_met.keys())):
                    n = n + 1
                    bigg_met1 = []
                    for met1 in orig_met1:
                        if met1 == orig_id:
                            bigg_met1.append(value[1])
                        else:
                            bigg_met1.append(first_structural_met.get(met1)[1])
                    bigg_met2 = []
                    for met2 in orig_met2:
                        if met2 == orig_id:
                            bigg_met2.append(value[1])
                        else:
                            bigg_met2.append(first_structural_met.get(met2)[1])
                    tmp_bigg_r = findReactionViaPureEquation(bigg_met1, bigg_met2, BiGG_network_r, "found")
                    if tmp_bigg_r: found = found + 1
            if n != 0: found_group.update({orig_id: found/n})
            if n == 0: found_group.update({orig_id: "met_not_struct"})
        occur = Counter(list(found_group.values()))
        if (occur[1.0] == 1) & (set(occur.keys()) <= set([1.0, 0.0, "met_not_struct"])):
            recommend = general.findKeysByValue(found_group, 1, operator.eq)[0]
            many_to_one_suggestions.update({recommend: group.get(recommend)})
    return many_to_one_suggestions

def test(model_types: [str], final_metabolites: dict, first_stage_selected_m: dict, first_stage_selected_r: dict,
                            all_models: dict, bigg_network: dict,
                            obj_type: "metabolites" or "reactions"):
    structural_conversion = {}
    structural_conversion_to_one = {}
    for typ in model_types:
        structural_conversion.update({typ: {}})
        structural_conversion_to_one.update({typ: {}})
        if obj_type == "reactions":
            highest_m = first_stage_selected_m.get("intermediate_data").get("highest").get(typ)
            for key, additional in first_stage_selected_r.items():
                if typ in additional.keys():
                    for orig_id, select_bigg_id in additional.get(typ).items():
                        structural_bigg_id = convertReactionViaNetworkStructure(orig_id, all_models.get(typ),
                                                                                final_metabolites.get(typ), highest_m,
                                                                                bigg_network.get(obj_type))
                        structural_conversion.get(typ).update({orig_id: [select_bigg_id, structural_bigg_id, {"first_selection_type": key}]})
                        if (len(structural_bigg_id.keys()) == 1) & (list(structural_bigg_id.keys())[0] != "NOT_found"):
                            structural_conversion_to_one.get(typ).update({orig_id: [select_bigg_id[0], list(structural_bigg_id.keys())]})
    return structural_conversion, structural_conversion_to_one
