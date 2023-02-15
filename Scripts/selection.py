import copy
import itertools
from itertools import groupby
import operator
from copy import deepcopy
import pandas as pd
import cobra
from collections import Counter
import general


def getHighestConversion(model_types: [str], all_converted: dict):
    highest_conversion_level = {}
    for typ in model_types:
        highest_conversion_level.update({typ: {}})
        for id, conv in all_converted.get(typ).items():
            for level, converted in conv[1].items():
                if len(converted) > 0:
                    highest_conversion_level.get(typ).update({id: [conv[0], converted]})
                    break
                if level == "6-NOconv":
                    highest_conversion_level.get(typ).update({id: [conv[0], []]})
    return highest_conversion_level


def checkSameConversion(models_same_db: dict, highest_converted: dict, obj_type: "metabolites" or "reactions",
                        write_files=True, do_stat=True, useroutname=None):
    if useroutname is not None:
        filename_changed = "../Output/" + useroutname + "_" + obj_type + "_changed_for_consistency.tsv"
        filename_notconsistent = "../Output/" + useroutname + "_" + obj_type + "_not_consistent.tsv"
    else:
        filename_changed = "../Output/" + obj_type + "_changed_for_consistency.tsv"
        filename_notconsistent = "../Output/" + obj_type + "_not_consistent.tsv"
    consistent_highest = {}
    not_consistent_highest = {}
    none_consistent_highest = pd.DataFrame(
        columns=["Database", "Model", "original_id", "incommon_bigg_id", "individual_bigg_ids"])
    changed = pd.DataFrame(columns=["Database", "Model", "original_id", "incommon_bigg_id", "individual_bigg_ids"])
    changed_stat = {}
    for bd, models in models_same_db.items():
        bd_ids = []
        for model in models:
            consistent_highest.update({model: {}})
            not_consistent_highest.update({model: {}})
            changed_stat.update(
                {model: {"total_change": 0, "n_to_one": 0, "zero_to_one": 0, "n_to_n": 0, "zero_to_n": 0}})
            bd_ids = bd_ids + list(highest_converted.get(model).keys())
        bd_ids = list(set(bd_ids))
        for id in bd_ids:
            bigg_ids = []
            types_present = []
            for typ in models:
                if highest_converted.get(typ).get(id) is not None:
                    types_present.append(typ)
                    if highest_converted.get(typ).get(id)[1]:
                        bigg_ids.append(highest_converted.get(typ).get(id)[1])
            if not bigg_ids:
                for pres in types_present:
                    consistent_highest.get(pres).update({id: [highest_converted.get(pres).get(id)[0], []]})
            else:
                common_ids = list(set.intersection(*map(set, bigg_ids)))
                if len(common_ids) > 0:
                    for present in types_present:
                        if common_ids != highest_converted.get(present).get(id)[1]:
                            consistent_highest.get(present).update(
                                {id: [highest_converted.get(present).get(id)[0], common_ids]})
                            changed = pd.concat([changed,
                                                 pd.DataFrame(
                                                     [[bd, present, id, common_ids,
                                                       highest_converted.get(present).get(id)[1]]],
                                                     columns=changed.columns)], ignore_index=True)
                            changed_stat.get(present)["total_change"] = changed_stat.get(present)["total_change"] + 1
                            if (len(common_ids) == 1) & (len(highest_converted.get(present).get(id)[1]) > 1):
                                changed_stat.get(present)["n_to_one"] = changed_stat.get(present)["n_to_one"] + 1
                            if (len(common_ids) == 1) & (len(highest_converted.get(present).get(id)[1]) == 0):
                                changed_stat.get(present)["zero_to_one"] = changed_stat.get(present)["zero_to_one"] + 1
                            if (len(common_ids) > 1) & (len(highest_converted.get(present).get(id)[1]) > 1):
                                changed_stat.get(present)["n_to_n"] = changed_stat.get(present)["n_to_n"] + 1
                            if (len(common_ids) > 1) & (len(highest_converted.get(present).get(id)[1]) == 0):
                                changed_stat.get(present)["zero_to_n"] = changed_stat.get(present)["zero_to_n"] + 1
                        else:
                            consistent_highest.get(present).update({id: [highest_converted.get(present).get(id)[0],
                                                                         highest_converted.get(present).get(id)[1]]})
                else:
                    for pr in types_present:
                        not_consistent_highest.get(pr).update(
                            {id: [highest_converted.get(pr).get(id)[0], highest_converted.get(pr).get(id)[1]]})
                        none_consistent_highest = pd.concat(
                            [none_consistent_highest,
                             pd.DataFrame([[bd, pr, id, common_ids, highest_converted.get(pr).get(id)[1]]],
                                          columns=none_consistent_highest.columns)], ignore_index=True)
    if write_files:
        changed.to_csv(filename_changed, sep='\t')
        none_consistent_highest.to_csv(filename_notconsistent, sep='\t')
    if do_stat:
        for key, value in changed_stat.items():  # TODO: write this stat in file
            print(f"For {key} model in total {value['total_change']} {obj_type} were changed")
            print(f"For {key} model {value['n_to_one']} {obj_type} were changed from several ids to one")
            print(f"For {key} model {value['zero_to_one']} {obj_type} were changed from no found ids to one")
            print(f"For {key} model {value['n_to_n']} {obj_type} were changed from several ids to smaller amount of ids")
            print(f"For {key} model {value['zero_to_n']} {obj_type} were changed from no found ids to several")
    return consistent_highest, not_consistent_highest


def checkOneToManyHighest(model_types: [str], consistent_highest: dict):
    to_one = {}
    to_many = {}
    not_converted = {}
    for typ in model_types:
        to_one.update({typ: {}})
        to_many.update({typ: {}})
        not_converted.update({typ: {}})
        for orig_id, bigg_id in consistent_highest.get(typ).items():
            if len(bigg_id[1]) == 1:
                to_one.get(typ).update({orig_id: [bigg_id[0], bigg_id[1][0]]})
            if len(bigg_id[1]) > 1:
                to_many.get(typ).update({orig_id: bigg_id})
            if len(bigg_id[1]) == 0:
                not_converted.get(typ).update({orig_id: bigg_id})
    return to_one, to_many, not_converted


def checkManyToOne(model_types: [str], to_one: dict):
    one_to_one = {}
    many_to_one = {}
    for typ in model_types:
        reversed_to_one = {bigg_id[1]: [] for bigg_id in to_one.get(typ).values()}
        for (k, v) in to_one.get(typ).items(): reversed_to_one[v[1]].append(k)
        one_to_one.update(
            {typ: {orig_ids[0]: [to_one.get(typ).get(orig_ids[0])[0], conv_id] for conv_id, orig_ids in
                   reversed_to_one.items() if len(orig_ids) == 1}})
        many_to_one.update({typ: {}})
        for value in reversed_to_one.values():
            if len(value) > 1:
                for val in value:
                    many_to_one.get(typ).update({val: to_one.get(typ).get(val)})
    return one_to_one, many_to_one


def checkManyToMany(model_types: [str], to_many: dict):
    to_many_str = {}
    one_to_many = {}
    many_to_many = {}
    for typ in model_types:
        to_many_str.update(
            {typ: {orig_id: [bigg_ids[0], " ".join(sorted(bigg_ids[1]))] for orig_id, bigg_ids in
                   to_many.get(typ).items()}})
    one_to_many_str, many_to_many_str = checkManyToOne(model_types, to_many_str)
    for t in model_types:
        one_to_many.update({t: {key: [value[0], value[1].split()] for key, value in one_to_many_str.get(t).items()}})
        many_to_many.update({t: {key: [value[0], value[1].split()] for key, value in many_to_many_str.get(t).items()}})
    return one_to_many, many_to_many


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


def trySeveralMetabolites(reaction_id: str, bigg_met1, bigg_met2, met1_to_many, met2_to_many, compart1, compart2, BiGG_network_r):
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
            bigg_met1.append(first_confidence_conversion[met1][1])
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
            bigg_met2.append(first_confidence_conversion[met2][1])
        if one_to_many_conversion:
            if met2 in one_to_many_conversion:
                met2_to_many.update({met2: one_to_many_conversion[met2]})
    if (len(bigg_met1) == len(orig_met1)) & (len(bigg_met2) == len(orig_met2)):
        bigg_r = findReactionViaPureEquation(bigg_met1, bigg_met2, BiGG_network_r,
                                                 "Found_via_pure_reaction_equation")
        if not bigg_r:
            bigg_r = workWithH(bigg_met1, bigg_met2, compart1, compart2, BiGG_network_r)
            if not bigg_r:
                bigg_r = {"NOT_found": "Not_found_in_BiGG_network"}
            # TODO: If not found try 3) change c and e compartment to c and p or e and p
    elif (len(bigg_met1) + len(met1_to_many.keys()) == len(orig_met1)) & (
            len(bigg_met2) + len(met2_to_many.keys()) == len(orig_met2)):
        bigg_r = trySeveralMetabolites(reaction_id, bigg_met1, bigg_met2, met1_to_many, met2_to_many, compart1, compart2, BiGG_network_r)
        if not bigg_r:
            bigg_r = {"NOT_found": "Not_found_via_one_to_many_metabolites"}
    else:
        bigg_r = {"NOT_found": "Not_all_metabolites_for_reaction_are_converted"}
    return bigg_r



def selectBasedOnConversionQuality(model_types: [str], converted_obj: dict, obj_type: "metabolites" or "reactions",
                                   models_same_db: dict):
    highest = getHighestConversion(model_types, converted_obj)
    consist, not_consistent = checkSameConversion(models_same_db, highest, obj_type)
    consistent = copy.deepcopy(highest)
    for same_models in models_same_db.values():
        for model in same_models:
            consistent[model] = consist[model]
    to_one, to_many, not_converted = checkOneToManyHighest(model_types, consistent)
    one_to_one, many_to_one = checkManyToOne(model_types, to_one)
    one_to_many, many_to_many = checkManyToMany(model_types, to_many)
    return {"one_to_one": one_to_one, "one_to_many": one_to_many, "many_to_one": many_to_one, "many_to_many": many_to_many,
                                          "not_converted": not_converted, "not_consistent": not_consistent,
            "intermediate_data": {"highest": highest, "consistent": consist, "to_one": to_one, "to_many": to_many}}


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
                # structural_conversion.get(typ).update({key: {}})
                if typ in additional.keys():
                    for orig_id, select_bigg_id in additional.get(typ).items():
                        structural_bigg_id = convertReactionViaNetworkStructure(orig_id, all_models.get(typ),
                                                                                one_one_m, highest_m,
                                                                                bigg_network.get(obj_type), one_to_many_m)
                        structural_conversion.get(typ).update({orig_id: [select_bigg_id, structural_bigg_id, {"first_selection_type": key}]})
                        if (len(structural_bigg_id.keys()) == 1) & (list(structural_bigg_id.keys())[0] != "NOT_found"):
                            structural_conversion_to_one.get(typ).update({orig_id: [select_bigg_id[0], list(structural_bigg_id.keys())]})
                        # if structural_bigg_id[1] in ["Not_found_in_BiGG_network",
                        #                              "Not_all_metabolites_for_reaction_are_converted"]:
                        #     functionA(orig_id, select_bigg_id[1], all_models.get(typ), first_confidence_converted_m,
                        #               bigg_network.get("reactions"))

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
                # structural_conversion.get(typ).update({key: {}})
                if typ in additional.keys():
                    for orig_id, select_bigg_id in additional.get(typ).items():
                        structural_bigg_id = convertReactionViaNetworkStructure(orig_id, all_models.get(typ),
                                                                                final_metabolites.get(typ), highest_m,
                                                                                bigg_network.get(obj_type))
                        structural_conversion.get(typ).update({orig_id: [select_bigg_id, structural_bigg_id, {"first_selection_type": key}]})
                        if (len(structural_bigg_id.keys()) == 1) & (list(structural_bigg_id.keys())[0] != "NOT_found"):
                            structural_conversion_to_one.get(typ).update({orig_id: [select_bigg_id[0], list(structural_bigg_id.keys())]})
                        # if structural_bigg_id[1] in ["Not_found_in_BiGG_network",
                        #                              "Not_all_metabolites_for_reaction_are_converted"]:
                        #     functionA(orig_id, select_bigg_id[1], all_models.get(typ), first_confidence_converted_m,
                        #               bigg_network.get("reactions"))

    return structural_conversion, structural_conversion_to_one
