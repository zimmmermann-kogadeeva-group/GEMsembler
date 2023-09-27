import itertools
from itertools import groupby, combinations
import operator
from copy import deepcopy
import pandas as pd
import cobra
from collections import Counter
import general
import selection


def getReaction(bigg_met1, bigg_met2, BiGG_network_r, comment):
    """ Find reaction id from reaction's metabolites """

    bigg_met1_str = " ".join(sorted(bigg_met1))
    bigg_met2_str = " ".join(sorted(bigg_met2))
    bigg_equation = "<->".join(sorted([bigg_met1_str, bigg_met2_str]))
    bigg_r = BiGG_network_r[bigg_equation]
    if bigg_r:
        return {bigg_r: comment}
    else:
        return {}


def Hydrogens(
    bigg_met1, bigg_met2, compart1, compart2, BiGG_network_r, addit_comment=""
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
                BiGG_network_r,
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
                        BiGG_network_r,
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
                BiGG_network_r,
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
                BiGG_network_r,
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
                        BiGG_network_r,
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
                BiGG_network_r,
                addit_comment + "Found_via_adding_H_to_products",
            )
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
    return bigg_r


def Periplasmic(
    bigg_met1,
    bigg_met2,
    bigg_orig1,
    bigg_orig2,
    compart1,
    compart2,
    BiGG_network_r,
    addit_comment="",
):
    """ Trying to convert reaction after replacing compartments to periplasmic for all possible set of metabolites """
    bigg_r = {}
    # bigg_r_H = {}
    bigg1_comb = []
    for r1 in range(len(bigg_met1) + 1):
        bigg1_comb = bigg1_comb + list(combinations(bigg_met1, r1))
    bigg2_comb = []
    for r2 in range(len(bigg_met2) + 1):
        bigg2_comb = bigg2_comb + list(combinations(bigg_met2, r2))
    for comb1 in bigg1_comb:
        for comb2 in bigg2_comb:
            bigg1 = list(set(bigg_met1) - set(comb1))
            c1 = list(set([b1[-1] for b1 in bigg1]))
            orig1 = [v for k, v in bigg_orig1.items() if k in comb1]
            if comb1:
                bigg1_p = [m1[:-1] + "p" for m1 in comb1]
                c1.append("p")
            else:
                bigg1_p = []
            bigg2 = list(set(bigg_met2) - set(comb2))
            c2 = list(set([b2[-1] for b2 in bigg2]))
            orig2 = [vv for kk, vv in bigg_orig2.items() if kk in comb2]
            if comb2:
                bigg2_p = [m2[:-1] + "p" for m2 in comb2]
                c2.append("p")
            else:
                bigg2_p = []
            tmp_bigg_r = getReaction(
                bigg1 + bigg1_p,
                bigg2 + bigg2_p,
                BiGG_network_r,
                f"Found_via_adding_periplasmic_compartment-{' '.join(orig1)}-{' '.join(comb1)}-{' '.join(bigg1_p)}-{' '.join(orig2)}-{' '.join(comb2)}-{' '.join(bigg2_p)}-{' '.join(list(set(c1 + c2)))}",
            )
            # tmp_bigg_r_H = Hydrogens(bigg1+bigg1_p, bigg2+bigg2_p, compart1+["p"], compart2+["p"], BiGG_network_r,
            #                          f"Found_via_adding_periplasmic_compartment-{bigg1_p}-{bigg2_p}")
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
            # if tmp_bigg_r_H: bigg_r_H.update(tmp_bigg_r_H)
    # if bigg_r:
    return bigg_r
    # else:
    #     return bigg_r_H


def SeveralMetabolites(
    bigg_met1, bigg_met2, met1_to_many, met2_to_many, compart1, compart2, BiGG_network_r
):
    """ Trying to convert reaction with several variants coming from 1-n metabolites by variants from 1-n metabolites to 1-1 converted metabolites """
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
                tmp_bigg_r = getReaction(
                    bigg_met1_mod,
                    bigg_met2_mod,
                    BiGG_network_r,
                    f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}",
                )
                tmp_bigg_r_H = Hydrogens(
                    bigg_met1_mod,
                    bigg_met2_mod,
                    compart1,
                    compart2,
                    BiGG_network_r,
                    f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}",
                )

                if tmp_bigg_r:
                    bigg_r.update(tmp_bigg_r)
                if tmp_bigg_r_H:
                    bigg_r_H.update(tmp_bigg_r_H)
    elif bigg_to_many_variants1 and (not bigg_to_many_variants2):
        for variant1 in bigg_to_many_variants1:
            bigg_met1_mod = bigg_met1 + list(variant1)
            tmp_bigg_r = getReaction(
                bigg_met1_mod,
                bigg_met2,
                BiGG_network_r,
                f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}",
            )
            tmp_bigg_r_H = Hydrogens(
                bigg_met1_mod,
                bigg_met2,
                compart1,
                compart2,
                BiGG_network_r,
                f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants1)}-{' '.join(variant1)}",
            )
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
            if tmp_bigg_r_H:
                bigg_r_H.update(tmp_bigg_r_H)
    elif bigg_to_many_variants2 and (not bigg_to_many_variants1):
        for variant2 in list(bigg_to_many_variants2):
            bigg_met2_mod = bigg_met2 + list(variant2)
            tmp_bigg_r = getReaction(
                bigg_met1,
                bigg_met2_mod,
                BiGG_network_r,
                f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}",
            )
            tmp_bigg_r_H = Hydrogens(
                bigg_met1,
                bigg_met2_mod,
                compart1,
                compart2,
                BiGG_network_r,
                f"Found_via_one_to_many_metabolites-{' '.join(orig_to_many_variants2)}-{' '.join(variant2)}",
            )
            if tmp_bigg_r:
                bigg_r.update(tmp_bigg_r)
            if tmp_bigg_r_H:
                bigg_r_H.update(tmp_bigg_r_H)
    if bigg_r:
        return bigg_r
    else:
        return bigg_r_H


def convertReactionViaNetworkStructure(
    reaction_id: str,
    model: cobra.core.model.Model,
    first_confidence_conversion: dict,
    highest: dict,
    BiGG_network_r: pd.core.frame.DataFrame,
    do_priplasmic: bool,
    one_to_many_conversion=None,
):
    """ Converting one reaction with several strategies. 1) Try just via it's metabolites if all them are converted 1-1.
     2) If not successful try to add/remove H 3) If not successful try to consider periplasmic compartment 4) If all
     metabolites in the reaction are converted 1-1 or 1-n try different variants from this n options. Writing comments
     with information how it was converted and what happened to metabolites for further suggestions"""
    orig_met1 = [react.id for react in model.reactions.get_by_id(reaction_id).reactants]
    bigg_met1 = []
    compart1 = []
    met1_to_many = {}
    bigg_orig1 = {}
    for met1 in orig_met1:
        compart1.append(highest[met1][0][0])
        if met1 in first_confidence_conversion:
            bigg_met1.append(first_confidence_conversion[met1][1][0])
            bigg_orig1.update({first_confidence_conversion[met1][1][0]: met1})
        if one_to_many_conversion:
            if met1 in one_to_many_conversion:
                met1_to_many.update({met1: one_to_many_conversion[met1]})
    orig_met2 = [pro.id for pro in model.reactions.get_by_id(reaction_id).products]
    bigg_met2 = []
    compart2 = []
    met2_to_many = {}
    bigg_orig2 = {}
    for met2 in orig_met2:
        compart2.append(highest[met2][0][0])
        if met2 in first_confidence_conversion:
            bigg_met2.append(first_confidence_conversion[met2][1][0])
            bigg_orig2.update({first_confidence_conversion[met2][1][0]: met2})
        if one_to_many_conversion:
            if met2 in one_to_many_conversion:
                met2_to_many.update({met2: one_to_many_conversion[met2]})
    if (len(bigg_met1) == len(orig_met1)) & (len(bigg_met2) == len(orig_met2)):
        bigg_r = getReaction(
            bigg_met1, bigg_met2, BiGG_network_r, "Found_via_pure_reaction_equation"
        )
        if not bigg_r:
            bigg_r = Hydrogens(bigg_met1, bigg_met2, compart1, compart2, BiGG_network_r)
            if not bigg_r:
                if do_priplasmic:
                    bigg_r = Periplasmic(
                        bigg_met1,
                        bigg_met2,
                        bigg_orig1,
                        bigg_orig2,
                        compart1,
                        compart2,
                        BiGG_network_r,
                    )
                if not bigg_r:
                    bigg_r = {"NOT_found": "Not_found_in_BiGG_network"}
    elif (len(bigg_met1) + len(met1_to_many.keys()) == len(orig_met1)) & (
        len(bigg_met2) + len(met2_to_many.keys()) == len(orig_met2)
    ):
        bigg_r = SeveralMetabolites(
            bigg_met1,
            bigg_met2,
            met1_to_many,
            met2_to_many,
            compart1,
            compart2,
            BiGG_network_r,
        )
        if not bigg_r:
            bigg_r = {"NOT_found": "Not_found_via_one_to_many_metabolites"}
    else:
        bigg_r = {"NOT_found": "Not_all_metabolites_for_reaction_are_converted"}
    return bigg_r


def runStructuralConversion(
    model_types: [str],
    met_for_struct: dict,
    first_stage_selected_m: dict,
    first_stage_selected_r: dict,
    all_models: dict,
    bigg_network: dict,
    models_wo_periplasmic: list,
    met_for_many_sug=None,
):
    """ Running structural conversion for all reactions. Selection reactions that have only 1 id as result """
    structural_conversion = {}
    structural_conversion_to_one = {}
    for typ in model_types:
        structural_conversion.update({typ: {}})
        structural_conversion_to_one.update({typ: {}})
        one_one_m = met_for_struct.get(typ)
        highest_m = (
            first_stage_selected_m.get("intermediate_data").get("highest").get(typ)
        )
        if met_for_many_sug:
            one_to_many_m = met_for_many_sug.get(typ)
        else:
            one_to_many_m = None
        for key, additional in first_stage_selected_r.items():
            if typ in additional.keys():
                for orig_id, select_bigg_id in additional.get(typ).items():
                    structural_bigg_id = convertReactionViaNetworkStructure(
                        orig_id,
                        all_models.get(typ),
                        one_one_m,
                        highest_m,
                        bigg_network.get("reactions"),
                        typ in models_wo_periplasmic,
                        one_to_many_m,
                    )
                    structural_conversion.get(typ).update(
                        {
                            orig_id: [
                                select_bigg_id,
                                structural_bigg_id,
                                {"first_selection_type": key},
                            ]
                        }
                    )
                    if (len(structural_bigg_id.keys()) == 1) & (
                        list(structural_bigg_id.keys())[0] != "NOT_found"
                    ):
                        if list(structural_bigg_id.values())[0].startswith(
                            "Found_via_adding_periplasmic_compartment"
                        ):
                            structural_conversion_to_one.get(typ).update(
                                {
                                    orig_id: [
                                        list(structural_bigg_id.values())[0]
                                        .split("-")[7]
                                        .split(),
                                        list(structural_bigg_id.keys()),
                                    ]
                                }
                            )
                        else:
                            structural_conversion_to_one.get(typ).update(
                                {
                                    orig_id: [
                                        select_bigg_id[0],
                                        list(structural_bigg_id.keys()),
                                    ]
                                }
                            )
    return structural_conversion, structural_conversion_to_one


def getSuggestionForOneToManyMet(
    model_types: [str], structural_r_one_one: dict, structural_r_all: dict
):
    """ Parse comments from structural conversion successful with 1-n metabolites and get which variant from this n was successful """
    suggestions_m = {}
    suggestions_m_sel = {}
    for typ in model_types:
        suggestions_m.update({typ: {}})
        suggestions_m_sel.update({typ: {}})
        for orig_id, struct_id in structural_r_one_one.get(typ).items():
            struct_id = struct_id[1][0]
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
            suggestions_m.get(typ)[key] = {
                variant: number for variant, number in occurrence.items()
            }
            if len(occurrence.keys()) == 1:
                suggestions_m_sel.get(typ).update({key: list(occurrence.keys())})
    return suggestions_m, suggestions_m_sel


def getSuggestionForManyToOneMet(
    model_types: [str],
    many_to_one: dict,
    first_structural_met: dict,
    models: dict,
    BiGG_network_r: pd.core.frame.DataFrame,
    models_same_db: dict,
):
    """ Getting group of original ids that were converted to the same bigg id. For each metabolite in the group trying
    to convert all its reactions via reaction equation if for reaction all its other metabolites are approved
    (converted 1-1 or from previous suggestions) """
    many_to_one_suggestions = {}
    for models_db in models_same_db.values():
        for model in models_db:
            many_to_one_suggestions.update({model: {}})
            others_models = list(set(models_db) - set(model))
            for other in others_models:
                sug_orig = list(
                    set(many_to_one.get(model).keys())
                    & set(first_structural_met.get(other).keys())
                )
                for i in sug_orig:
                    if many_to_one.get(model).get(i) == first_structural_met.get(
                        other
                    ).get(i):
                        many_to_one_suggestions.get(model).update(
                            {i: many_to_one.get(model).get(i)[1]}
                        )

    by_value = operator.itemgetter(1)
    for typ in model_types:
        if typ not in many_to_one_suggestions.keys():
            many_to_one_suggestions.update({typ: {}})
        grouped_many_one = [
            dict(g)
            for k, g in groupby(
                sorted(many_to_one.get(typ).items(), key=by_value), by_value
            )
        ]
        for group in grouped_many_one:
            found_group = {}
            for orig_id, value in group.items():
                reactions = models.get(typ).metabolites.get_by_id(orig_id).reactions
                found = 0
                n = 0
                for reaction in reactions:
                    orig_met1 = [m1.id for m1 in reaction.reactants]
                    orig_met2 = [m2.id for m2 in reaction.products]
                    orig_met = [m.id for m in reaction.metabolites if m.id != orig_id]
                    if (orig_met != []) & (
                        set(orig_met) <= set(first_structural_met.get(typ).keys())
                    ):
                        n = n + 1
                        bigg_met1 = []
                        for met1 in orig_met1:
                            if met1 == orig_id:
                                bigg_met1.append(value[1][0])
                            else:
                                bigg_met1.append(
                                    first_structural_met.get(typ).get(met1)[1][0]
                                )
                        bigg_met2 = []
                        for met2 in orig_met2:
                            if met2 == orig_id:
                                bigg_met2.append(value[1][0])
                            else:
                                bigg_met2.append(
                                    first_structural_met.get(typ).get(met2)[1][0]
                                )
                        tmp_bigg_r = getReaction(
                            bigg_met1, bigg_met2, BiGG_network_r, "found"
                        )
                        if tmp_bigg_r:
                            found = found + 1
                if n != 0:
                    found_group.update({orig_id: found / n})
                if n == 0:
                    found_group.update({orig_id: "met_not_struct"})
            occur = Counter(list(found_group.values()))
            if (occur[1.0] == 1) & (
                set(occur.keys()) <= set([1.0, 0.0, "met_not_struct"])
            ):
                recommend = general.findKeysByValue(found_group, 1, operator.eq)[0]
                if recommend not in many_to_one_suggestions.get(typ).keys():
                    many_to_one_suggestions.get(typ).update(
                        {recommend: group.get(recommend)[1]}
                    )
    return many_to_one_suggestions


def getSuggestionPeriplasmic(
    model_types: [str],
    structural_r_one_one: dict,
    structural_r_all: dict,
    bigg_network: dict,
    models: dict,
):
    """ Getting dictionary of metabolites that gave reaction equation if their compartment is changed to periplasmic and
     dictionary with corresponding reactions. Also, getting ids for transport reactions for metabolite from original
     compartment to periplasmic. Checking which metabolites are changed to periplasmic in all their reactions (replace)
     and wich in part of their reactions (not_replace - split). """
    met_periplasmic = {}
    react_periplasmic = {}
    for typ in model_types:
        met_periplasmic.update({typ: {}})
        react_periplasmic.update({typ: {}})
        for orig_id, struct_id in structural_r_one_one.get(typ).items():
            struct_id = struct_id[1][0]
            comment = structural_r_all.get(typ).get(orig_id)[1].get(struct_id)
            if comment.startswith("Found_via_adding_periplasmic_compartment"):
                react_periplasmic.get(typ).update({orig_id: {}})
                orig_m1 = comment.split("-")[1].split(" ")
                b_m1 = comment.split("-")[2].split(" ")
                p_m1 = comment.split("-")[3].split(" ")
                orig_m2 = comment.split("-")[4].split(" ")
                b_m2 = comment.split("-")[5].split(" ")
                p_m2 = comment.split("-")[6].split(" ")
                for i in range(len(orig_m1)):
                    if orig_m1[i]:
                        react_periplasmic.get(typ).get(orig_id).update(
                            {orig_m1[i]: p_m1[i]}
                        )
                        if orig_m1[i] not in met_periplasmic.get(typ).keys():
                            tr_p = getReaction(
                                [b_m1[i]],
                                [p_m1[i]],
                                bigg_network.get("reactions"),
                                "transport_r_for_periplasmic",
                            )
                            if tr_p:
                                tr_r_p = list(tr_p.keys())[0]
                            elif p_m1[i][:-1] + "e" in models.get(typ).metabolites:
                                tr_p = getReaction(
                                    [p_m1[i][:-1] + "e"],
                                    [p_m1[i]],
                                    bigg_network.get("reactions"),
                                    "transport_r_for_periplasmic_e",
                                )
                                if tr_p:
                                    tr_r_p = list(tr_p.keys())[0]
                                else:
                                    tr_r_p = "notExistTPforP"
                            else:
                                tr_r_p = "notExistTPforP"
                            met_periplasmic.get(typ).update(
                                {orig_m1[i]: [b_m1[i], p_m1[i], tr_r_p, 1]}
                            )
                        else:
                            met_periplasmic.get(typ).get(orig_m1[i])[3] = (
                                met_periplasmic.get(typ).get(orig_m1[i])[3] + 1
                            )
                for j in range(len(orig_m2)):
                    if orig_m2[j]:
                        react_periplasmic.get(typ).get(orig_id).update(
                            {orig_m2[j]: p_m2[j]}
                        )
                        if orig_m2[j] not in met_periplasmic.get(typ).keys():
                            tr_p = getReaction(
                                [b_m2[j]],
                                [p_m2[j]],
                                bigg_network.get("reactions"),
                                "transport_r_for_periplasmic",
                            )
                            if tr_p:
                                tr_r_p = list(tr_p.keys())[0]
                            elif p_m2[j][:-1] + "e" in models.get(typ).metabolites:
                                tr_p = getReaction(
                                    [p_m2[j][:-1] + "e"],
                                    [p_m2[j]],
                                    bigg_network.get("reactions"),
                                    "transport_r_for_periplasmic_e",
                                )
                                if tr_p:
                                    tr_r_p = list(tr_p.keys())[0]
                                else:
                                    tr_r_p = "notExistTPforP"
                            else:
                                tr_r_p = "notExistTPforP"
                            met_periplasmic.get(typ).update(
                                {orig_m2[j]: [b_m2[j], p_m2[j], tr_r_p, 1]}
                            )
                        else:
                            met_periplasmic.get(typ).get(orig_m2[j])[3] = (
                                met_periplasmic.get(typ).get(orig_m2[j])[3] + 1
                            )
    for t in model_types:
        for orig_m, p_sugg in met_periplasmic.get(t).items():
            if len(models.get(t).metabolites.get_by_id(orig_m).reactions) > p_sugg[3]:
                met_periplasmic.get(t).get(orig_m).append("not_replace")
            else:
                met_periplasmic.get(t).get(orig_m).append("replace")
    return met_periplasmic, react_periplasmic


def completeSuggestions(
    model_types: dict, suggestions: dict, selected: dict, models_same_db: dict
):
    """ Adding  metabolites from the same database but different models for further consistency test. """
    met_same_db = set()
    for same_models in models_same_db.values():
        for model in same_models:
            met_same_db = met_same_db | set(suggestions.get(model).keys())
    complete_suggestions = {}
    for typ in model_types:
        complete_suggestions[typ] = {
            key: [
                selected.get("intermediate_data").get("highest").get(typ).get(key)[0],
                val,
            ]
            for key, val in suggestions.get(typ).items()
        }
        if typ in selected.get("intermediate_data").get("consistent").keys():
            for met_to_add in met_same_db:
                if (met_to_add not in list(suggestions.get(typ).keys())) & (
                    met_to_add
                    in list(
                        selected.get("intermediate_data")
                        .get("consistent")
                        .get(typ)
                        .keys()
                    )
                ):
                    complete_suggestions.get(typ).update(
                        {
                            met_to_add: selected.get("intermediate_data")
                            .get("consistent")
                            .get(typ)
                            .get(met_to_add)
                        }
                    )
    return complete_suggestions


def runSuggestionsMet(
    model_types: [str],
    structural_r_info: dict,
    struct_r_uniq: dict,
    allmet_selected: dict,
    models_same_db: dict,
    models: dict,
    BiGG_network_r: dict,
):
    """ Getting suggestion for metabolite from 1st round of structural conversion for one_to_many and many_to_many
    metabolites. """
    suggestions_one_many_m, suggestions_one_many_m_sel = getSuggestionForOneToManyMet(
        model_types, struct_r_uniq, structural_r_info
    )
    compl_sugg_one_many = completeSuggestions(
        model_types, suggestions_one_many_m_sel, allmet_selected, models_same_db
    )
    (
        m_one_many_sug_consist,
        tmp_m_om_sug_consist,
        m_om_sug_not_consist,
    ) = selection.checkDBConsistency(
        models_same_db,
        compl_sugg_one_many,
        "metabolites",
        write_files=False,
        do_stat=False,
    )
    tmp_first_structural_met = deepcopy(m_one_many_sug_consist)
    for typ in model_types:
        tmp_first_structural_met.get(typ).update(
            allmet_selected.get("one_to_one").get(typ)
        )
        # tmp_first_structural_met.update({typ: {k: [v[0], v[1], {"selection_type": "structural_suggestions_for_one_many"}]
        #                                for k, v in m_one_many_sug_consist[typ].items()}})
        # tmp_first_structural_met.get(typ).update(
        #     {key: [val[0], val[1], {"selection_type": "selection_one_one"}] for key, val in
        #      allmet_selected.get("one_to_one").get(typ).items()})
    (
        first_structural_met,
        first_structural_met_many_one,
    ) = selection.checkFromOneFromMany(model_types, tmp_first_structural_met)
    many_to_one_suggestions = getSuggestionForManyToOneMet(
        model_types,
        allmet_selected.get("many_to_one"),
        first_structural_met,
        models,
        BiGG_network_r.get("reactions"),
        models_same_db,
    )
    compl_sugg_many_one = completeSuggestions(
        model_types, many_to_one_suggestions, allmet_selected, models_same_db
    )
    (
        m_many_one_sug_consist,
        tmp_m_mo_sug_consist,
        m_mo_sug_not_consist,
    ) = selection.checkDBConsistency(
        models_same_db,
        compl_sugg_many_one,
        "metabolites",
        write_files=False,
        do_stat=False,
    )
    tmp_final_met = deepcopy(m_many_one_sug_consist)
    for typ in model_types:
        tmp_final_met.get(typ).update(first_structural_met.get(typ))
        # tmp_final_met.update({typ: {k: [v[0], v[1][0], {"selection_type": "structural_suggestions_for_one_many"}] for k, v
        #                         in m_many_one_sug_consist[typ].items()}})
        # tmp_final_met.get(typ).update(first_structural_met.get(typ))
    final_metabolites, final_metabolites_many_one = selection.checkFromOneFromMany(
        model_types, tmp_final_met
    )
    output = {
        "one_one_sugg_met": final_metabolites,
        "checks": [
            final_metabolites_many_one,
            m_mo_sug_not_consist,
            first_structural_met_many_one,
            m_om_sug_not_consist,
        ],
        "suggestions": [suggestions_one_many_m_sel, many_to_one_suggestions],
    }
    return output


def runStructuralCheck(
    model_types: [str],
    allreact_checked: dict,
    allreact_not_pass: dict,
    all_models: dict,
    bigg_network: dict,
):
    """ Checking reactions equations for models with no conversion need (with BiGG ids originally). Should be in BiGG
    database if not exchange or biomass reaction. """
    allreact_struct_checked = deepcopy(allreact_checked)
    allreact_not_pass_struct = deepcopy(allreact_not_pass)
    for typ in model_types:
        for id_r in list(allreact_not_pass.get(typ).keys()):
            bigg_met1 = [
                met1.id
                for met1 in all_models.get(typ).reactions.get_by_id(id_r).reactants
            ]
            bigg_met2 = [
                met2.id
                for met2 in all_models.get(typ).reactions.get_by_id(id_r).products
            ]
            bigg_r = getReaction(
                bigg_met1,
                bigg_met2,
                bigg_network.get("reactions"),
                "found_structural_for_none_converted",
            )
            if bigg_r:
                allreact_struct_checked.get(typ).update(
                    {
                        id_r: [
                            allreact_not_pass.get(typ).get(id_r)[0],
                            list(bigg_r.keys()),
                        ]
                    }
                )
                allreact_not_pass_struct.get(typ)[id_r] = [
                    allreact_not_pass.get(typ).get(id_r)[0],
                    list(bigg_r.values()),
                ]
            else:
                if (
                    (len(bigg_met1) == 1) & (len(bigg_met2) == 0)
                    | (len(bigg_met1) == 0) & (len(bigg_met2) == 1)
                ) & (id_r.startswith("EX_")):
                    allreact_struct_checked.get(typ).update(
                        {id_r: [allreact_not_pass.get(typ).get(id_r)[0], [id_r]]}
                    )
                    allreact_not_pass_struct.get(typ)[id_r] = [
                        allreact_not_pass.get(typ).get(id_r)[0],
                        "not_found_but_Exchange_reaction",
                    ]
                elif (len(bigg_met1) > 20) | (len(bigg_met2) > 20):
                    allreact_not_pass_struct.get(typ)[id_r] = [
                        allreact_not_pass.get(typ).get(id_r)[0],
                        "Growth_reaction",
                    ]
    return allreact_struct_checked, allreact_not_pass_struct
