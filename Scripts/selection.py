import copy
import pandas as pd
import cobra
import BiGGnetwork


def getHighestConversion(model_types: [str], all_converted: dict):
    highest_conversion_level = {}
    for typ in model_types:
        highest_conversion_level.update({typ: {}})
        for id, conv in all_converted.get(typ).items():
            for level, converted in conv[1].items():
                if len(converted) > 0:
                    highest_conversion_level.get(typ).update({id: converted})
                    break
                if level == "6-NOconv":
                    highest_conversion_level.get(typ).update({id: []})
    return highest_conversion_level


def checkSameConversion(models_same_db: dict, highest_converted: dict, obj_type: "metabolites" or "reactions",
                        useroutname=None):
    # TODO: check why uncosistent reaction in modelseed was found via structural info and for gapseq not
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
    for bd, models in models_same_db.items():
        bd_ids = []
        for model in models:
            consistent_highest.update({model: {}})
            not_consistent_highest.update({model: {}})
            bd_ids = bd_ids + list(highest_converted.get(model).keys())
        bd_ids = list(set(bd_ids))
        for id in bd_ids:
            bigg_ids = []
            types_to_change = []
            types_present = []
            for typ in models:
                if highest_converted.get(typ).get(id) is not None:
                    types_present.append(typ)
                    if highest_converted.get(typ).get(id) != []:
                        bigg_ids.append(highest_converted.get(typ).get(id))
                    if len(highest_converted.get(typ).get(id)) != 1:
                        types_to_change.append(typ)
            if bigg_ids == []:
                for pres in types_present:
                    consistent_highest.get(pres).update({id: []})
            else:
                common_ids = list(set.intersection(*map(set, bigg_ids)))
                if (len(common_ids) == 1) & (types_to_change != []):
                    for t in types_to_change:
                        changed = pd.concat([changed,
                                             pd.DataFrame([[bd, t, id, common_ids, highest_converted.get(t).get(id)]],
                                                          columns=changed.columns)], ignore_index=True)
                if len(common_ids) > 0:
                    for present in types_present:
                        consistent_highest.get(present).update({id: common_ids})
                else:
                    for pr in types_present:
                        not_consistent_highest.get(pr).update({id: highest_converted.get(pr).get(id)})
                        none_consistent_highest = pd.concat(
                            [none_consistent_highest,
                             pd.DataFrame([[bd, pr, id, common_ids, highest_converted.get(pr).get(id)]],
                                          columns=none_consistent_highest.columns)], ignore_index=True)
    changed.to_csv(filename_changed, sep='\t')
    none_consistent_highest.to_csv(filename_notconsistent, sep='\t')
    return consistent_highest, not_consistent_highest


def checkOneToManyHighest(model_types: [str], consistent_highest: dict):
    # TODO: check why numbers for to_many don't fit to conversion
    to_one = {}
    to_many = {}
    not_converted = {}
    for typ in model_types:
        to_one.update({typ: {}})
        to_many.update({typ: {}})
        not_converted.update({typ: {}})
        for orig_id, bigg_id in consistent_highest.get(typ).items():
            if len(bigg_id) == 1:
                to_one.get(typ).update({orig_id: bigg_id[0]})
            if len(bigg_id) > 1:
                to_many.get(typ).update({orig_id: bigg_id})
            if len(bigg_id) == 0:
                not_converted.get(typ).update({orig_id: bigg_id})
    return to_one, to_many, not_converted


def checkManyToOneHighest(model_types: [str], to_one: dict):
    # TODO: check why for gapseq many_to_one is less then duplicated reactions in original model
    one_to_one = {}
    many_to_one = {}
    for typ in model_types:
        reversed_to_one = {bigg_id: [] for bigg_id in to_one.get(typ).values()}
        for (k, v) in to_one.get(typ).items(): reversed_to_one[v].append(k)
        one_to_one.update(
            {typ: {orig_ids[0]: conv_id for conv_id, orig_ids in reversed_to_one.items() if len(orig_ids) == 1}})
        many_to_one.update({typ: {}})
        for key, value in reversed_to_one.items():
            if len(value) > 1:
                for val in value:
                    many_to_one.get(typ).update({val: key})
    return one_to_one, many_to_one


def convertReactionViaNetworkStructure(reaction_id: str, model: cobra.core.model.Model,
                                       first_confidence_conversion: dict,
                                       BiGG_network_r: pd.core.frame.DataFrame):
    orig_met1 = [react.id for react in model.reactions.get_by_id(reaction_id).reactants]
    bigg_met1 = []
    for met1 in orig_met1:
        if met1 in first_confidence_conversion:
            bigg_met1.append(first_confidence_conversion[met1])
    orig_met2 = [pro.id for pro in model.reactions.get_by_id(reaction_id).products]
    bigg_met2 = []
    for met2 in orig_met2:
        if met2 in first_confidence_conversion:
            bigg_met2.append(first_confidence_conversion[met2])
    print(bigg_met1)
    print(orig_met1)
    print(bigg_met2)
    print(orig_met2)
    if (len(bigg_met1) == len(orig_met1)) & (len(bigg_met2) == len(orig_met2)):
        bigg_met1 = " ".join(sorted(bigg_met1))
        bigg_met2 = " ".join(sorted(bigg_met2))
        bigg_r = BiGG_network_r[
            ((BiGG_network_r["1metabolites"] == bigg_met1) & (BiGG_network_r["2metabolites"] == bigg_met2)) | (
                    (BiGG_network_r["1metabolites"] == bigg_met2) & (BiGG_network_r["2metabolites"] == bigg_met1))][
            "reaction"]
        if not bigg_r.empty:
            bigg_r = bigg_r.values[0]
        else:
            bigg_r = "Not_found_in_BiGG_network"
    else:
        bigg_r = "Not_all_metabolites_for_reaction_are_converted"
    return bigg_r


def convertMetaboliteViaNetworkStructure():
    pass


def selectFirstConfidenceConversion(model_types: [str], converted_obj: dict, obj_type: "metabolites" or "reactions",
                                    models_same_db: dict):
    highest = getHighestConversion(model_types, converted_obj)
    consist, not_consistent = checkSameConversion(models_same_db, highest, obj_type)
    consistent = copy.deepcopy(highest)
    for same_models in models_same_db.values():
        for model in same_models:
            consistent[model] = consist[model]
    to_one, to_many, not_converted = checkOneToManyHighest(model_types, consistent)
    one_to_one, many_to_one = checkManyToOneHighest(model_types, to_one)
    return {"first_confidence": {"one_to_one": one_to_one},
            "for_additional_conversion": {"to_many": to_many, "many_to_one": many_to_one,
                                          "not_converted": not_converted, "not_consistent": not_consistent},
            "intermediate_data": {"highest": highest, "consistent": consist, "to_one": to_one}}


def runAdditionalConversion(model_types: [str], first_stage_selected_m: dict, first_stage_selected_r: dict,
                            all_models: dict, bigg_network: dict,
                            obj_type: "metabolites" or "reactions"):
    structural_conversion = {}
    for typ in model_types:
        structural_conversion.update({typ: {}})
        if obj_type == "reactions":
            first_confidence_converted = first_stage_selected_m.get("first_confidence").get(typ)
            for key, additional in first_stage_selected_r.get("for_additional_conversion").items():
                structural_conversion.get(typ).update({key: {}})
                if typ in additional.keys():
                    for orig_id in additional.get(typ).keys():
                        structural_bigg_id = convertReactionViaNetworkStructure(orig_id, all_models.get(typ),
                                                                                first_confidence_converted,
                                                                                bigg_network.get("reactions"))
                        print(structural_bigg_id)
                        structural_conversion.get(typ).get(key).update({orig_id: structural_bigg_id})
    return structural_conversion


def selectOneID():
    pass
