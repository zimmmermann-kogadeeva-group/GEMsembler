import copy
import operator

import pandas as pd
from collections import Counter
import general


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


def convertReactionViaNetworkStructure():
    pass


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
    return {"first_confidence": one_to_one,
            "for_additional_conversion": {"to_many": to_many, "many_to_one": many_to_one,
                                          "not_converted": not_converted, "not_consistent": not_consistent},
            "intermediate_data": [highest, consist, to_one]}


def runAdditionalConversion():
    pass


def selectOneID():
    pass
