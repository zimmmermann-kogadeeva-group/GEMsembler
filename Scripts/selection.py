import pandas as pd


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
    none_consistent_highest = pd.DataFrame(
        columns=["Database", "Model", "original_id", "incommon_bigg_id", "individual_bigg_ids"])
    changed = pd.DataFrame(columns=["Database", "Model", "original_id", "incommon_bigg_id", "individual_bigg_ids"])
    for bd, models in models_same_db.items():
        bd_ids = []
        for model in models:
            consistent_highest.update({model: {}})
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
                        none_consistent_highest = pd.concat(
                            [none_consistent_highest,
                             pd.DataFrame([[bd, pr, id, common_ids, highest_converted.get(pr).get(id)]],
                                          columns=none_consistent_highest.columns)], ignore_index=True)
    changed.to_csv(filename_changed, sep='\t')
    none_consistent_highest.to_csv(filename_notconsistent, sep='\t')
    return consistent_highest


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
            if len(bigg_id) != 1:
                to_many.get(typ).update({orig_id: bigg_id})
                not_converted.get(typ).update({orig_id: bigg_id})
    return to_one, to_many, not_converted

def checkManyToOneHighest(to_one: dict):
    pass


def convertReactionViaNetworkStructure():
    pass


def convertMetaboliteViaNetworkStructure():
    pass


def selectFirstConfidenceConversion(highest_converted: dict):
    pass


def runAdditionalConversion():
    pass


def selectOneID():
    pass
