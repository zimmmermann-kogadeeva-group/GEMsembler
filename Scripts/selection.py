import copy
from copy import deepcopy
import pandas as pd


def getHighestConversion(model_types: [str], all_converted: dict):
    """ Getting list of converted ID with highest level of conversion (1 - the heighest). """
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


def checkDBConsistency(models_same_db: dict, highest_converted: dict, obj_type: "metabolites" or "reactions",
                       write_files=True, do_stat=True, useroutname=None):
    """ Checking ID for different models with IDs from the same database. If one original ID from one database
    is converted separately for different models, looking for intersection of those lists. If no intersection found,
     conversion is inconsistent. If intersection is less that whole converted list, remove ids outside intersection. """
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
                            elif (len(common_ids) == 1) & (len(highest_converted.get(present).get(id)[1]) == 0):
                                changed_stat.get(present)["zero_to_one"] = changed_stat.get(present)["zero_to_one"] + 1
                            elif (len(common_ids) > 1) & (len(highest_converted.get(present).get(id)[1]) > 1):
                                changed_stat.get(present)["n_to_n"] = changed_stat.get(present)["n_to_n"] + 1
                            elif (len(common_ids) > 1) & (len(highest_converted.get(present).get(id)[1]) == 0):
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
            print(f"For {key} model {value['n_to_n']} {obj_type} were changed from several ids to smaller amount of ids") #TODO check why numbers heir are inconcistent from run to run
            print(f"For {key} model {value['zero_to_n']} {obj_type} were changed from no found ids to several")
    consistent = copy.deepcopy(highest_converted)
    for same_models in models_same_db.values():
        for model in same_models:
            consistent[model] = consistent_highest[model]
    return consistent, consistent_highest, not_consistent_highest


def checkToOneToMany(model_types: [str], consistent_highest: dict):
    """ Selecting converted uniquely: with one ID in converted list to_one/to_many. """
    to_one = {}
    to_many = {}
    to_none = {}
    for typ in model_types:
        to_one.update({typ: {}})
        to_many.update({typ: {}})
        to_none.update({typ: {}})
        for orig_id, bigg_id in consistent_highest.get(typ).items():
            if len(bigg_id[1]) == 1:
                to_one.get(typ).update({orig_id: bigg_id})
            elif len(bigg_id[1]) > 1:
                to_many.get(typ).update({orig_id: bigg_id})
            else:
                to_none.get(typ).update({orig_id: bigg_id})
    return to_one, to_many, to_none


def checkFromOneFromMany(model_types: [str], to_smth: dict):
    """ Selecting converted with one or several original IDs giving the same converted results: from_one/from_many. """
    from_one = {}
    from_many = {}
    for typ in model_types:
        to_smth_str = {orig_id: [bigg_ids[0], " ".join(sorted(bigg_ids[1]))] for orig_id, bigg_ids in to_smth.get(typ).items()}
        reversed_to_smth = {bigg_id[1]: [] for bigg_id in to_smth_str.values()}
        for (k, v) in to_smth_str.items(): reversed_to_smth[v[1]].append(k)
        from_one.update(
            {typ: {orig_ids[0]: [to_smth.get(typ).get(orig_ids[0])[0], conv_id.split(" ")] for conv_id, orig_ids in
                   reversed_to_smth.items() if len(orig_ids) == 1}})
        from_many.update({typ: {}})
        for value in reversed_to_smth.values():
            if len(value) > 1:
                for val in value:
                    from_many.get(typ).update({val: to_smth.get(typ).get(val)})
    return from_one, from_many



def runSelection(model_types: [str], converted_obj: dict, obj_type: "metabolites" or "reactions",
                 models_same_db: dict):
    """ Running selection for converted original_to_bigg by uniqueness:
     one_to_one, one_to_many, many_to_one, many_to_many. """
    highest = getHighestConversion(model_types, converted_obj)
    consistent, consist, not_consistent = checkDBConsistency(models_same_db, highest, obj_type)
    to_one, to_many, not_converted = checkToOneToMany(model_types, consistent)
    one_to_one, many_to_one = checkFromOneFromMany(model_types, to_one)
    one_to_many, many_to_many = checkFromOneFromMany(model_types, to_many)
    return {"one_to_one": one_to_one, "one_to_many": one_to_many, "many_to_one": many_to_one, "many_to_many": many_to_many,
                                          "not_converted": not_converted, "not_consistent": not_consistent,
            "intermediate_data": {"highest": highest, "consistent": consist, "to_one": to_one, "to_many": to_many}}


def runNotSelectedMet(model_types: [str], final_obj: dict, selected: dict):
    """ Getting finally not selected metabolites"""
    not_selected = {}
    for typ in model_types:
        not_selected.update({typ: {}})
        for selected_type in selected.values():
            if typ in selected_type.keys():
                ids_notsel = list(set(selected_type.get(typ).keys())-set(final_obj.get("one_one_sugg_met").get(typ).keys()))
                if ids_notsel:
                    not_selected.get(typ).update({i: selected_type.get(typ).get(i) for i in ids_notsel})
    return not_selected


def runNotSelectedR(model_types: [str], final_obj: dict, not_consist: dict, not_uniq: dict, r_info: dict, models: dict):
    """ Getting finally not selected reactions
    and adding periplasmic compartment to them if applicable after structural conversion. """
    not_selected = {}
    for typ in model_types:
        not_selected.update({typ: {}})
        if typ in not_uniq: not_selected.get(typ).update(not_uniq.get(typ))
        if typ in not_consist: not_selected.get(typ).update(not_consist.get(typ))
        ids_notsel = list(set(r_info.get(typ).keys()) - set(final_obj.get(typ).keys()) - set(not_selected.get(typ).keys()))
        if ids_notsel:
            for i in ids_notsel:
                if len(models.get(typ).reactions.get_by_id(i).metabolites) > 24:
                    not_selected.get(typ).update({i: [r_info.get(typ).get(i)[0][0], ["Biomass"]]})
                else:
                    if list(r_info.get(typ).get(i)[1].keys())[0] != "NOT_found":
                        p = []
                        for value in list(r_info.get(typ).get(i)[1].values()):
                            if value.startswith("Found_via_adding_periplasmic_compartment"):
                                p.append("p")
                        if p:
                            not_selected.get(typ).update(
                                {i: [r_info.get(typ).get(i)[0][0]+list(set(p)), list(r_info.get(typ).get(i)[1].keys())]})
                        else:
                            not_selected.get(typ).update({i: [r_info.get(typ).get(i)[0][0], list(r_info.get(typ).get(i)[1].keys())]})
                    else:
                        not_selected.get(typ).update({i: r_info.get(typ).get(i)[0]})
    return not_selected