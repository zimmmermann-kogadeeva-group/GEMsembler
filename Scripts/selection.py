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
                    highest_conversion_level.get(typ).update({id: [conv[0], converted]})
                    break
                if level == "6-NOconv":
                    highest_conversion_level.get(typ).update({id: [conv[0], []]})
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
    changed.to_csv(filename_changed, sep='\t')
    for key, value in changed_stat.items():  # TODO: write this stat in file
        print(f"For {key} model in total {value['total_change']} {obj_type} were changed")
        print(f"For {key} model {value['n_to_one']} {obj_type} were changed from several ids to one")
        print(f"For {key} model {value['zero_to_one']} {obj_type} were changed from no found ids to one")
        print(f"For {key} model {value['n_to_n']} {obj_type} were changed from several ids to smaller amount of ids")
        print(f"For {key} model {value['zero_to_n']} {obj_type} were changed from no found ids to several")
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


def convertReactionViaNetworkStructure(reaction_id: str, model: cobra.core.model.Model,
                                       first_confidence_conversion: dict, one_to_many_conversion: dict, consistent: dict,
                                       BiGG_network_r: pd.core.frame.DataFrame):
    orig_met1 = [react.id for react in model.reactions.get_by_id(reaction_id).reactants]
    bigg_met1 = []
    compart1 = []
    met1_to_many = {}
    for met1 in orig_met1:
        # compart1.append(consistent[met1][0][0])
        if met1 in first_confidence_conversion:
            bigg_met1.append(first_confidence_conversion[met1][1])
        elif met1 in one_to_many_conversion:
            met1_to_many.update({met1: one_to_many_conversion[met1]})
    orig_met2 = [pro.id for pro in model.reactions.get_by_id(reaction_id).products]
    bigg_met2 = []
    compart2 = []
    met2_to_many = {}
    for met2 in orig_met2:
        if met2 in first_confidence_conversion:
            # compart2.append(consistent[met2][0][0])
            bigg_met2.append(first_confidence_conversion[met2][1])
        elif met2 in one_to_many_conversion:
            met2_to_many.update({met2: one_to_many_conversion[met2]})
    if (len(bigg_met1) == len(orig_met1)) & (len(bigg_met2) == len(orig_met2)):
        bigg_met1 = " ".join(sorted(bigg_met1))
        bigg_met2 = " ".join(sorted(bigg_met2))
        bigg_r = BiGG_network_r[
            ((BiGG_network_r["1metabolites"] == bigg_met1) & (BiGG_network_r["2metabolites"] == bigg_met2)) | (
                    (BiGG_network_r["1metabolites"] == bigg_met2) & (BiGG_network_r["2metabolites"] == bigg_met1))][
            "reaction"]
        if not bigg_r.empty:
            bigg_r = bigg_r.values[0]
            comment = "Found_via_pure_reaction_equation"
        else:
            # for c1 in compart1:
            #     bigg_met1_mod = bigg_met1.split()
            #     if "h_" + c1 in bigg_met1_mod:
            #         bigg_met1_mod.remove("h_"+c1)
            #     else:
            #         bigg_met1_mod.append("h_"+c1)
            bigg_r = None
            comment = "Not_found_in_BiGG_network"
            # TODO: If not found try 1) add h 2) remove h 3) change c and e compartment to c and p or e and p
            # TODO: Write modifications made at the begging/end of id
    # elif (len(bigg_met1) + len(met1_to_many.keys()) == len(orig_met1)) & (
    #     len(bigg_met2) + len(met2_to_many.keys()) == len(orig_met2)):
    #     bigg_variants = []
    #     if met1_to_many != {}:
    #         for key1, val1 in met1_to_many.items():
    #             for met_variant1 in val1[1]:
    #                 bigg_met1.append(met_variant1)
    #                 print("a")
    #                 bigg_met1 = " ".join(sorted(bigg_met1))
    #                 if met2_to_many != {}:
    #                     for key2, val2 in met2_to_many.items():
    #                         for met_variant2 in val2:
    #                             bigg_met2.append(met_variant2)
    #                             bigg_met2 = " ".join(sorted(bigg_met2))
    #                             bigg_r = BiGG_network_r[
    #                                 ((BiGG_network_r["1metabolites"] == bigg_met1) & (
    #                                         BiGG_network_r["2metabolites"] == bigg_met2)) | (
    #                                         (BiGG_network_r["1metabolites"] == bigg_met2) & (
    #                                         BiGG_network_r["2metabolites"] == bigg_met1))][
    #                                 "reaction"]
    #                             if not bigg_r.empty:
    #                                 bigg_r = bigg_r.values[0]
    #                                 comment = f"Found_via_one_to_many_metabolites: "
    else:
        bigg_r = None
        comment = "Not_all_metabolites_for_reaction_are_converted"
        # TODO: If some metabolites are not converted with first confidence try options from one_to_many for metabolites
        # TODO: If successful, write suggestions for metabolite conversion according to particular reaction (connect m and r)
    return [bigg_r, comment]


# TODO: write functionA that will work if previous (conversion by only structure) didn't work out
#   It will work for all except not_converted to check available conversion variants and select the best
#   Compere reactions equations for all options + look for conversion possibilities for metabolites => come up with strategy

def functionA(reaction_id: str, bigg_r_ids: [str], model: cobra.core.model.Model,
              first_confidence_conversion: dict,
              BiGG_network_r: pd.core.frame.DataFrame):
    print(reaction_id)
    orig_met1 = [react.id for react in model.reactions.get_by_id(reaction_id).reactants]
    bigg_met1 = []
    for met1 in orig_met1:
        if met1 in first_confidence_conversion:
            bigg_met1.append(first_confidence_conversion[met1][1])
        else:
            bigg_met1.append(met1)
    orig_met2 = [pro.id for pro in model.reactions.get_by_id(reaction_id).products]
    bigg_met2 = []
    for met2 in orig_met2:
        if met2 in first_confidence_conversion:
            bigg_met2.append(first_confidence_conversion[met2][1])
        else:
            bigg_met2.append(met2)
    bigg_met1 = " ".join(bigg_met1)
    bigg_met2 = " ".join(bigg_met2)
    print(bigg_met1)
    print(bigg_met2)
    for r in bigg_r_ids:
        bigg_equation = BiGG_network_r[BiGG_network_r["reaction"] == r]
        print(bigg_equation)


def convertMetaboliteViaNetworkStructure():
    pass
    # for each reaction in reactions for metabolites check structure in bigg_network and conversion for other metabolites


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
    return {"first_confidence": {"one_to_one": one_to_one},
            "for_additional_conversion": {"one_to_many": one_to_many, "many_to_one": many_to_one,
                                          "many_to_many": many_to_many,
                                          "not_converted": not_converted, "not_consistent": not_consistent},
            "intermediate_data": {"highest": highest, "consistent": consist, "to_one": to_one, "to_many": to_many}}


def runAdditionalConversion(model_types: [str], first_stage_selected_m: dict, first_stage_selected_r: dict,
                            all_models: dict, bigg_network: dict,
                            obj_type: "metabolites" or "reactions"):
    structural_conversion = {}
    for typ in model_types:
        structural_conversion.update({typ: {}})
        if obj_type == "reactions":
            structural_conversion.get(typ).update({"one_to_one": {}})
            first_confidence_converted_m = first_stage_selected_m.get("first_confidence").get("one_to_one").get(typ)
            first_confidence_converted_r = first_stage_selected_r.get("first_confidence").get("one_to_one").get(typ)
            one_to_many_m = first_stage_selected_m.get("for_additional_conversion").get("one_to_many").get(typ)
            consist_m = first_stage_selected_m.get("intermediate_data").get("consistent").get(typ)
            # TODO write final checking of r first confidence: if not converted or not found leave original; if structural id is found and different take structural
            for orig_first_r in first_confidence_converted_r.keys():
                struct_first_r = convertReactionViaNetworkStructure(orig_first_r, all_models.get(typ),
                                                                    first_confidence_converted_m, one_to_many_m,
                                                                    consist_m,
                                                                    bigg_network.get("reactions"))
                structural_conversion.get(typ).get("one_to_one").update(
                    {orig_first_r: [first_confidence_converted_r[orig_first_r], struct_first_r]})
            for key, additional in first_stage_selected_r.get("for_additional_conversion").items():
                structural_conversion.get(typ).update({key: {}})
                if typ in additional.keys():
                    for orig_id, select_bigg_id in additional.get(typ).items():
                        structural_bigg_id = convertReactionViaNetworkStructure(orig_id, all_models.get(typ),
                                                                                first_confidence_converted_m,
                                                                                one_to_many_m, consist_m,
                                                                                bigg_network.get("reactions"))
                        structural_conversion.get(typ).get(key).update({orig_id: [select_bigg_id, structural_bigg_id]})
                        if structural_bigg_id[1] in ["Not_found_in_BiGG_network",
                                                     "Not_all_metabolites_for_reaction_are_converted"]:
                            functionA(orig_id, select_bigg_id[1], all_models.get(typ), first_confidence_converted_m,
                                      bigg_network.get("reactions"))
            # TODO: check for not found and many_to_one problem after convertReactionViaNetworkStructure functiom
            #   and pass all this to functionA

        # TODO: write for metabolites
        #   1) one_to_many with suggestions 2) others
    return structural_conversion


def selectOneID():
    pass
