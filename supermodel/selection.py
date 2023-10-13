import copy
import pandas as pd


class Selected(object):
    def __init__(self, highest: list, consistent: list, compartments: list):
        self.converted = True
        self.to_one_id = None
        self.from_one_id = None
        self.compartments = compartments
        self.highest_consistent = consistent
        if not consistent:
            self.converted = False
            if not highest:
                self.consistent = f"Not converted"
            else:
                self.consistent = f"No: {' '.join(highest)}"
        elif consistent == highest:
            self.consistent = "Yes"
        else:
            self.consistent = f"Changed: {' '.join(list(set(highest)-set(consistent)))}"
        if len(self.highest_consistent) == 1:
            self.to_one_id = True
        if len(self.highest_consistent) > 1:
            self.to_one_id = False


def checkDBConsistency(models_same_db: dict, converted_model: dict, attr_to_check: str):
    """ Checking ID for different models with IDs from the same database. If one original ID from one database
    is converted separately for different models, looking for intersection of those lists. If no intersection found,
     conversion is inconsistent. If intersection is less that whole converted list, remove ids outside intersection. """
    consistent = {}
    for bd, models in models_same_db.items():
        bd_ids = {"metabolites": [], "reactions": []}
        if len(set(list(models.values()))) <= 1 or bd == "bigg":
            for m in models.keys():
                for obj_t in bd_ids.keys():
                    consistent.update(
                        {
                            obj_t: {
                                m: {
                                    idd: Selected(
                                        getattr(
                                            converted_model.get(m).get(obj_t).get(idd),
                                            attr_to_check,
                                        ),
                                        getattr(
                                            converted_model.get(m).get(obj_t).get(idd),
                                            attr_to_check,
                                        ),
                                        converted_model.get(m)
                                        .get(obj_t)
                                        .get(idd)
                                        .compartment,
                                    )
                                    for idd in converted_model.get(obj_t).keys()
                                }
                            }
                        }
                    )
        else:
            for model in models.keys():
                consistent.update({model: {"metabolites": {}, "reactions": {}}})
                bd_ids["metabolites"] = bd_ids["metabolites"] + list(
                    converted_model[model]["metabolites"].keys()
                )
                bd_ids["reactions"] = bd_ids["reactions"] + list(
                    converted_model[model]["reactions"].keys()
                )
            for obj_type, ids in bd_ids.items():
                for iid in list(set(ids)):
                    bigg_ids = []
                    mod_present = []
                    for mod in models.keys():
                        if converted_model.get(mod).get(obj_type).get(iid) is not None:
                            mod_present.append(mod)
                            if getattr(
                                converted_model.get(mod).get(obj_type).get(iid),
                                attr_to_check,
                            ):
                                bigg_ids.append(
                                    getattr(
                                        converted_model.get(mod).get(obj_type).get(iid),
                                        attr_to_check,
                                    )
                                )
                    if not bigg_ids:
                        for pres in mod_present:
                            consistent.get(pres).get(obj_type).update(
                                {
                                    iid: Selected(
                                        [],
                                        [],
                                        converted_model.get(pres)
                                        .get(obj_type)
                                        .get(iid)
                                        .compartment,
                                    )
                                }
                            )
                    else:
                        common_ids = list(set.intersection(*map(set, bigg_ids)))
                        for present in mod_present:
                            consistent.get(present).update(
                                {
                                    iid: Selected(
                                        getattr(
                                            converted_model.get(present)
                                            .get(obj_type)
                                            .get(iid),
                                            attr_to_check,
                                        ),
                                        common_ids,
                                        converted_model.get(present)
                                        .get(obj_type)
                                        .get(iid)
                                        .compartment,
                                    )
                                }
                            )
    return consistent


def checkFromOneFromMany(selected: dict):
    """ Selecting converted with one or several original IDs giving the same converted results: from_one/from_many. """
    for sel_obj_typ in selected.values():
        to_smth_str = {
            orig_id: " ".join(sorted(bigg_ids.highest_consistent))
            for orig_id, bigg_ids in sel_obj_typ.items()
        }
        reversed_to_smth = {bigg_id: [] for bigg_id in to_smth_str.values()}
        for (k, v) in to_smth_str.items():
            reversed_to_smth[v].append(k)
        for value in reversed_to_smth.values():
            if len(value) == 1:
                sel_obj_typ[value[0]].from_one_id = True
            else:
                for val in value:
                    sel_obj_typ[val].from_one_id = False


def runNotSelectedMet(model_types: [str], final_obj: dict, selected: dict):
    """ Getting finally not selected metabolites"""
    not_selected = {}
    for typ in model_types:
        not_selected.update({typ: {}})
        for selected_type in selected.values():
            if typ in selected_type.keys():
                ids_notsel = list(
                    set(selected_type.get(typ).keys())
                    - set(final_obj.get("one_one_sugg_met").get(typ).keys())
                )
                if ids_notsel:
                    not_selected.get(typ).update(
                        {i: selected_type.get(typ).get(i) for i in ids_notsel}
                    )
    return not_selected


def runNotSelectedR(
    model_types: [str],
    final_obj: dict,
    not_consist: dict,
    not_uniq: dict,
    r_info: dict,
    models: dict,
):
    """ Getting finally not selected reactions
    and adding periplasmic compartment to them if applicable after structural conversion. """
    not_selected = {}
    for typ in model_types:
        not_selected.update({typ: {}})
        if typ in not_uniq:
            not_selected.get(typ).update(not_uniq.get(typ))
        if typ in not_consist:
            not_selected.get(typ).update(not_consist.get(typ))
        ids_notsel = list(
            set(r_info.get(typ).keys())
            - set(final_obj.get(typ).keys())
            - set(not_selected.get(typ).keys())
        )
        if ids_notsel:
            for i in ids_notsel:
                if len(models.get(typ).reactions.get_by_id(i).metabolites) > 24:
                    not_selected.get(typ).update(
                        {i: [r_info.get(typ).get(i)[0][0], ["Biomass"]]}
                    )
                else:
                    if list(r_info.get(typ).get(i)[1].keys())[0] != "NOT_found":
                        p = []
                        for value in list(r_info.get(typ).get(i)[1].values()):
                            if value.startswith(
                                "Found_via_adding_periplasmic_compartment"
                            ):
                                p.append("p")
                        if p:
                            not_selected.get(typ).update(
                                {
                                    i: [
                                        r_info.get(typ).get(i)[0][0] + list(set(p)),
                                        list(r_info.get(typ).get(i)[1].keys()),
                                    ]
                                }
                            )
                        else:
                            not_selected.get(typ).update(
                                {
                                    i: [
                                        r_info.get(typ).get(i)[0][0],
                                        list(r_info.get(typ).get(i)[1].keys()),
                                    ]
                                }
                            )
                    else:
                        not_selected.get(typ).update({i: r_info.get(typ).get(i)[0]})
    return not_selected
