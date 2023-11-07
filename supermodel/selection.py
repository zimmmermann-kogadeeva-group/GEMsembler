import copy
import pandas as pd
from collections import defaultdict


class Selected(object):
    def __init__(
        self,
        highest: list,
        consistent: list,
        compartments: list,
        replace_with_consistent: bool,
    ):
        self.converted = True
        self.to_one_id = None
        self.from_one_id = None
        self.compartments = compartments
        if replace_with_consistent:
            self.highest_consistent = consistent
        else:
            self.highest_consistent = highest
        if not consistent:
            if not highest:
                self.converted = False
                self.consistent = f"Not converted"
            elif replace_with_consistent:
                self.converted = False
                self.consistent = f"No: {' '.join(highest)}"
            else:
                self.consistent = f"Can be considered as No"
        elif sorted(consistent) == sorted(highest):
            self.consistent = "Yes"
        elif replace_with_consistent:
            self.consistent = (
                f"Changed: from {', '.join(highest)} to {', '.join(consistent)}"
            )
        else:
            self.consistent = f"Could be Changed: from {', '.join(highest)} to {', '.join(consistent)}"
            if not highest:
                self.converted = False
        if len(self.highest_consistent) == 1:
            self.to_one_id = True
        if len(self.highest_consistent) > 1:
            self.to_one_id = False


def checkDBConsistency(
    models_same_db: dict,
    converted_model: dict,
    attr_to_check: str,
    replace_with_consistent=True,
):
    """
    Checking ID for different models with IDs from the same database. If one
    original ID from one database is converted separately for different models,
    looking for intersection of those lists. If no intersection found,
    conversion is inconsistent. If intersection is less that whole converted
    list, remove ids outside intersection. 
    """
    consistent = defaultdict(dict)
    for bd, models in models_same_db.items():
        bd_ids = []
        # no consistency check
        if len(set(list(models.values()))) <= 1 or bd == "bigg":
            for m in models.keys():
                consistent[m] = {
                    idd: Selected(
                        getattr(converted_model.get(m).get(idd), attr_to_check,),
                        getattr(converted_model.get(m).get(idd), attr_to_check,),
                        converted_model.get(m).get(idd).compartments,
                        replace_with_consistent,
                    )
                    for idd in converted_model.get(m).keys()
                }
        # consistency check
        else:
            # collecting all original ids
            for model in models.keys():
                bd_ids = bd_ids + list(converted_model[model].keys())
            for iid in list(set(bd_ids)):
                bigg_ids = []
                mod_present = []
                for mod in models.keys():
                    if converted_model.get(mod).get(iid) is not None:
                        # collecting in models ids is present
                        mod_present.append(mod)
                        if getattr(converted_model.get(mod).get(iid), attr_to_check,):
                            bigg_ids.append(
                                getattr(
                                    converted_model.get(mod).get(iid), attr_to_check,
                                )
                            )
                if not bigg_ids:
                    for pres in mod_present:
                        consistent[pres].update(
                            {
                                iid: Selected(
                                    [],
                                    [],
                                    converted_model.get(pres).get(iid).compartments,
                                    replace_with_consistent,
                                )
                            }
                        )
                else:
                    common_ids = list(set.intersection(*map(set, bigg_ids)))
                    for present in mod_present:
                        consistent[present].update(
                            {
                                iid: Selected(
                                    getattr(
                                        converted_model.get(present).get(iid),
                                        attr_to_check,
                                    ),
                                    common_ids,
                                    converted_model.get(present).get(iid).compartments,
                                    replace_with_consistent,
                                )
                            }
                        )
    return consistent


def checkFromOneFromMany(selected: dict):
    """ Selecting converted with one or several original IDs giving the same converted results: from_one/from_many. """
    to_smth_str = {
        orig_id: " ".join(sorted(bigg_ids.highest_consistent))
        for orig_id, bigg_ids in selected.items()
    }
    reversed_to_smth = {bigg_id: [] for bigg_id in to_smth_str.values()}
    for (k, v) in to_smth_str.items():
        reversed_to_smth[v].append(k)
    for value in reversed_to_smth.values():
        if len(value) == 1:
            selected[value[0]].from_one_id = True
        else:
            for val in value:
                selected[val].from_one_id = False


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
