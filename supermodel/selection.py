from collections import defaultdict


class Selected(object):
    """ Class for selected object that checks results of different stages. After conversion stage
    replace_with_consistent = True and results, going further (highest_consistent), have to be consistent. After
    structural stage replace_with_consistent = False and results, going further (highest_consistent), don't change (
    highest). How results of consistency check differ from original is written in consistency attr. Compartments are
    just passed through for later. To_one_id shows amount of result ids (highest_consistent): True if 1,
    False if more than 1 and None if there is highest_consistent is empty. From_one_id is set as None and will become
    True or False in checkFromOneFromMany function. From_many_other_ids is set as empty list and will be populated
    later in checkFromOneFromMany function. """

    def __init__(
        self,
        highest: list,
        consistent: list,
        compartments: list,
        replace_with_consistent: bool,
        other_present=None,
        other_highest=None,
    ):
        self.to_one_id = None
        self.from_one_id = None
        self.from_many_other_ids = []
        self.compartments = compartments
        if other_present is None:
            self.in_other_models = {}
        else:
            self.in_other_models = {other: [] for other in other_present}
        if replace_with_consistent:
            self.highest_consistent = consistent
        else:
            self.highest_consistent = highest
        if not consistent:
            if not highest:
                self.consistent = f"Not converted"
            elif replace_with_consistent:
                self.consistent = (
                    f"Not consistent. Originally bigg ids were {' '.join(highest)}. But in other models "
                    f"bigg ids were {other_highest}"
                )
            else:
                self.consistent = (
                    f"Can be considered as not consistent. Originally bigg ids were {' '.join(highest)}. "
                    f"But in other models bigg ids were {other_highest}"
                )
        elif set(consistent) == set(highest):
            self.consistent = "Yes"
        elif replace_with_consistent:
            self.consistent = (
                f"Changed: from {', '.join(highest)} to {', '.join(consistent)}. In other models bigg "
                f"ids were {other_highest}"
            )
        else:
            self.consistent = (
                f"Could be Changed: from {', '.join(highest)} to {', '.join(consistent)}. In other "
                f"models bigg ids were {other_highest}"
            )
        if len(self.highest_consistent) == 1:
            self.to_one_id = True
        if len(self.highest_consistent) > 1:
            self.to_one_id = False

    # adding results of mapping for selected objects with the same original id, but from other models
    def check_others(self, original_id, selected: dict):
        for other_model in self.in_other_models.keys():
            self.in_other_models[other_model] = [
                selected[other_model][original_id].to_one_id,
                selected[other_model][original_id].from_one_id,
            ]


def checkDBConsistency(
    models_same_db: dict,
    converted_model: dict,
    attr_to_check: str,
    replace_with_consistent: bool,
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
            # looping per original id through all models
            for iid in list(set(bd_ids)):
                # list for intersecting present ids from models with different db #[[b1, b2], [b1]]
                bigg_ids = []
                # dict connecting ids and models for making .in_other_models attr
                bigg_ids_dict = {}
                mod_present = []
                for mod in models.keys():
                    if converted_model.get(mod).get(iid) is not None:
                        # collecting in models ids is present
                        mod_present.append(mod)
                        # adding only not empty lists because intersecting #[[b1, b2], []] isn't reasonable
                        if getattr(converted_model.get(mod).get(iid), attr_to_check,):
                            bigg_ids.append(
                                getattr(
                                    converted_model.get(mod).get(iid), attr_to_check,
                                )
                            )
                            bigg_ids_dict[mod] = getattr(
                                converted_model.get(mod).get(iid), attr_to_check
                            )
                # if list for intersection is totally empty no intersection is needed
                if not bigg_ids:
                    for pres in mod_present:
                        consistent[pres].update(
                            {
                                iid: Selected(
                                    [],
                                    [],
                                    converted_model.get(pres).get(iid).compartments,
                                    replace_with_consistent,
                                    set(mod_present) - {pres},
                                )
                            }
                        )
                # intersecting result ids for the same original id from models with different db
                else:
                    common_ids = list(set.intersection(*map(set, bigg_ids)))
                    for present in mod_present:
                        other_ids = {
                            k: v for k, v in bigg_ids_dict.items() if k != present
                        }
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
                                    set(mod_present) - {present},
                                    other_ids,
                                )
                            }
                        )
    return consistent


def checkFromOneFromMany(selected: dict):
    """ Checking in selected objects which original id where converted uniquely to one or several bigg ids
    (from_one_id = True) and which original ids have the same results (from_one_id = False).
    Writing down original ids with the same results in from_many_other_ids."""
    # making dict with original ids and bigg results ids as str instead of list
    to_smth_str = {
        orig_id: " ".join(sorted(bigg_ids.highest_consistent))
        for orig_id, bigg_ids in selected.items()
    }
    # creating empty dict with bigg results as keys
    reversed_to_smth = {bigg_id: [] for bigg_id in to_smth_str.values()}
    # populating bigg keys in the empty dict with original ids, for which these bigg ids are results
    for (k, v) in to_smth_str.items():
        reversed_to_smth[v].append(k)
    # writing check results depending on amount of original ids for bigg results
    for value in reversed_to_smth.values():
        if len(value) == 1:
            selected[value[0]].from_one_id = True
        else:
            for val in value:
                selected[val].from_one_id = False
                selected[val].from_many_other_ids = list(set(value) - {val})


def run_selection(
    same_db_models: dict,
    previous_stage: dict,
    attr_to_check: str,
    replace_with_consistent=True,
):
    """Wrapping function to run checking consistency in resulting bigg ids for the same original id but different models
    with different DB. Then checking whether for particular model there is only one original id gives particular results
    or there are several original ids with the same results. And checking how the same original id performs in different
    models, can it give 1 to 1 conversion somewhere else.If there is a group of original ids with the same results,
    connecting them together. """
    first_stage_selected = checkDBConsistency(
        same_db_models, previous_stage, attr_to_check, replace_with_consistent
    )
    for s in first_stage_selected.values():
        checkFromOneFromMany(s)
    for ss in first_stage_selected.values():
        for or_id, sel in ss.items():
            sel.check_others(or_id, first_stage_selected)
    return first_stage_selected


'''

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

'''
