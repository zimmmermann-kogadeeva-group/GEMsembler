from collections import defaultdict


class Selected(object):
    """
    Class for selected object that checks results of different stages. After
    conversion stage replace_with_consistent = True and results, going further
    (highest_consistent), have to be consistent. After structural stage
    replace_with_consistent = False and results, going further
    (highest_consistent), don't change (highest). How results of consistency
    check differ from original is written in consistency attr. Compartments are
    just passed through for later. To_one_id shows amount of result ids
    (highest_consistent): True if 1, False if more than 1 and None if there is
    highest_consistent is empty. From_one_id is set as None and will become
    True or False in checkFromOneFromMany function. From_many_other_ids is set
    as empty list and will be populated later in checkFromOneFromMany function.
    """

    def __init__(
        self,
        compartments: list,
        replace_with_consistent: bool,
        highest: list = None,
        consistent: list = None,
        other_present=None,
        other_highest=None,
    ):
        if highest is None:
            highest = []
        if consistent is None:
            consistent = highest

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

    # adding results of mapping for selected objects with the same original id,
    # but from other models
    def check_others(self, original_id, selected: dict):
        for other_model in self.in_other_models.keys():
            self.in_other_models[other_model] = [
                selected[other_model][original_id].to_one_id,
                selected[other_model][original_id].from_one_id,
            ]


def checkDBConsistency(
    models_same_db: dict,
    converted_models: dict,
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
    for db_type, models in models_same_db.items():
        # models is a dictionary of model_id - model_type
        num_model_types = len(set(models.values()))

        # no consistency check
        if num_model_types == 1 or db_type == "bigg":
            for model_id in models.keys():
                consistent[model_id] = {
                    idd: Selected(
                        species.compartments,
                        replace_with_consistent,
                        getattr(species, attr_to_check),
                    )
                    for idd, species in converted_models.get(model_id).items()
                }
        # consistency check
        else:
            # Get a unique set of ids for metabolites or reactions from all
            # models associated with the same database
            all_species_ids = set.union(
                *[
                    set(species.keys())  # mb or rcts ids
                    for model_id, species in converted_models.items()
                    if model_id in models.keys()
                ]
            )

            # looping through all id through all models with the same db type
            for iid in all_species_ids:
                # list for intersecting present ids from models with different
                # db #[[b1, b2], [b1]]
                bigg_ids = []
                # dict connecting ids and models for making .in_other_models attr
                bigg_ids_dict = {}
                models_present = []
                for model_id in models.keys():
                    species = converted_models.get(model_id).get(iid)
                    if species is not None:
                        species_attr = getattr(species, attr_to_check)

                        # collecting in models ids is present
                        models_present.append(model_id)

                        # adding only not empty lists because intersecting
                        # #[[b1, b2], []] isn't reasonable
                        if species_attr:
                            bigg_ids.append(species_attr)
                            bigg_ids_dict[model_id] = species_attr

                # if list for intersection is totally empty no intersection is needed
                if not bigg_ids:
                    for pres in models_present:
                        consistent[pres].update(
                            {
                                iid: Selected(
                                    converted_models.get(pres).get(iid).compartments,
                                    replace_with_consistent,
                                    other_present=set(models_present) - {pres},
                                )
                            }
                        )
                # intersecting result ids for the same original id from models with different db
                else:
                    common_ids = list(set.intersection(*map(set, bigg_ids)))
                    for present in models_present:
                        other_ids = {
                            k: v for k, v in bigg_ids_dict.items() if k != present
                        }
                        consistent[present].update(
                            {
                                iid: Selected(
                                    converted_models.get(present).get(iid).compartments,
                                    replace_with_consistent,
                                    getattr(
                                        converted_models.get(present).get(iid),
                                        attr_to_check,
                                    ),
                                    common_ids,
                                    set(models_present) - {present},
                                    other_ids,
                                )
                            }
                        )
    return consistent


def checkFromOneFromMany(selected: dict):
    """
    Checking in selected objects which original id where converted uniquely to
    one or several bigg ids (from_one_id = True) and which original ids have
    the same results (from_one_id = False).  Writing down original ids with the
    same results in from_many_other_ids.
    """
    # making dict with original ids and bigg results ids as str instead of list
    to_smth_str = {
        orig_id: " ".join(sorted(bigg_ids.highest_consistent))
        for orig_id, bigg_ids in selected.items()
    }
    # creating empty dict with bigg results as keys
    reversed_to_smth = {bigg_id: [] for bigg_id in to_smth_str.values()}
    # populating bigg keys in the empty dict with original ids, for which these bigg ids are results
    for k, v in to_smth_str.items():
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
    """
    Wrapping function to run checking consistency in resulting bigg ids for the
    same original id but different models with different DB. Then checking
    whether for particular model there is only one original id gives particular
    results or there are several original ids with the same results. And
    checking how the same original id performs in different models, can it give
    1 to 1 conversion somewhere else.If there is a group of original ids with
    the same results, connecting them together.
    """
    first_stage_selected = checkDBConsistency(
        same_db_models, previous_stage, attr_to_check, replace_with_consistent
    )
    for s in first_stage_selected.values():
        checkFromOneFromMany(s)
    for ss in first_stage_selected.values():
        for or_id, sel in ss.items():
            sel.check_others(or_id, first_stage_selected)
    return first_stage_selected
