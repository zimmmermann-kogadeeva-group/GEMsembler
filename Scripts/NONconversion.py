import pandas as pd
def tryOldtoNew(id_to_check : str, obj_type : "metabolites" or "reactions",
                newBiGG_ids : [str], OldNew_table : pd.core.frame.DataFrame, compartment : str) -> str:
    """# for not metabolites and reactions from model with no need to convert update old ids to new if aplicable"""

    if id_to_check in newBiGG_ids:
        return id_to_check
    else:
        new = OldNew_table[OldNew_table["old"] == id_to_check]["new"]
        if new.empty:
            return id_to_check
        else:
            if obj_type == "metabolites":
                new_id = new.values[0] + "_" + compartment[0]
            if obj_type == "reactions":
                new_id = new.values[0]
            return new_id