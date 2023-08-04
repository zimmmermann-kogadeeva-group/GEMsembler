import pandas as pd


# Helper functions for pandas data manipulation
def separate(data, col, into, sep=" ", **kwargs):
    data[into] = data[col].str.split(sep, **kwargs)
    return data


def set_value(data, col, value):
    data[col] = value
    return data


def apply(data, col, func, new_col=None, **kwargs):
    new_col = new_col or col
    data[new_col] = data[col].apply(func, **kwargs)
    return data


# Function to get a unique list of reactions. 
# `leave_from_mixed_directions` ...
def get_reaction_ids(data, leave_from_mixed_directions=True):
    mixed_reactions = data["reaction"].unique().tolist()
    if leave_from_mixed_directions:
        mixed_reactions = set(mixed_reactions) - set(
            data.drop_duplicates(["1metabolites", "2metabolites"], keep="first")[
                "reaction"
            ].tolist()
        )
    return mixed_reactions


def getBiGGnetwork(
    bigg_database_r: pd.core.frame.DataFrame, leave_from_mixed_directions=True
):
    """Getting tables with reaction ids unique reaction equations with
    metabolites sorted for the whole BiGG database. If reaction equations are
    duplicated one used in most amount of models selected."""

    # Get equations for all reactions removing amounts of each compound and
    # padding with spaces for later steps
    reactions = (
        bigg_database_r["reaction_string"]
        .str.replace(r"(\d+\.\d*(e-)?\d*|\d+e-\d*)|\+", "", regex=True)
        .apply(lambda x: f" {x} ")
    )

    # Get unique metabolites from set of all reactions
    uniq_all_met = (
        reactions.str.replace(r"(<->)", "", regex=True).str.split().explode().unique()
    )

    # Make a table mapping each metabolite to list of reactions that contain it
    m_connections = pd.DataFrame(
        [
            (x, bigg_database_r["bigg_id"][reactions.str.contains(f" {x} ")].tolist())
            for x in uniq_all_met
        ],
        columns=["metabolite", "reactions"],
    )

    # Create a new table that contains number of models for each reaction. Add
    # columns for metabolites in each side of the equation and sort the
    # metabolites in each of these columns. Then sort the whole table by both
    # metabolite columns and number of models for each reaction. Sorting will
    # help later.
    r_connections_uniq = (
        pd.DataFrame(
            {
                "reaction": bigg_database_r["bigg_id"],
                "models_number": bigg_database_r["model_list"].str.split().apply(len),
            }
        )
        .pipe(
            set_value,
            ["1metabolites", "2metabolites"],
            reactions.str.split("<->", n=1, expand=True),
        )
        .pipe(apply, "1metabolites", lambda x: " ".join(sorted(x.split())))
        .pipe(apply, "2metabolites", lambda x: " ".join(sorted(x.split())))
        .sort_values(
            ["1metabolites", "2metabolites", "models_number"],
            ascending=[True, True, False],
        )
        .drop_duplicates(["1metabolites", "2metabolites"], keep="first")
        .reset_index(drop=True)
    )

    # Get a list of reaction IDs by making a new table that contains reverse of
    # each reaction ...
    mixed_reactions = (
        r_connections_uniq.copy()
        .rename(
            columns={"1metabolites": "2metabolites", "2metabolites": "1metabolites"}
        )
        .pipe(lambda x: pd.concat([x, r_connections_uniq], ignore_index=True))
        .sort_values(
            ["1metabolites", "2metabolites", "models_number", "reaction"],
            ascending=[True, True, False, True],
        )
        .pipe(lambda x: x[x.duplicated(["1metabolites", "2metabolites"], keep=False)])
        .pipe(get_reaction_ids, leave_from_mixed_directions)
    )

    # Get a subset of the table of reactions based on the above list of
    # reaction IDs. 
    r_connections_no_mix = (
        r_connections_uniq[~r_connections_uniq["reaction"].isin(mixed_reactions)]
        .copy()
        .pipe(
            apply,
            ["1metabolites", "2metabolites"],
            lambda x: "<->".join(sorted(x)),
            "equation",
            axis=1,
        )
        .reset_index(drop=True)
    )

    return {"reactions": r_connections_no_mix, "metabolites": m_connections}
