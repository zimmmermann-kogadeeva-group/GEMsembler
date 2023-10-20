import pandas as pd
from pathlib import Path
import re


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


def remove_duplicates(data, leave_from_mixed_directions=True):
    duplicated = (
        data.copy()
        .rename(
            columns={"1metabolites": "2metabolites", "2metabolites": "1metabolites"}
        )
        .pipe(lambda x: pd.concat([x, data], ignore_index=True))
        .sort_values(
            ["1metabolites", "2metabolites", "models_number", "reaction"],
            ascending=[True, True, False, True],
        )
        .pipe(lambda x: x[x.duplicated(["1metabolites", "2metabolites"], keep=False)])
    )

    reaction_list = duplicated["reaction"].unique().tolist()
    reaction_list_wo_dupl_metabolites = (
        duplicated.drop_duplicates(["1metabolites", "2metabolites"], keep="first")[
            "reaction"
        ]
        .unique()
        .tolist()
    )
    if leave_from_mixed_directions:
        reaction_list = set(reaction_list) - set(reaction_list_wo_dupl_metabolites)
    return data[~data["reaction"].isin(reaction_list)]


def get_bigg_network(path_to_dbs=None, leave_from_mixed_directions=True):
    """Getting dictionary for BiGG topology: unique reaction equation as key and id as value.
    Maybe writing down a table with reaction ids unique reaction equations with
    metabolites sorted for the whole BiGG database. If reaction equations are
    duplicated one used in most amount of models selected."""

    if not path_to_dbs:
        path_to_dbs = Path(__file__).parent.parent / "data_package"
    bigg_database_r = pd.read_csv(
        path_to_dbs / "bigg_models_reactions.txt.gz", sep="\t"
    )

    reaction_regex = re.compile(r"(\d+\.\d*(e-)?\d*|\d+e-\d*)|\+|<->")

    r_connections = (
        pd.DataFrame(
            {
                "reaction": bigg_database_r["universal_bigg_id"],
                "models_number": bigg_database_r["model_list"].str.split().apply(len),
            }
        )
        # Add two columns for set of metabolites on both sides of the equation
        .pipe(
            set_value,
            ["1metabolites", "2metabolites"],
            (
                bigg_database_r["reaction_string"]
                .str.replace(r"(\d+\.\d*(e-)?\d*|\d+e-\d*)|\+", "", regex=True)
                .apply(lambda x: f" {x} ")
                .str.split("<->", n=1, expand=True)
            ),
        )
        # Sort the lists of metabolites in both columns
        .pipe(apply, "1metabolites", lambda x: " ".join(sorted(x.split())))
        .pipe(apply, "2metabolites", lambda x: " ".join(sorted(x.split())))
        # Sort the dataframe by the metabolites
        .sort_values(
            ["1metabolites", "2metabolites", "models_number"],
            ascending=[True, True, False],
        )
        # Remove duplicates based on metabolites
        .drop_duplicates(["1metabolites", "2metabolites"], keep="first")
        # Remove duplicates based on reactions
        .pipe(remove_duplicates, leave_from_mixed_directions)
        # Add an equation column
        .pipe(
            apply,
            ["1metabolites", "2metabolites"],
            lambda x: "<->".join(sorted(x)),
            "equation",
            axis=1,
        )
        .reset_index(drop=True)
    )

    bigg_r_network = (
        r_connections.groupby("equation")
        .apply(lambda x: x["reaction"].values[0])
        .to_dict()
    )
    return bigg_r_network


def get_db(db_name, path_to_dbs=None):
    """Loading conversion tables for different databases db_name: old_new_bigg_m, old_new_bigg_r, seed_orig_m, seed_orig_r,
     seed_addit_m, seed_addit_r, kegg_bigg_m, kegg_bigg_r. """
    if not path_to_dbs:
        path_to_dbs = Path(__file__).parent.parent / "data_package"
    typ = db_name[-1]
    source = db_name.split("_")[0]

    data_table = (
        pd.read_csv(path_to_dbs / str(db_name[:-2] + ".tsv.gz"), sep="\t")
        .rename(columns={source + "_ids": "old", "bigg_ids": "new"})
        .dropna()
    )

    typ_conv = data_table.query(f"type == '{typ}'")
    typ_conv["new"] = typ_conv["new"].str.split(",", expand=False)
    conv_dict = dict(typ_conv.drop(columns="type").values)

    return conv_dict


def get_BiGG_lists(metabolites: bool, path_to_dbs=None):
    if not path_to_dbs:
        path_to_dbs = Path(__file__).parent.parent / "data_package"
    if metabolites:
        bigg_data = pd.read_csv(
            path_to_dbs / "bigg_models_metabolites.txt.gz", sep="\t"
        )
    else:
        bigg_data = pd.read_csv(path_to_dbs / "bigg_models_reactions.txt.gz", sep="\t")
        bigg_data["universal_bigg_id"] = bigg_data["bigg_id"]
    bigg_list = list(set(bigg_data["universal_bigg_id"]))
    return bigg_list
