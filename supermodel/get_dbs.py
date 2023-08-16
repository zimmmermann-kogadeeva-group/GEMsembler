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


def getBiGGnetwork(
    bigg_database_r: pd.core.frame.DataFrame, leave_from_mixed_directions=True
):
    """Getting tables with reaction ids unique reaction equations with
    metabolites sorted for the whole BiGG database. If reaction equations are
    duplicated one used in most amount of models selected."""

    reaction_regex = re.compile(r"(\d+\.\d*(e-)?\d*|\d+e-\d*)|\+|<->")

    m_connections = (
        bigg_database_r[["bigg_id", "reaction_string"]]
        .copy()
        # Remove +, <-> and num. moles in chem. eq. to get list of metabolites
        .pipe(apply, "reaction_string", lambda x: reaction_regex.sub("", x).split())
        # Expand the list of metabolites to one per row
        .explode("reaction_string")
        # Group by the metabolites (column yet to be renamed)
        .groupby("reaction_string")
        # Collect together reaction IDs per each metabolite
        .apply(lambda x: list(x["bigg_id"]))
        .reset_index()
        .rename(columns={0: "reactions", "reaction_string": "metabolite"})
    )

    r_connections = (
        pd.DataFrame(
            {
                "reaction": bigg_database_r["bigg_id"],
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

    return r_connections, m_connections


def get_dbs(path_to_dbs: str):
    path_to_dbs = Path(path_to_dbs)
    # SEEDmodel
    seed_orig = pd.read_csv(path_to_dbs / "seed_to_bigg.tsv.gz", sep="\t").rename(
        columns={"seed_ids": "old", "bigg_ids": "new"}
    )

    seed_orig_m = (
        seed_orig.query("type == 'm'").drop(columns="type").reset_index(drop=True)
    )
    seed_orig_r = (
        seed_orig.query("type == 'r'").drop(columns="type").reset_index(drop=True)
    )

    # SEEDmodel additional
    seed_addit = pd.read_csv(
        path_to_dbs / "seed_to_bigg_metanetx.tsv.gz", sep="\t"
    ).rename(columns={"seed_ids": "old", "bigg_ids": "new"})

    seed_addit_m = (
        seed_addit.query("type == 'm'").drop(columns="type").reset_index(drop=True)
    )
    seed_addit_r = (
        seed_addit.query("type == 'r'").drop(columns="type").reset_index(drop=True)
    )

    # KEGG to BIGG
    kegg_bigg = pd.read_csv(
        path_to_dbs / "kegg_to_bigg_metanetx.tsv.gz", sep="\t"
    ).rename(columns={"kegg_ids": "old", "bigg_ids": "new"})

    kegg_bigg_m = (
        kegg_bigg.query("type == 'm'").drop(columns="type").reset_index(drop=True)
    )
    kegg_bigg_r = (
        kegg_bigg.query("type == 'r'").drop(columns="type").reset_index(drop=True)
    )

    # Old to new BIGG
    old_new_bigg = pd.read_csv(path_to_dbs / "old_to_new_bigg.tsv.gz", sep="\t").rename(
        columns={"old_bigg_ids": "old", "bigg_ids": "new"}
    )

    old_new_bigg_m = (
        old_new_bigg.query("type == 'm'").drop(columns="type").reset_index(drop=True)
    )
    old_new_bigg_r = (
        old_new_bigg.query("type == 'r'").drop(columns="type").reset_index(drop=True)
    )
    old_new_bigg_m["new"] = old_new_bigg_m["new"].str[:-2]

    # BIGG
    bigg_all_m = pd.read_csv(path_to_dbs / "bigg_models_metabolites.tsv.gz", sep="\t")

    bigg_all_r = pd.read_csv(path_to_dbs / "bigg_models_reactions.tsv.gz", sep="\t")
    bigg_all_r["universal_bigg_id"] = bigg_all_r["bigg_id"]

    bigg_db_network_r, bigg_db_network_m = getBiGGnetwork(bigg_all_r)

    return {
        "seed_orig_m": seed_orig_m,
        "seed_orig_r": seed_orig_r,
        "seed_addit_m": seed_addit_m,
        "seed_addit_r": seed_addit_r,
        "kegg_bigg_m": kegg_bigg_m,
        "kegg_bigg_r": kegg_bigg_r,
        "bigg_all_m": bigg_all_m,
        "bigg_all_r": bigg_all_r,
        "bigg_db_network_m": bigg_db_network_m,
        "bigg_db_network_r": bigg_db_network_r,
    }
