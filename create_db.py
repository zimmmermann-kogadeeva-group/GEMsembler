#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path
import sqlite3


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
    reactions = (
        bigg_database_r["reaction_string"]
        .str.replace(r"(\d+\.\d*(e-)?\d*|\d+e-\d*)|\+", "", regex=True)
        .apply(lambda x: f" {x} ")
    )

    uniq_all_met = (
        reactions.str.replace(r"(<->)", "", regex=True).str.split().explode().unique()
    )

    m_connections = pd.DataFrame(
        [
            (x, bigg_database_r["bigg_id"][reactions.str.contains(f" {x} ")].tolist())
            for x in uniq_all_met
        ],
        columns=["metabolite", "reactions"],
    ).pipe(apply, "reactions", lambda x: ",".join(x))

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

    return r_connections_no_mix, m_connections


def main(tsv_dir, sql_file):
    # Read in the tables from tsv files
    tsv_dir = Path(tsv_dir)

    # SEEDmodel
    seed_orig = pd.read_csv(tsv_dir / "seed_to_bigg.tsv.gz", sep="\t").rename(
        columns={"seed_ids": "old", "bigg_ids": "new"}
    )

    seed_orig_m = (
        seed_orig.query("type == 'm'").drop(columns="type").reset_index(drop=True)
    )
    seed_orig_r = (
        seed_orig.query("type == 'r'").drop(columns="type").reset_index(drop=True)
    )

    # SEEDmodel additional
    seed_addit = pd.read_csv(tsv_dir / "seed_to_bigg_metanetx.tsv.gz", sep="\t").rename(
        columns={"seed_ids": "old", "bigg_ids": "new"}
    )

    seed_addit_m = (
        seed_addit.query("type == 'm'").drop(columns="type").reset_index(drop=True)
    )
    seed_addit_r = (
        seed_addit.query("type == 'r'").drop(columns="type").reset_index(drop=True)
    )

    # KEGG to BIGG
    kegg_bigg = pd.read_csv(tsv_dir / "kegg_to_bigg_metanetx.tsv.gz", sep="\t").rename(
        columns={"kegg_ids": "old", "bigg_ids": "new"}
    )

    kegg_bigg_m = (
        kegg_bigg.query("type == 'm'").drop(columns="type").reset_index(drop=True)
    )
    kegg_bigg_r = (
        kegg_bigg.query("type == 'r'").drop(columns="type").reset_index(drop=True)
    )

    # Old to new BIGG
    old_new_bigg = pd.read_csv(tsv_dir / "old_to_new_bigg.tsv.gz", sep="\t").rename(
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
    bigg_all_m = pd.read_csv(tsv_dir / "bigg_models_metabolites.tsv.gz", sep="\t")

    bigg_all_r = pd.read_csv(tsv_dir / "bigg_models_reactions.tsv.gz", sep="\t")
    bigg_all_r["universal_bigg_id"] = bigg_all_r["bigg_id"]

    # Put all the tables in one sqlite database
    con = sqlite3.connect(sql_file)
    cur = con.cursor()

    seed_orig_m.to_sql("seed_orig_m", con)
    seed_orig_r.to_sql("seed_orig_r", con)

    seed_addit_m.to_sql("seed_addit_m", con)
    seed_addit_r.to_sql("seed_addit_r", con)

    kegg_bigg_m.to_sql("kegg_bigg_m", con)
    kegg_bigg_r.to_sql("kegg_bigg_r", con)

    old_new_bigg_m.to_sql("old_new_bigg_m", con)
    old_new_bigg_r.to_sql("old_new_bigg_r", con)

    bigg_all_m.to_sql("bigg_all_m", con)
    bigg_all_r.to_sql("bigg_all_r", con)

    bigg_db_network_r, bigg_db_network_m = getBiGGnetwork(bigg_all_r)
    bigg_db_network_r.to_sql("bigg_db_network_r", con)
    bigg_db_network_m.to_sql("bigg_db_network_m", con)

    con.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tsv_dir", help="Path to all tsv files (default=data)", default="data"
    )
    parser.add_argument(
        "--sql_file",
        help="Path to output sqlite file (default=data/all.sqlite",
        default="data/all.sqlite",
    )
    args = parser.parse_args()

    main(**vars(args))
