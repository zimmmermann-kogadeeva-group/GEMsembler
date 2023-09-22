import os
from os.path import join
from pathlib import Path
import pandas as pd


def getBiGGnetwork(bigg_database_r: pd.core.frame.DataFrame, leave_from_mixed_directions=True):
    """ Getting tables with reaction ids unique reaction equations with metabolites sorted for the whole BiGG database.
     If reaction equations are duplicated one used in most amount of models selected. """
    r_connections = pd.DataFrame(columns=["reaction", "1metabolites", "2metabolites", "models_number"])
    r_connections["reaction"] = bigg_database_r["bigg_id"]
    r_connections["models_number"] = bigg_database_r['model_list'].str.split().apply(len)
    reactions = (bigg_database_r["reaction_string"]
                 .str.replace(r"(\d+\.\d*(e-)?\d*|\d+e-\d*)|\+", "", regex=True)
                 .apply(lambda x: f" {x} "))
    r_connections[["1metabolites", "2metabolites"]] = reactions.str.split("<->", n=1, expand=True)
    r_connections["1metabolites"] = r_connections["1metabolites"].str.split().apply(lambda x: " ".join(sorted(x)))
    r_connections["2metabolites"] = r_connections["2metabolites"].str.split().apply(lambda x: " ".join(sorted(x)))
    r_connections = r_connections.sort_values(["1metabolites", "2metabolites", "models_number"],
                                              ascending=[True, True, False])
    r_connections_uniq = r_connections.drop_duplicates(["1metabolites", "2metabolites"], keep="first")
    mixed = pd.DataFrame()
    mixed[["reaction", "1metabolites", "2metabolites", "models_number"]] = r_connections_uniq[
        ["reaction", "2metabolites", "1metabolites", "models_number"]]
    mixed = pd.concat([r_connections_uniq, mixed], ignore_index=True)
    mixed = mixed.sort_values(["1metabolites", "2metabolites", "models_number", "reaction"],
                              ascending=[True, True, False, True])
    dupl = mixed[mixed.duplicated(["1metabolites", "2metabolites"], keep=False)]
    mixed_reactions = list(set(dupl["reaction"].tolist()))
    if leave_from_mixed_directions:
        dupl_drop = dupl.drop_duplicates(["1metabolites", "2metabolites"], keep="first")
        mixed_reactions = list(set(dupl["reaction"].tolist())-set(dupl_drop["reaction"].tolist()))
    r_connections_no_mix = r_connections_uniq[~r_connections_uniq["reaction"].isin(mixed_reactions)]
    r_connections_no_mix["equation"] = r_connections_no_mix[["1metabolites", "2metabolites"]].values.tolist()
    r_connections_no_mix["equation"] = r_connections_no_mix["equation"].apply(lambda x: "<->".join(sorted(x)))
    return r_connections_no_mix


def get_dbs(path_to_dbs: str):
    path_to_dbs = Path(path_to_dbs)
    # SEEDmodel
    seed_orig = pd.read_csv(path_to_dbs / "seed_to_bigg.tsv.gz", sep="\t").rename(
        columns={"seed_ids": "old", "bigg_ids": "new"})
    seed_orig_m = (seed_orig.query("type == 'm'").drop(columns="type").reset_index(drop=True))
    seed_orig_r = (seed_orig.query("type == 'r'").drop(columns="type").reset_index(drop=True))

    # SEEDmodel additional
    seed_addit = pd.read_csv(
        path_to_dbs / "seed_to_bigg_metanetx.tsv.gz", sep="\t").rename(columns={"seed_ids": "old", "bigg_ids": "new"})

    seed_addit_m = (seed_addit.query("type == 'm'").drop(columns="type").reset_index(drop=True))
    seed_addit_r = (seed_addit.query("type == 'r'").drop(columns="type").reset_index(drop=True))

    # KEGG to BIGG
    kegg_bigg = pd.read_csv(
        path_to_dbs / "kegg_to_bigg_metanetx.tsv.gz", sep="\t").rename(columns={"kegg_ids": "old", "bigg_ids": "new"})
    kegg_bigg_m = (kegg_bigg.query("type == 'm'").drop(columns="type").reset_index(drop=True))
    kegg_bigg_r = (kegg_bigg.query("type == 'r'").drop(columns="type").reset_index(drop=True))

    # Old to new BIGG
    old_new_bigg = pd.read_csv(path_to_dbs / "old_to_new_bigg.tsv.gz", sep="\t").rename(
        columns={"old_bigg_ids": "old", "bigg_ids": "new"})
    old_new_bigg_m = (old_new_bigg.query("type == 'm'").drop(columns="type").reset_index(drop=True))
    old_new_bigg_r = (old_new_bigg.query("type == 'r'").drop(columns="type").reset_index(drop=True))
    old_new_bigg_m["new"] = old_new_bigg_m["new"].str[:-2]

    # BIGG
    bigg_all_m = pd.read_csv(path_to_dbs / "bigg_models_metabolites.tsv.gz", sep="\t")
    bigg_all_r = pd.read_csv(path_to_dbs / "bigg_models_reactions.tsv.gz", sep="\t")
    bigg_all_r["universal_bigg_id"] = bigg_all_r["bigg_id"]

    bigg_db_network_r = getBiGGnetwork(bigg_all_r)

    return {
        "seed_orig_m": seed_orig_m,
        "seed_orig_r": seed_orig_r,
        "seed_addit_m": seed_addit_m,
        "seed_addit_r": seed_addit_r,
        "kegg_bigg_m": kegg_bigg_m,
        "kegg_bigg_r": kegg_bigg_r,
        "bigg_all_m": bigg_all_m,
        "bigg_all_r": bigg_all_r,
        "old_new_bigg_m": old_new_bigg_m,
        "old_new_bigg_r": old_new_bigg_r,
        "bigg_db_network_r": bigg_db_network_r,
    }


if __name__ == '__main__':
    # region Open conversion tables
    fileDir = os.path.dirname(os.path.realpath('__file__'))  # getting directory of the script for pathes to files
    bigg_all_r = pd.read_csv(join(fileDir, "../../../Databases/BiGG/bigg_models_reactions.txt"), sep="\t")
    bigg_network = getBiGGnetwork(bigg_all_r)
