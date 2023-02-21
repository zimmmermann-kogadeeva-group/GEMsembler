import os
from os.path import join
import pandas as pd


def getBiGGnetwork(bigg_database_r: pd.core.frame.DataFrame, leave_from_mixed_directions=True):
    r_connections = pd.DataFrame(columns=["reaction", "1metabolites", "2metabolites", "models_number"])
    r_connections["reaction"] = bigg_database_r["bigg_id"]
    r_connections["models_number"] = bigg_database_r['model_list'].str.split().apply(len)
    reactions = (bigg_database_r["reaction_string"]
                 .str.replace(r"(\d+\.\d*(e-)?\d*|\d+e-\d*)|\+", "", regex=True)
                 .apply(lambda x: f" {x} "))
    r_connections[["1metabolites", "2metabolites"]] = reactions.str.split("<->", n=1, expand=True)
    a = r_connections
    uniq_all_met = (reactions.str.replace(r"(<->)", "", regex=True)
                    .str.split().explode().unique())
    m_connections = pd.DataFrame(
        [(x, bigg_database_r["bigg_id"][reactions.str.contains(f" {x} ")].tolist())
         for x in uniq_all_met],
        columns=["metabolite", "reactions"]
    )
    r_connections["1metabolites"] = r_connections["1metabolites"].str.split().apply(lambda x: " ".join(sorted(x)))
    r_connections["2metabolites"] = r_connections["2metabolites"].str.split().apply(lambda x: " ".join(sorted(x)))
    b = r_connections
    r_connections = r_connections.sort_values(["1metabolites", "2metabolites", "models_number"],
                                              ascending=[True, True, False])
    r_connections_uniq = r_connections.drop_duplicates(["1metabolites", "2metabolites"], keep="first")
    c = r_connections_uniq
    mixed = pd.DataFrame()
    mixed[["reaction", "1metabolites", "2metabolites", "models_number"]] = r_connections_uniq[
        ["reaction", "2metabolites", "1metabolites", "models_number"]]
    mixed = pd.concat([r_connections_uniq, mixed], ignore_index=True)
    mixed = mixed.sort_values(["1metabolites", "2metabolites", "models_number", "reaction"],
                              ascending=[True, True, False, True])
    dupl = mixed[mixed.duplicated(["1metabolites", "2metabolites"], keep=False)]
    d = dupl
    mixed_reactions = list(set(dupl["reaction"].tolist()))
    if leave_from_mixed_directions:
        dupl_drop = dupl.drop_duplicates(["1metabolites", "2metabolites"], keep="first")
        mixed_reactions = list(set(dupl["reaction"].tolist())-set(dupl_drop["reaction"].tolist()))
    r_connections_no_mix = r_connections_uniq[~r_connections_uniq["reaction"].isin(mixed_reactions)]
    r_connections_no_mix["equation"] = r_connections_no_mix[["1metabolites", "2metabolites"]].values.tolist()
    r_connections_no_mix["equation"] = r_connections_no_mix["equation"].apply(lambda x: "<->".join(sorted(x)))
    return {"reactions": r_connections_no_mix, "metabolites": m_connections, "additional": [a, b, c, d]}


if __name__ == '__main__':
    # region Open conversion tables
    fileDir = os.path.dirname(os.path.realpath('__file__'))  # getting directory of the script for pathes to files
    bigg_all_r = pd.read_csv(join(fileDir, "../../../Databases/BiGG/bigg_models_reactions.txt"), sep="\t")
    bigg_network = getBiGGnetwork(bigg_all_r)
