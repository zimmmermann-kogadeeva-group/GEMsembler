import bisect
import os
from os.path import join
import pandas as pd
import general


def removeNumbersAndSort(equation: str):
    equation_list = equation.replace("+", "").split()
    metabolites_list = []
    for element in equation_list:
        if (not general.is_float(element)) and (element != "+"):
            bisect.insort(metabolites_list, element)
    return metabolites_list


def checkForMetabolite(metabolites: list, metabolite: str):
    return metabolite in metabolites

def getBiGGnetwork(bigg_database_r: pd.core.frame.DataFrame):
    r_connections = pd.DataFrame(columns=["reaction", "1metabolites", "2metabolites"])
    r_connections["reaction"] = bigg_database_r["bigg_id"].str.replace("_[cep]", "", regex=True)
    r_connections[["1metabolites", "2metabolites"]] = bigg_database_r['reaction_string'].str.split("<->", 1,
                                                                                                   expand=True)
    r_connections["1metabolites"] = r_connections["1metabolites"].apply(removeNumbersAndSort)
    r_connections["2metabolites"] = r_connections["2metabolites"].apply(removeNumbersAndSort)
    all_met = r_connections["1metabolites"].tolist() + r_connections["2metabolites"].tolist()
    uniq_all_met = list(set([item for sublist in all_met for item in sublist]))
    m_reactions = {}
    for i, met in enumerate(uniq_all_met):
        print(f"{i} out from {len(uniq_all_met)}")
        reactions = r_connections[
            (r_connections["1metabolites"].apply(checkForMetabolite, metabolite=met)) | r_connections[
                "2metabolites"].apply(checkForMetabolite, metabolite=met)]['reaction'].tolist()
        m_reactions[met]=reactions
    m_connections = pd.DataFrame(m_reactions.items(), columns=["metabolite", "reactions"])
    return r_connections, m_connections


if __name__ == '__main__':
    # region Open conversion tables
    fileDir = os.path.dirname(os.path.realpath('__file__'))  # getting directory of the script for pathes to files
    bigg_all_r = pd.read_csv(join(fileDir, "../../../Databases/BiGG/bigg_models_reactions.txt"), sep="\t")
    r_table, m_table = getBiGGnetwork(bigg_all_r)
    print(m_table)
