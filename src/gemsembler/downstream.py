import operator
import re
import warnings
from collections import defaultdict
from copy import deepcopy
from pathlib import Path

import dill
import networkx as nx
import pandas as pd
from cobra.flux_analysis import pfba
from .creation import SuperModel
from .comparison import getCoreGPR, getCoreConnections
import seaborn as sns
from .drawing import (
    draw_one_known_pathway,
    draw_one_synt_path,
    get_pyvis_from_nx,
    define_edge_features,
    define_node_features,
    get_color_palette,
    MET_NOT_INT_GLOBAL,
)


def write_metabolites_production_output(
    out_bp_production_tab,
    write_output_to_folder: str,
    file_name=None,
    yticklabels=True,
    cmap="mako",
    metric="jaccard",
    method="single",
    center=0.3,
    linewidths=0.003,
    **kwargs,
):
    """Function to plot heatmap of produced or not produced metabolites.
    Parameters: table with production, output fold and
    attributes of seaborn clustermap function."""
    if file_name is None:
        out_bp_production_tab.to_csv(
            f"{write_output_to_folder}/all_metabolites_production.tsv",
            sep="\t",
            index=False,
        )
    else:
        out_bp_production_tab.to_csv(
            f"{write_output_to_folder}/{file_name}.tsv", sep="\t", index=False,
        )
    out_bp_production_tab.index = out_bp_production_tab["Metabolites"]
    num = out_bp_production_tab[
        out_bp_production_tab.columns.difference(["Metabolites"])
    ]
    num_no_grow = num.drop(index="overall_growth", errors="ignore")
    num_no_grow = num_no_grow.astype(
        {mod_name: "float" for mod_name in num_no_grow.columns}
    )
    bp_heatmap = sns.clustermap(
        num_no_grow,
        yticklabels=yticklabels,
        metric=metric,
        method=method,
        cmap=cmap,
        center=center,
        linewidths=linewidths,
        **kwargs,
    )
    if file_name is None:
        bp_heatmap.savefig(f"{write_output_to_folder}/all_metabolites_production.png")
    else:
        bp_heatmap.savefig(f"{write_output_to_folder}/{file_name}" f"")


def table_reactions_confidence(
    supermodel: SuperModel,
    output_name: str,
    pathway_r=None,
    path_r_dist=None,
    yes_range=1,
    no_range=1,
    genes=True,
    and_as_solid=False,
    add_original_models=True,
    write_table=True,
):
    if pathway_r is None:
        pathway_r = list(supermodel.reactions.assembly.keys())
    output = {"Reaction": pathway_r, "R_confidence": [], "Status": [], "R_models": []}
    if path_r_dist is not None:
        if len(set(pathway_r) & set(list(path_r_dist.keys()))) != len(pathway_r):
            raise ValueError(
                f"Not all reactions ids from pathway are present"
                f" in the distance dictionary. "
                f"Please check pathway_r and path_r_dist input"
            )
        output.update({f"Distance from {path_r_dist['metabolite']}": []})
    if genes:
        output.update({"GPR_confidence": [], "GPR_core": [], "GPR_assembly": []})
    if add_original_models:
        for source in supermodel.sources:
            output.update({source: []})
    for r_id in pathway_r:
        if r_id in supermodel.reactions.assembly.keys():
            r = supermodel.reactions.assembly[r_id]
            output["R_confidence"].append(f"{r.in_models['models_amount']}")
            if r.in_models["models_amount"] >= len(supermodel.sources) - yes_range:
                output["Status"].append("yes")
            elif r.in_models["models_amount"] <= no_range:
                output["Status"].append("no")
            else:
                output["Status"].append("q")
            output["R_models"].append(" ".join(sorted(r.in_models["models_list"])))
            if genes:
                if not r.gene_reaction_rule["assembly"]:
                    output["GPR_confidence"].append("0")
                    output["GPR_core"].append("")
                    output["GPR_assembly"].append("")
                else:
                    for i in range(r.in_models["models_amount"], 0, -1):
                        gpr_core = getCoreGPR(
                            r.gene_reaction_rule,
                            i,
                            operator.ge,
                            r.in_models["models_list"],
                            and_as_solid,
                        )
                        if gpr_core:
                            output["GPR_confidence"].append(f"{i}")
                            output["GPR_core"].append(gpr_core[0])
                            output["GPR_assembly"].append(
                                r.gene_reaction_rule["assembly"][0]
                            )
                            break
            if add_original_models:
                for source in supermodel.sources:
                    if source in r.in_models["models_list"]:
                        output[source].append("+")
                    else:
                        output[source].append("-")
            if path_r_dist is not None:
                output[f"Distance from {path_r_dist['metabolite']}"].append(
                    path_r_dist[r_id]
                )
        else:
            output["R_confidence"].append("0")
            output["Status"].append("no")
            output["R_models"].append("")
            if genes:
                output["GPR_confidence"].append("0")
                output["GPR_core"].append("")
                output["GPR_assembly"].append("")
            if add_original_models:
                for source in supermodel.sources:
                    output[source].append("-")
            if path_r_dist is not None:
                output[f"Distance from {path_r_dist['metabolite']}"].append(
                    path_r_dist[r_id]
                )
    output_tb = pd.DataFrame(output)
    if write_table:
        output_tb.to_csv(output_name, sep="\t", index=False)
    return output_tb


def pathway_of_interest(
    supermodel: SuperModel,
    pathway_r: dict or list,
    draw_pathway_to_file=None,
    write_table_to_file=None,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
    yes_range=1,
    no_range=1,
):
    if type(pathway_r) == list:
        if met_not_int is None:
            met_not_int = deepcopy(MET_NOT_INT_GLOBAL)
        pathway_r_dict = {}
        for r_id in pathway_r:
            pathway_r_dict[r_id] = []
            r = supermodel.reactions.assembly[r_id]
            for rea in r.reactants["assembly"]:
                for pro in r.products["assembly"]:
                    if (re.sub("_([cep])$", "", rea.id) not in met_not_int.keys()) and (
                        re.sub("_([cep])$", "", pro.id) not in met_not_int.keys()
                    ):
                        pathway_r_dict[r_id].append((rea.id, pro.id))
        pathway_r_list = pathway_r
    else:
        pathway_r_dict = pathway_r
        pathway_r_list = list(pathway_r.keys())
    if draw_pathway_to_file is not None:
        g = draw_one_known_pathway(
            supermodel,
            pathway_r_dict,
            draw_pathway_to_file,
            additional_met,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
    if write_table_to_file is not None:
        t = table_reactions_confidence(
            supermodel,
            write_table_to_file,
            pathway_r_list,
            yes_range=yes_range,
            no_range=no_range,
            genes=genes,
            and_as_solid=and_as_solid,
        )
    if draw_pathway_to_file and write_table_to_file:
        return g, t
    elif draw_pathway_to_file and not write_table_to_file:
        return g
    elif not draw_pathway_to_file and write_table_to_file:
        return t


def biosynthesis_pathway_with_media_and_metabolite(
    supermodel: SuperModel,
    pathway: list,
    medium: list,
    met_to_synt: str,
    draw_pathway_to_file=None,
    write_pathway_table_to_file=None,
    calc_dist_from_synt_met=True,
    check_distance=5,
    **kwargs,
):
    if draw_pathway_to_file is not None:
        g = draw_one_synt_path(
            supermodel, pathway, medium, [met_to_synt], draw_pathway_to_file, **kwargs
        )
    if write_pathway_table_to_file is not None:
        if calc_dist_from_synt_met:
            r_dist_dict = calc_dist_for_synt_path(
                pathway, met_to_synt, supermodel, check_distance
            )
            t = table_reactions_confidence(
                supermodel, write_pathway_table_to_file, pathway, r_dist_dict, **kwargs
            )
        else:
            t = table_reactions_confidence(
                supermodel, write_pathway_table_to_file, pathway, **kwargs
            )
    if draw_pathway_to_file and write_pathway_table_to_file:
        return g, t
    elif draw_pathway_to_file and not write_pathway_table_to_file:
        return g
    elif not draw_pathway_to_file and write_pathway_table_to_file:
        return t


def glycolysis(
    supermodel: SuperModel,
    draw_pathway_to_file=None,
    write_table_to_file=None,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
    yes_range=1,
    no_range=1,
):
    glycolysis = {
        "HEX1": [("glc__D_c", "g6p_c")],
        "PGI": [("g6p_c", "f6p_c")],
        "PFK": [("f6p_c", "fdp_c")],
        "FBA": [("fdp_c", "g3p_c")],
        "GAPD": [("g3p_c", "13dpg_c")],
        "PGK": [("13dpg_c", "3pg_c")],
        "PGM": [("3pg_c", "2pg_c")],
        "ENO": [("2pg_c", "pep_c")],
        "PYK": [("pep_c", "pyr_c")],
    }
    if draw_pathway_to_file is not None:
        g = draw_one_known_pathway(
            supermodel,
            glycolysis,
            draw_pathway_to_file,
            additional_met,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
    if write_table_to_file:
        t = table_reactions_confidence(
            supermodel,
            write_table_to_file,
            list(glycolysis.keys()),
            yes_range=yes_range,
            no_range=no_range,
            genes=genes,
            and_as_solid=and_as_solid,
        )
    if draw_pathway_to_file and write_table_to_file:
        return g, t
    elif draw_pathway_to_file and not write_table_to_file:
        return g
    elif not draw_pathway_to_file and write_table_to_file:
        return t


def tca(
    supermodel: SuperModel,
    draw_pathway_to_file=None,
    write_table_to_file=None,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
    yes_range=1,
    no_range=1,
):
    tca = {
        "PYK": [("pep_c", "pyr_c")],
        "PPC": [("pep_c", "oaa_c")],
        "PPCK": [("pep_c", "oaa_c")],
        "PEPCK_re": [("pep_c", "oaa_c")],
        "PC": [("pyr_c", "oaa_c")],
        "PDH": [("pyr_c", "accoa_c")],
        "PFL": [("pyr_c", "accoa_c")],
        "CS": [("accoa_c", "cit_c"), ("oaa_c", "cit_c")],
        "ACONT": [("cit_c", "icit_c")],
        "ACONTa": [("cit_c", "acon_C_c")],
        "ACONTb": [("acon_C_c", "icit_c")],
        "ICL": [("icit_c", "succ_c"), ("icit_c", "glx_c")],
        "ICDHyr": [("icit_c", "akg_c")],
        "ICDHx": [("icit_c", "akg_c")],
        "ICITRED": [("icit_c", "osuc_c")],
        "OSUCCL": [("osuc_c", "akg_c")],
        "AKGDH": [("akg_c", "succoa_c")],
        "OOR2r": [("akg_c", "succoa_c"), ("fdxo_42_c", "fdxr_42_c")],
        "AKGDa": [("akg_c", "sdhlam_c"), ("lpam_c", "sdhlam_c")],
        "AKGDb": [("sdhlam_c", "succoa_c"), ("sdhlam_c", "dhlam_c")],
        "PDHcr": [("dhlam_c", "lpam_c")],
        "SUCOAS": [("succoa_c", "succ_c")],
        "SUCDi": [("succ_c", "fum_c")],
        "FRD7": [("fum_c", "succ_c")],
        "FUM": [("fum_c", "mal__L_c")],
        "MALS": [("glx_c", "mal__L_c")],
        "MDH": [("mal__L_c", "oaa_c")],
        "MDH2": [("mal__L_c", "oaa_c")],
        "MDH3": [("mal__L_c", "oaa_c")],
    }
    if draw_pathway_to_file:
        g = draw_one_known_pathway(
            supermodel,
            tca,
            draw_pathway_to_file,
            additional_met,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
    if write_table_to_file:
        t = table_reactions_confidence(
            supermodel,
            write_table_to_file,
            list(tca.keys()),
            yes_range=yes_range,
            no_range=no_range,
            genes=genes,
            and_as_solid=and_as_solid,
        )
    if draw_pathway_to_file and write_table_to_file:
        return g, t
    elif draw_pathway_to_file and not write_table_to_file:
        return g
    elif not draw_pathway_to_file and write_table_to_file:
        return t


def pentose_phosphate(
    supermodel: SuperModel,
    draw_pathway_to_file=None,
    write_table_to_file=None,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
    yes_range=1,
    no_range=1,
):
    pentose_phosphate_pathway = {
        "G6PBDH": [("g6p_B_c", "6pgl_c")],
        "G6PDH2r": [("g6p_c", "6pgl_c")],
        "PGL": [("6pgl_c", "6pgc_c")],
        "GND": [("6pgc_c", "ru5p__D_c")],
        "RPE": [("ru5p__D_c", "xu5p__D_c")],
        "RPI": [("ru5p__D_c", "r5p_c")],
        "TKT1": [("xu5p__D_c", "g3p_c"), ("r5p_c", "s7p_c")],
        "TALA": [("g3p_c", "f6p_c"), ("s7p_c", "e4p_c")],
        "TKT2": [("e4p_c", "f6p_c"), ("xu5p__D_c", "g3p_c")],
        "PGI": [("f6p_c", "g6p_c")],
        "G6PI": [("g6p_c", "g6p_B_c")],
    }
    if draw_pathway_to_file:
        g = draw_one_known_pathway(
            supermodel,
            pentose_phosphate_pathway,
            draw_pathway_to_file,
            additional_met,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
    if write_table_to_file:
        t = table_reactions_confidence(
            supermodel,
            write_table_to_file,
            list(pentose_phosphate_pathway.keys()),
            yes_range=yes_range,
            no_range=no_range,
            genes=genes,
            and_as_solid=and_as_solid,
        )
    if draw_pathway_to_file and write_table_to_file:
        return g, t
    elif draw_pathway_to_file and not write_table_to_file:
        return g
    elif not draw_pathway_to_file and write_table_to_file:
        return t


def biomass(
    supermodel: SuperModel,
    only_difference=False,
    draw_plot_to_file=None,
    write_table_to_file=None,
    directed=True,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
    yes_range=1,
    no_range=1,
):
    if (draw_plot_to_file is not None) and (not draw_plot_to_file.endswith(".html")):
        raise ValueError(
            "Output file for the plot is of wrong format. Please use html file name."
        )
    color_brewer = get_color_palette(len(supermodel.sources))
    if n_letter is None:
        n_letter = supermodel.get_short_name_len()
    if write_table_to_file is not None:
        output = {
            "Metabolite": [],
            "Type": [],
            "Met_confidence": [],
            "Biomass_confidence": [],
            "Status": [],
        }
        for source in supermodel.sources:
            output.update({source: []})
    g = nx.DiGraph()
    biomass_r = supermodel.reactions.assembly.get("Biomass")
    colname_r = define_node_features(
        color_brewer, "single_path_r", biomass_r, n_letter,
    )
    g.add_node(
        colname_r[0],
        label=colname_r[1],
        shape="box",
        color=colname_r[2],
        title=colname_r[3],
    )
    core_rea = []
    if only_difference:
        core_rea = getCoreConnections(
            biomass_r.reactants,
            len(supermodel.sources),
            operator.ge,
            supermodel.sources,
        )
    for rea in biomass_r.reactants.get("assembly"):
        if rea not in core_rea:
            colname_rea = define_node_features(
                color_brewer, "metabolites", rea, n_letter
            )
            rea_edge = define_edge_features(
                color_brewer,
                "single_path_r",
                "metabolites",
                rea,
                colname_rea[1],
                biomass_r,
                colname_r[1],
                "reactants",
                n_letter,
            )
            g.add_node(
                colname_rea[0],
                shape="o",
                label=colname_rea[1],
                color=colname_rea[2],
                title=colname_rea[3],
            )
            for e in rea_edge:
                g.add_edge(
                    e[0], e[1], color=e[2], font_color="black", title=e[3],
                )
            if write_table_to_file is not None:
                output["Metabolite"].append(colname_rea[0])
                output["Type"].append("reactant")
                output["Met_confidence"].append(f"Core{rea.in_models['models_amount']}")
                cr = int(rea_edge[0][3].split(":\n")[0].split(" ")[-1])
                output["Biomass_confidence"].append(f"Core{cr}")
                if cr >= len(supermodel.sources) - yes_range:
                    output["Status"].append("yes")
                elif cr <= no_range:
                    output["Status"].append("no")
                else:
                    output["Status"].append("q")
                short_s_bp = rea_edge[0][3].split(":\n")[1]
                l = int(len(short_s_bp) / cr)
                sources_bp = [
                    short_s_bp[i : i + l] for i in range(0, len(short_s_bp), l)
                ]
                for source in supermodel.sources:
                    if source[:l] in sources_bp:
                        output[source].append("+")
                    else:
                        output[source].append("-")
    core_pro = []
    if only_difference:
        core_pro = getCoreConnections(
            biomass_r.products,
            len(supermodel.sources),
            operator.ge,
            supermodel.sources,
        )
    for pro in biomass_r.products.get("assembly"):
        if pro not in core_pro:
            colname_pro = define_node_features(
                color_brewer, "metabolites", pro, n_letter
            )
            pro_edge = define_edge_features(
                color_brewer,
                "single_path_r",
                "metabolites",
                pro,
                colname_pro[1],
                biomass_r,
                colname_r[1],
                "products",
                n_letter,
            )
            g.add_node(
                colname_pro[0],
                shape="o",
                label=colname_pro[1],
                color=colname_pro[2],
                title=colname_pro[3],
            )
            for e in pro_edge:
                g.add_edge(
                    e[0], e[1], color=e[2], font_color="black", title=e[3],
                )
            if write_table_to_file is not None:
                output["Metabolite"].append(colname_pro[0])
                output["Type"].append("product")
                output["Met_confidence"].append(f"Core{pro.in_models['models_amount']}")
                cp = int(pro_edge[0][3].split(":\n")[0].split(" ")[-1])
                output["Biomass_confidence"].append(f"Core{cp}")
                if cp >= len(supermodel.sources) - yes_range:
                    output["Status"].append("yes")
                elif cp <= no_range:
                    output["Status"].append("no")
                else:
                    output["Status"].append("q")
                short_s_bp = pro_edge[0][3].split(":\n")[1]
                l = int(len(short_s_bp) / cp)
                sources_bp = [
                    short_s_bp[i : i + l] for i in range(0, len(short_s_bp), l)
                ]
                for source in supermodel.sources:
                    if source[:l] in sources_bp:
                        output[source].append("+")
                    else:
                        output[source].append("-")
    if draw_plot_to_file:
        pyvis_graph = get_pyvis_from_nx(g, directed, size, wid, hei)
        pyvis_graph.write_html(draw_plot_to_file, notebook=False)
    if write_table_to_file is not None:
        output_tb = pd.DataFrame(output)
        output_tb.to_csv(write_table_to_file, sep="\t", index=False)
        if draw_plot_to_file:
            return g, output_tb
        else:
            return output_tb
    elif draw_plot_to_file:
        return g


def get_met_neighborhood(
    supermodel: SuperModel,
    metabolite_id: str,
    neighborhood_dist=2,
    highly_connected_t=10,
    draw_neiborhood_to_file=None,
    draw_met_not_int=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=1920,
    hei=1080,
    size=25,
):
    if metabolite_id not in supermodel.metabolites.assembly.keys():
        raise ValueError(f"{metabolite_id} is not in the supermodel")
    if met_not_int is None:
        met_not_int = deepcopy(MET_NOT_INT_GLOBAL)
    all_r = set()
    all_g = set()
    med_high_connect = []
    all_m = {metabolite_id}
    dist_m = [metabolite_id]
    for i in range(neighborhood_dist):
        new_all_m = set()
        for m_id in dist_m:
            tmp_m = re.sub("_([cep])$", "", m_id)
            if (tmp_m in met_not_int.keys()) & (not draw_met_not_int):
                continue
            met = supermodel.metabolites.assembly.get(m_id)
            reactions = [r.id for r in met.reactions["assembly"]]
            if len(reactions) > highly_connected_t:
                reactions = []
                med_high_connect.append(re.sub("_([cep])$", "", met.id))
            all_r = all_r | set(reactions)
            for r_id in reactions:
                if r_id == "Biomass":
                    continue
                new_all_m = new_all_m | set(
                    [
                        m.id
                        for m in supermodel.reactions.assembly[r_id]
                        .metabolites["assembly"]
                        .keys()
                    ]
                )
                all_g = all_g | set(
                    [
                        g.id
                        for g in supermodel.reactions.assembly[r_id].genes["assembly"]
                    ]
                )
        dist_m = new_all_m
        all_m = all_m | new_all_m
    if draw_neiborhood_to_file is not None:
        graph = draw_one_synt_path(
            supermodel,
            list(all_r),
            med_high_connect,
            [metabolite_id],
            draw_neiborhood_to_file,
            draw_met_not_int,
            genes,
            and_as_solid,
            directed,
            met_not_int,
            n_letter,
            wid,
            hei,
            size,
        )
        return graph, all_r, all_g, all_m
    else:
        return all_r, all_g, all_m


def calc_dist_for_synt_path(
    pathway_rs: list, met_of_interest: str, supermodel: SuperModel, check_distance=5,
):
    path_r_dist = {r: ">5" for r in pathway_rs}
    path_r_dist.update({"metabolite": met_of_interest})
    for i in range(1, check_distance + 1):
        all_r, all_g, all_m = get_met_neighborhood(
            supermodel, met_of_interest, neighborhood_dist=i, highly_connected_t=50
        )
        r_intersect = list(set(pathway_rs) & set(all_r))
        for ri in r_intersect:
            if path_r_dist[ri] == ">5":
                path_r_dist[ri] = i
    return path_r_dist


def fba_growth_met_production(
    models_to_analyse: dict,
    medium: dict,
    biomass_precursors=True,
    met_of_interest=None,
    biomass_r_id=None,
):
    medium_no_comp = [re.sub("_([cep])$", "", m) for m in medium.keys()]
    if biomass_r_id is None:
        biomass_r_id = "Biomass"
    if met_of_interest is None and biomass_precursors:
        met_of_interest = []
        for m in models_to_analyse.values():
            met_of_interest = list(
                set(
                    met_of_interest
                    + [
                        react.id
                        for react in m.reactions.get_by_id(biomass_r_id).reactants
                    ]
                )
            )
    out_bp_production = {"Metabolites": ["overall_growth"] + met_of_interest}
    flux_res_out = defaultdict(dict)
    path_pfba_out = defaultdict(dict)
    stat_out = {}
    bp_synt = []
    for k, model in models_to_analyse.items():
        med_mod = {}
        for e in model.exchanges:
            if list(e.metabolites.keys())[0].id in medium.keys():
                med_mod.update({e.id: medium[list(e.metabolites.keys())[0].id]})
            else:
                med_mod.update({e.id: 0})
        old_medium = model.medium
        model.medium = med_mod
        res = model.optimize()
        flux_res_out[k] = {"overall_fba": res}
        model_data = [res.objective_value]
        model_stat = 0
        model_med_stat = 0
        all_met = [mm.id for mm in model.metabolites]
        all_r = [rr.id for rr in model.reactions]
        bp_model = [
            react.id for react in model.reactions.get_by_id(biomass_r_id).reactants
        ]
        for bp in met_of_interest:
            if bp not in all_met:
                model_data.append(-0.5)
                continue
            demand_added = False
            if ("DM_" + bp) not in all_r:
                model.add_boundary(model.metabolites.get_by_id(bp), type="demand")
                demand_added = True
            model.objective = model.demands.get_by_id("DM_" + bp)
            res_bp = model.optimize()
            flux_res_out[k].update({bp + "_fba": res_bp})
            if res_bp.objective_value > 0.001:
                model_data.append(1)
                if re.sub("_([cep])$", "", bp) in medium_no_comp:
                    model_med_stat = model_med_stat + 1
                else:
                    model_stat = model_stat + 1
                if bp not in bp_synt:
                    bp_synt.append(bp)
                pfba_res = pfba(model)
                flux_res_out[k].update({bp + "_pfba": pfba_res})
                reactions = list(
                    pfba_res.to_frame()[
                        (pfba_res.to_frame()["fluxes"] > 0.001)
                        | (pfba_res.to_frame()["fluxes"] < -0.001)
                    ].index
                )
                path_pfba_out[k].update({bp + "_path_pfba": reactions})
            else:
                if (bp not in bp_model) and biomass_precursors:
                    model_data.append(0.5)
                else:
                    model_data.append(0)
            if demand_added:
                model.remove_reactions(["DM_" + bp])
        model.medium = old_medium
        model.objective = model.reactions.get_by_id(biomass_r_id)
        out_bp_production.update({k: model_data})
        stat_out.update({k: model_stat})
        stat_out.update({"medium_" + k: model_med_stat})
    stat_out.update({"not_synthesized": len(set(met_of_interest) - set(bp_synt))})
    return out_bp_production, flux_res_out, path_pfba_out, stat_out


def write_pfba_results(
    path_pfba_out: dict,
    supermodel: SuperModel,
    medium: list,
    write_output_to_folder: str,
    draw_pfba_for_models=None,
    draw_met_not_int=False,
    draw_pfba=True,
    table_pfba=True,
    calc_r_dist=True,
    check_distance=5,
    **kwargs,
):
    if draw_pfba_for_models is None:
        met_model = {}
        for mod, res in path_pfba_out.items():
            for met in res.keys():
                if met in met_model.keys():
                    if mod.startswith("core") and (
                        (not met_model[met][0].startswith("core"))
                        or (
                            int(mod.removeprefix("core"))
                            > int(met_model[met][0].removeprefix("core"))
                        )
                    ):
                        met_model[met] = [mod]
                    elif mod == "assembly" and (
                        not met_model[met][0].startswith("core")
                    ):
                        met_model[met] = [mod]
                    elif (not met_model[met][0].startswith("core")) and (
                        met_model[met][0] != "assembly"
                    ):
                        met_model[met].append(mod)
                else:
                    met_model[met] = [mod]
        pd.DataFrame(
            met_model.items(), columns=["Metabolite", "Most confident pfba"]
        ).to_csv(
            f"{write_output_to_folder}/bp_most_confident_pfba_production.tsv",
            index=False,
            sep="\t",
        )
        for m, mod_conf in met_model.items():
            for model_id in mod_conf:
                v = [
                    vv for vv in path_pfba_out[model_id][m] if not vv.startswith("DM_")
                ]
                if draw_pfba:
                    g = draw_one_synt_path(
                        supermodel,
                        v,
                        medium,
                        [m.removesuffix("_path_pfba")],
                        f"{write_output_to_folder}/{m}_{model_id}.html",
                        draw_met_not_int=draw_met_not_int,
                    )
                if table_pfba:
                    if calc_r_dist:
                        r_dist_dict = calc_dist_for_synt_path(
                            v, m.removesuffix("_path_pfba"), supermodel, check_distance
                        )
                        t = table_reactions_confidence(
                            supermodel,
                            f"{write_output_to_folder}/{m}_{model_id}.tsv",
                            v,
                            r_dist_dict,
                            **kwargs,
                        )

                    else:
                        t = table_reactions_confidence(
                            supermodel,
                            f"{write_output_to_folder}/{m}_{model_id}.tsv",
                            v,
                            **kwargs,
                        )
    else:
        for model_id in draw_pfba_for_models:
            for k, v in path_pfba_out[model_id].items():
                v = [vv for vv in v if not vv.startswith("DM_")]
                if draw_pfba:
                    g = draw_one_synt_path(
                        supermodel,
                        v,
                        medium,
                        [k.removesuffix("_path_pfba")],
                        f"{write_output_to_folder}/{k}_{model_id}.html",
                        draw_met_not_int=draw_met_not_int,
                    )
                if table_pfba:
                    if calc_r_dist:
                        r_dist_dict = calc_dist_for_synt_path(
                            v, k.removesuffix("_path_pfba"), supermodel, check_distance
                        )
                        t = table_reactions_confidence(
                            supermodel,
                            f"{write_output_to_folder}/{k}_{model_id}.tsv",
                            v,
                            r_dist_dict,
                            **kwargs,
                        )

                    else:
                        t = table_reactions_confidence(
                            supermodel,
                            f"{write_output_to_folder}/{k}_{model_id}.tsv",
                            v,
                            **kwargs,
                        )


def run_growth_full_flux_analysis(
    models_to_analyse: dict,
    medium: dict,
    biomass_precursors=True,
    metabolites_of_interest=None,
    write_output_to_folder=None,
    file_name=None,
    draw_pfba_for_models=None,
    supermodel_for_write_pfba=None,
    draw_pfba=True,
    table_pfba=True,
    draw_met_not_int=False,
    biomass_r_id=None,
    **kwargs,
):
    (
        out_bp_production,
        flux_res_out,
        path_pfba_out,
        stat_out,
    ) = fba_growth_met_production(
        models_to_analyse,
        medium,
        biomass_precursors,
        metabolites_of_interest,
        biomass_r_id,
    )
    out_bp_production_tab = pd.DataFrame(out_bp_production)
    stat_out_tab = pd.DataFrame(
        stat_out.items(),
        columns=["Metabolites confidence production", "Metabolites amount"],
    )
    if write_output_to_folder is not None:
        write_metabolites_production_output(
            out_bp_production_tab, write_output_to_folder, file_name, **kwargs
        )
        if file_name is None:
            file_name = write_output_to_folder + "/production_confidence_stat.tsv"
        stat_out_tab.to_csv(file_name, sep="\t", index=False)
    if supermodel_for_write_pfba is not None:
        if write_output_to_folder is None:
            raise ValueError(
                "Output folder is not provided. "
                "Please provide output folder "
                "if you want to draw or write pfba pathways."
            )
        if draw_pfba_for_models is None:
            warnings.warn(
                "Models for which to draw pfba of biomass precursors is not provided. "
                "So for each biomass precursor the model with "
                "the highest confidence level will be used."
            )
        write_pfba_results(
            path_pfba_out,
            supermodel_for_write_pfba,
            list(medium.keys()),
            write_output_to_folder,
            draw_pfba_for_models,
            draw_met_not_int,
            draw_pfba,
            table_pfba,
        )
    return out_bp_production_tab, flux_res_out, path_pfba_out, stat_out_tab


def run_metquest_results_analysis(
    folder_with_mq_res_folders: str,
    model_list: list,
    metabolites_ids: list,
    medium: list,
    supermodel=None,
    cofactors=None,
    output_folder=None,
    output_file_name=None,
    draw_mq_path=False,
    table_mq_path=False,
    calc_dist_path=True,
    check_distance=5,
    output_folder_mq_paths_plots=None,
    output_folder_mq_paths_tables=None,
    **kwargs,
):
    if output_folder_mq_paths_plots is not None:
        draw_mq_path = True
    if output_folder_mq_paths_tables is not None:
        table_mq_path = True
    if (draw_mq_path is True) and (supermodel is None):
        raise ValueError(
            "Supermodel is needed to draw pathways, but not provided. "
            "Please, provide supermodel."
        )
    if (table_mq_path is True) and (supermodel is None):
        raise ValueError(
            "Supermodel is needed to write tables of the pathways, but not provided. "
            "Please, provide supermodel."
        )
    if (draw_mq_path is True) and (output_folder_mq_paths_plots is None):
        if output_folder_mq_paths_tables is not None:
            output_folder_mq_paths_plots = output_folder_mq_paths_tables
        else:
            output_folder_mq_paths_plots = output_folder
    if (table_mq_path is True) and (output_folder_mq_paths_tables is None):
        if output_folder_mq_paths_plots is not None:
            output_folder_mq_paths_tables = output_folder_mq_paths_plots
        else:
            output_folder_mq_paths_tables = output_folder
    if cofactors is None:
        cofactors = [
            "co2_c",
            "hco3_c",
            "pi_c",
            "ppi_c",
            "ACP_c",
            "atp_c",
            "adp_c",
            "amp_c",
            "nad_c",
            "nadh_c",
            "nadp_c",
            "nadph_c",
            "coa_c",
            "cmp_c",
            "cdp_c",
            "ctp_c",
            "gmp_c",
            "gdp_c",
            "gtp_c",
            "ump_c",
            "udp_c",
            "utp_c",
            "fadh2_c",
            "fad_c",
            "q8_c",
            "q8h2_c",
            "mqn8_c",
            "mql8_c",
            "mqn6_c",
            "mql6_c",
            "thf_c",
        ]
    metquest_all_res_paths = {}
    stat_out = {}
    met_synt = []
    met_med = []
    met_cof = []
    folder_with_out_folders = Path(folder_with_mq_res_folders)
    for mod in model_list:
        for f in folder_with_out_folders.iterdir():
            if f.is_dir() and f.name.startswith(mod):
                metquest_all_res_paths[mod] = f / "shortest_paths.pkl"
    synthes_out = {"Metabolites": metabolites_ids} | {mod: [] for mod in model_list}
    for mod in model_list:
        model_stat = 0
        mq_data = dill.load(open(metquest_all_res_paths[mod], "rb"))
        for bp in metabolites_ids:
            if bp not in mq_data["paths"].keys():
                to_append = -0.5
            else:
                if bp in medium:
                    to_append = 0.75
                    if bp not in met_med:
                        met_med.append(bp)
                elif bp in cofactors:
                    to_append = 0.5
                    if bp not in met_cof:
                        met_cof.append(bp)
                elif mq_data["paths"][bp]:
                    to_append = 1
                    model_stat = model_stat + 1
                    if bp not in met_synt:
                        met_synt.append(bp)
                else:
                    to_append = 0
            synthes_out[mod].append(to_append)
            if (to_append == 1) and (draw_mq_path is True):
                allalt = set()
                for altp_len in mq_data["paths"][bp].values():
                    for altp in altp_len:
                        allalt = allalt | altp
                draw_one_synt_path(
                    supermodel,
                    list(allalt),
                    medium,
                    [bp],
                    f"{output_folder_mq_paths_plots}/{bp}_{mod}.html",
                )
            if (to_append == 1) and (table_mq_path is True):
                allalt = set()
                for altp_len in mq_data["paths"][bp].values():
                    for altp in altp_len:
                        allalt = allalt | altp
                if calc_dist_path:
                    r_dist_dict = calc_dist_for_synt_path(
                        list(allalt), bp, supermodel, check_distance
                    )
                    table_reactions_confidence(
                        supermodel,
                        f"{output_folder_mq_paths_tables}/{bp}_{mod}.tsv",
                        list(allalt),
                        r_dist_dict,
                    )
                else:
                    table_reactions_confidence(
                        supermodel,
                        f"{output_folder_mq_paths_tables}/{bp}_{mod}.tsv",
                        list(allalt),
                    )
        stat_out.update({mod: model_stat})
    synthes_tab_out = pd.DataFrame(synthes_out)
    stat_out.update(
        {
            "Metabolites in medium": len(met_med),
            "Metabolites in cofactors": len(met_cof),
            "Metabolites not synthesized": len(
                set(metabolites_ids) - set(met_synt) - set(met_med) - set(met_cof)
            ),
        }
    )
    stat_out_tab = pd.DataFrame(
        stat_out.items(),
        columns=["Metabolites confidence production", "Metabolites amount"],
    )
    if output_folder is not None:
        write_metabolites_production_output(synthes_tab_out, output_folder, **kwargs)
        if output_file_name is None:
            output_file_name = output_folder + "/production_confidence_stat.tsv"
        stat_out_tab.to_csv(output_file_name, sep="\t", index=False)
    return synthes_tab_out, stat_out_tab
