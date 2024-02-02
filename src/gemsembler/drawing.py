import math
import re
import networkx as nx
import pyautogui
import seaborn as sns
from .creation import NewObject, SuperModel
from pyvis.network import Network


def getColorPalette(sources_number: int) -> dict:
    up = math.ceil(sources_number / 10)
    colors = sns.color_palette("bright").as_hex()
    out_colors = {
        "metabolites": sns.light_palette(colors[4], sources_number + up).as_hex()[up:],
        "intersection": sns.light_palette(colors[3], sources_number + up).as_hex()[up:],
        "1stPath": sns.light_palette(colors[1], sources_number + up).as_hex()[up:],
        "2dPath": sns.light_palette(colors[8], sources_number + up).as_hex()[up:],
        "metConnect": sns.light_palette(colors[0], sources_number + up).as_hex()[up:],
        "connect": sns.light_palette(colors[9], sources_number + up).as_hex()[up:],
        "interest": sns.light_palette(colors[2], sources_number + up).as_hex()[up:],
        "notFound": sns.light_palette(colors[7], sources_number + up).as_hex()[up:],
    }
    return out_colors


def defineNodeColor(
    colordata: dict, pallitra: str, object: NewObject, n_letter: int
) -> [str, str]:
    col = colordata.get(pallitra)[object.in_models["models_amount"] - 1]
    name = ""
    for tmp in object.in_models["models_list"]:
        name = name + tmp[:n_letter]
    return [col, object.id + "\n" + name]


def defineEdgeColor(
    colordata: dict,
    pallitra_forward: str,
    pallitra_reversed: str,
    mr: NewObject,
    colname_m: str,
    r: NewObject,
    colname_r: str,
    connections: str,
) -> [[str, str, str]]:
    n = 0
    if r.lower_bound.get("assembly")[0] >= 0:
        for tmp in r.in_models["models_list"]:
            if mr in getattr(r, connections).get(tmp):
                n = n + 1
        if connections == "reactants":
            col = colordata.get(pallitra_forward)[n - 1]
            output = [[colname_m, colname_r, col]]
        elif connections == "products":
            col = colordata.get(pallitra_reversed)[n - 1]
            output = [[colname_r, colname_m, col]]
    elif r.upper_bound.get("assembly")[0] <= 0:
        for tmp in r.in_models["models_list"]:
            if mr in getattr(r, connections).get(tmp):
                n = n + 1
        if connections == "reactants":
            col = colordata.get(pallitra_reversed)[n - 1]
            output = [[colname_r, colname_m, col]]
        elif connections == "products":
            col = colordata.get(pallitra_forward)[n - 1]
            output = [[colname_m, colname_r, col]]
    else:
        f = 0
        b = 0
        for tmp in r.in_models["models_list"]:
            if mr in getattr(r, connections).get(tmp) and r.upper_bound.get(tmp)[0] > 0:
                f = f + 1
            if mr in getattr(r, connections).get(tmp) and r.lower_bound.get(tmp)[0] < 0:
                b = b + 1
        if connections == "reactants":
            col_f = colordata.get(pallitra_forward)[f - 1]
            col_b = colordata.get(pallitra_forward)[b - 1]
            output = [
                [colname_m, colname_r, col_f],
                [colname_r, colname_m, col_b],
            ]
        elif connections == "products":
            col_f = colordata.get(pallitra_reversed)[f - 1]
            col_b = colordata.get(pallitra_reversed)[b - 1]
            output = [
                [colname_r, colname_m, col_f],
                [colname_m, colname_r, col_b],
            ]
    return output


def draw_one_known_pathway(
    supermodel: SuperModel,
    pathway: dict,
    output_name: str,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=None,
    hei=None,
):
    if not output_name.endswith(".html"):
        raise ValueError(
            "Output file for the plot is of wrong format. Please use html file name."
        )
    if met_not_int is None:
        met_not_int = {
            "h": 0,
            "h2o": 0,
            "h2": 0,
            "oh1": 0,
            "o2": 0,
            "co2": 0,
            "coa": 0,
            "ppi": 0,
            "pi": 0,
            "amp": 0,
            "adp": 0,
            "atp": 0,
            "cmp": 0,
            "cdp": 0,
            "ctp": 0,
            "gmp": 0,
            "gdp": 0,
            "gtp": 0,
            "ump": 0,
            "udp": 0,
            "utp": 0,
            "nad": 0,
            "nadh": 0,
            "nadp": 0,
            "nadph": 0,
            "dadp": 0,
            "damp": 0,
            "nh3": 0,
            "nh4": 0,
            "fadh2": 0,
            "fad": 0,
            "ac": 0,
            "accoa": 0,
            "h2s": 0,
            "HC00250": 0,
        }
    color_brewer = getColorPalette(len(supermodel.sources))
    if n_letter is None:
        n_letter = supermodel.get_short_name_len()
    g = nx.DiGraph()
    path_met = []
    for value in pathway.values():
        for v in value:
            path_met.append(v[0])
            path_met.append(v[1])
    path_met = list(set(path_met))
    for r_id in pathway.keys():
        if r_id in supermodel.reactions.assembly.keys():
            r = supermodel.reactions.assembly.get(r_id)
            colname_r = defineNodeColor(color_brewer, "1stPath", r, n_letter,)
            g.add_node(colname_r[1], shape="box", color=colname_r[0])
            for rea in r.reactants.get("assembly"):
                tmp_rea = re.sub("_([cep])$", "", rea.id)
                if tmp_rea not in met_not_int.keys():
                    if rea.id in path_met:
                        rea_pal = "metabolites"
                    else:
                        rea_pal = "metConnect"
                    colname_rea = defineNodeColor(color_brewer, rea_pal, rea, n_letter,)
                    g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                    rea_edge = defineEdgeColor(
                        color_brewer,
                        "1stPath",
                        "metabolites",
                        rea,
                        colname_rea[1],
                        r,
                        colname_r[1],
                        "reactants",
                    )
                else:
                    g.add_node(
                        f"{tmp_rea} {met_not_int[tmp_rea]}",
                        label=tmp_rea,
                        shape="o",
                        color=color_brewer["metConnect"][
                            math.floor(len(supermodel.sources) / 2)
                        ],
                    )
                    rea_edge = defineEdgeColor(
                        color_brewer,
                        "1stPath",
                        "metabolites",
                        rea,
                        f"{tmp_rea} {met_not_int[tmp_rea]}",
                        r,
                        colname_r[1],
                        "reactants",
                    )
                    met_not_int[tmp_rea] = met_not_int[tmp_rea] + 1
                for e in rea_edge:
                    g.add_edge(
                        e[0], e[1], color=e[2], font_color="black",
                    )
            for pro in r.products.get("assembly"):
                tmp_pro = re.sub("_([cep])$", "", pro.id)
                if tmp_pro not in met_not_int.keys():
                    if pro.id in path_met:
                        pro_pal = "metabolites"
                    else:
                        pro_pal = "metConnect"
                    colname_pro = defineNodeColor(color_brewer, pro_pal, pro, n_letter,)

                    g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                    pro_edge = defineEdgeColor(
                        color_brewer,
                        "1stPath",
                        "metabolites",
                        pro,
                        colname_pro[1],
                        r,
                        colname_r[1],
                        "products",
                    )
                else:
                    g.add_node(
                        f"{tmp_pro} {met_not_int[tmp_pro]}",
                        label=tmp_pro,
                        shape="o",
                        color=color_brewer["metConnect"][
                            math.floor(len(supermodel.sources) / 2)
                        ],
                    )
                    pro_edge = defineEdgeColor(
                        color_brewer,
                        "1stPath",
                        "metabolites",
                        pro,
                        f"{tmp_pro} {met_not_int[tmp_pro]}",
                        r,
                        colname_r[1],
                        "products",
                    )
                    met_not_int[tmp_pro] = met_not_int[tmp_pro] + 1
                for e in pro_edge:
                    g.add_edge(
                        e[0], e[1], color=e[2], font_color="black",
                    )
        else:
            not_f_col = color_brewer["notFound"][
                math.floor(len(supermodel.sources) / 2)
            ]
            g.add_node(
                r_id, shape="box", color=not_f_col,
            )
            for pair in pathway.get(r_id):
                if pair[0] in supermodel.metabolites.assembly.keys():
                    colname_rea = defineNodeColor(
                        color_brewer,
                        "metabolites",
                        supermodel.metabolites.assembly[pair[0]],
                        n_letter,
                    )
                else:
                    colname_rea = [not_f_col, pair[0]]
                if pair[1] in supermodel.metabolites.assembly.keys():
                    colname_pro = defineNodeColor(
                        color_brewer,
                        "metabolites",
                        supermodel.metabolites.assembly[pair[1]],
                        n_letter,
                    )
                else:
                    colname_pro = [not_f_col, pair[1]]
                g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                g.add_edge(
                    colname_rea[1], r_id, color=not_f_col, font_color="black",
                )
                g.add_edge(
                    r_id, colname_pro[1], color=not_f_col, font_color="black",
                )
    if wid is None or hei is None:
        wid, hei = pyautogui.size()
    pyvis_graph = Network(
        width="{}px".format(wid),
        height="{}px".format(hei),
        directed=directed,
        notebook=False,
    )
    pyvis_graph.from_nx(g)
    pyvis_graph.show(output_name, notebook=False)
    return g


def draw_tca(
    supermodel: SuperModel,
    output_name: str,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=None,
    hei=None,
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
    g = draw_one_known_pathway(
        supermodel, tca, output_name, directed, met_not_int, n_letter, wid, hei
    )
    return g


def draw_glycolysis(
    supermodel: SuperModel,
    output_name: str,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=None,
    hei=None,
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
    g = draw_one_known_pathway(
        supermodel, glycolysis, output_name, directed, met_not_int, n_letter, wid, hei
    )
    return g


def draw_pentose_phosphate(
    supermodel: SuperModel,
    output_name: str,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=None,
    hei=None,
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
    g = draw_one_known_pathway(
        supermodel,
        pentose_phosphate_pathway,
        output_name,
        directed,
        met_not_int,
        n_letter,
        wid,
        hei,
    )
    return g


def drawTwoPathways(
    supermodel,
    pathway0,
    pathway1,
    met_not_int,
    colorBrewer,
    name,
    directed=False,
    Nletter=1,
):
    g = nx.DiGraph()
    for r_id in pathway1["reactions"]:
        r = supermodel.reactions.converted.get(r_id)
        for rea in r.reactants.get("union1"):
            tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
            for pro in r.products.get("union1"):
                tmp_pro = (
                    pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                )
                if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                    colname_r = defineNodeColor(
                        colorBrewer, "greens", r_id, supermodel.reactions, Nletter
                    )
                    if rea.id in pathway1["metabolites"]:
                        colname_rea = defineNodeColor(
                            colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter
                        )
                        cn_rea = defineEdgeColor(
                            colorBrewer, "reds", rea, r.reactants, Nletter
                        )
                    else:
                        colname_rea = defineNodeColor(
                            colorBrewer,
                            "blues",
                            rea.id,
                            supermodel.metabolites,
                            Nletter,
                        )
                        cn_rea = defineEdgeColor(
                            colorBrewer, "blues", rea, r.reactants, Nletter
                        )

                    if pro.id in pathway1["metabolites"]:
                        colname_pro = defineNodeColor(
                            colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter
                        )
                        cn_pro = defineEdgeColor(
                            colorBrewer, "reds", pro, r.products, Nletter
                        )
                    else:
                        colname_pro = defineNodeColor(
                            colorBrewer,
                            "blues",
                            pro.id,
                            supermodel.metabolites,
                            Nletter,
                        )
                        cn_pro = defineEdgeColor(
                            colorBrewer, "blues", pro, r.products, Nletter
                        )
                    g.add_node(colname_r[1], shape="box", color=colname_r[0])
                    g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                    g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                    g.add_edge(
                        colname_rea[1],
                        colname_r[1],
                        color=cn_rea[0],
                        font_color="black",
                    )
                    g.add_edge(
                        colname_r[1],
                        colname_pro[1],
                        color=cn_pro[0],
                        font_color="black",
                    )

    for r_id in pathway0["reactions"]:
        r = supermodel.reactions.converted.get(r_id)
        for rea in r.reactants.get("union1"):
            tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
            for pro in r.products.get("union1"):
                tmp_pro = (
                    pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                )
                if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                    colname_r = defineNodeColor(
                        colorBrewer, "purples", r_id, supermodel.reactions, Nletter
                    )
                    if rea.id in pathway0["metabolites"]:
                        colname_rea = defineNodeColor(
                            colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter
                        )
                        cn_rea = defineEdgeColor(
                            colorBrewer, "reds", rea, r.reactants, Nletter
                        )
                    else:
                        colname_rea = defineNodeColor(
                            colorBrewer,
                            "blues",
                            rea.id,
                            supermodel.metabolites,
                            Nletter,
                        )
                        cn_rea = defineEdgeColor(
                            colorBrewer, "blues", rea, r.reactants, Nletter
                        )

                    if pro.id in pathway0["metabolites"]:
                        colname_pro = defineNodeColor(
                            colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter
                        )
                        cn_pro = defineEdgeColor(
                            colorBrewer, "reds", pro, r.products, Nletter
                        )
                    else:
                        colname_pro = defineNodeColor(
                            colorBrewer,
                            "blues",
                            pro.id,
                            supermodel.metabolites,
                            Nletter,
                        )
                        cn_pro = defineEdgeColor(
                            colorBrewer, "blues", pro, r.products, Nletter
                        )
                    g.add_node(colname_r[1], shape="box", color=colname_r[0])
                    g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                    g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                    g.add_edge(
                        colname_rea[1],
                        colname_r[1],
                        color=cn_rea[0],
                        font_color="black",
                    )
                    g.add_edge(
                        colname_r[1],
                        colname_pro[1],
                        color=cn_pro[0],
                        font_color="black",
                    )

    wid, hei = pyautogui.size()
    pyvis_graph = Network(
        width="{}px".format(wid),
        height="{}px".format(hei),
        directed=directed,
        notebook=False,
    )
    pyvis_graph.from_nx(g)
    pyvis_graph.show("../Output/" + name + ".html")
    return g


def drawCore(
    supermodel,
    met_not_int,
    colorBrewer,
    name,
    union=False,
    directed=False,
    wid=2000,
    hei=1000,
    Nletter=1,
    union_size=1,
    core_size=None,
):
    g = nx.DiGraph()
    if union:
        attribute = "converted"
        attribute_connect = "union" + str(union_size)
    else:
        if core_size:
            attribute = "core" + str(core_size)
        else:
            attribute = "core" + str(len(supermodel.sources))
        attribute_connect = attribute
    for r in getattr(supermodel.reactions, attribute).values():
        if r.id != "Biomass":
            for rea in r.reactants.get(attribute_connect):
                tmp_rea = (
                    rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                )
                for pro in r.products.get(attribute_connect):
                    tmp_pro = (
                        pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                    )
                    if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                        colname_r = defineNodeColor(
                            colorBrewer, "purples", r.id, supermodel.reactions, Nletter
                        )
                        colname_rea = defineNodeColor(
                            colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter
                        )
                        colname_pro = defineNodeColor(
                            colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter
                        )
                        cn_rea = defineEdgeColor(
                            colorBrewer, "purples", rea, r.reactants, Nletter
                        )
                        cn_pro = defineEdgeColor(
                            colorBrewer, "purples", pro, r.products, Nletter
                        )
                        g.add_node(colname_r[1], shape="box", color=colname_r[0])
                        g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                        g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                        g.add_edge(
                            colname_rea[1],
                            colname_r[1],
                            color=cn_rea[0],
                            font_color="black",
                        )
                        g.add_edge(
                            colname_r[1],
                            colname_pro[1],
                            color=cn_pro[0],
                            font_color="black",
                        )

    pyvis_graph = Network(
        width="{}px".format(wid),
        height="{}px".format(hei),
        directed=directed,
        notebook=False,
    )
    pyvis_graph.from_nx(g)
    pyvis_graph.show("../Output/" + name + ".html")
    return g


def drawBiomass(
    supermodel,
    name,
    only_difference=False,
    not_converted=False,
    colorBrewer=None,
    directed=True,
    wid=2000,
    hei=1000,
    Nletter=1,
    union_size=1,
    core_size=None,
):
    if core_size:
        core = "core" + str(core_size)
    else:
        core = "core" + str(len(supermodel.sources))
    if only_difference:
        g = nx.DiGraph()
        biomass_r = supermodel.reactions.converted.get("Biomass")
        for rea in biomass_r.reactants.get("union" + str(union_size)):
            for pro in biomass_r.products.get("union" + str(union_size)):
                colname_r = defineNodeColor(
                    colorBrewer, "purples", biomass_r.id, supermodel.reactions, Nletter
                )
                g.add_node(colname_r[1], shape="box", color=colname_r[0])
                if rea not in biomass_r.reactants.get(core):
                    colname_rea = defineNodeColor(
                        colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter
                    )
                    cn_rea = defineEdgeColor(
                        colorBrewer, "purples", rea, biomass_r.reactants, Nletter
                    )
                    g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                    g.add_edge(
                        colname_rea[1],
                        colname_r[1],
                        label=cn_rea[1],
                        color=cn_rea[0],
                        font_color="black",
                    )
                if pro not in biomass_r.products.get(core):
                    colname_pro = defineNodeColor(
                        colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter
                    )
                    cn_pro = defineEdgeColor(
                        colorBrewer, "purples", pro, biomass_r.products, Nletter
                    )
                    g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                    g.add_edge(
                        colname_r[1],
                        colname_pro[1],
                        label=cn_pro[1],
                        color=cn_pro[0],
                        font_color="black",
                    )
    else:
        if not_converted:
            g = nx.DiGraph()
            biomass_r = supermodel.reactions.notconverted.get("Biomass")
            reactants = {}
            products = {}
            for typ in supermodel.sources:
                for rea in biomass_r.reactants.get(typ):
                    if rea.id not in reactants.keys():
                        if rea.id in rea.annotation.get(typ):
                            reactants.update({rea.id: f"{rea.id}\n{typ[0]}"})
                        else:
                            reactants.update(
                                {
                                    rea.id: f"{rea.id}\n{rea.annotation.get(typ)}\n{typ[0]}"
                                }
                            )
                    else:
                        if rea.id in rea.annotation.get(typ):
                            reactants.update(
                                {rea.id: f"{reactants.get(rea.id)}\n{rea.id}\n{typ[0]}"}
                            )
                        else:
                            reactants.update(
                                {
                                    rea.id: f"{reactants.get(rea.id)}\n{rea.id}\n{rea.annotation.get(typ)}\n{typ[0]}"
                                }
                            )
                for pro in biomass_r.products.get(typ):
                    if pro.id not in products.keys():
                        if pro.id in pro.annotation.get(typ):
                            products.update({pro.id: f"{pro.id}\n{typ[0]}"})
                        else:
                            reactants.update(
                                {
                                    pro.id: f"{pro.id}\n{pro.annotation.get(typ)}\n{typ[0]}"
                                }
                            )
                    else:
                        if pro.id in pro.annotation.get(typ):
                            products.update(
                                {pro.id: f"{products.get(pro.id)}\n{pro.id}\n{typ[0]}"}
                            )
                        else:
                            products.update(
                                {
                                    pro.id: f"{products.get(pro.id)}\n{pro.id}\n{pro.annotation.get(typ)}\n{typ[0]}"
                                }
                            )
            g.add_node(biomass_r.id, shape="box")
            for r in reactants.values():
                g.add_node(r, shape="o")
                g.add_edge(r, biomass_r.id, font_color="black")
            for p in products.values():
                g.add_node(p, shape="o")
                g.add_edge(biomass_r.id, p, font_color="black")
        else:
            g = nx.DiGraph()
            biomass_r = supermodel.reactions.converted.get("Biomass")
            for rea in biomass_r.reactants.get("union" + str(union_size)):
                for pro in biomass_r.products.get("union" + str(union_size)):
                    colname_r = defineNodeColor(
                        colorBrewer,
                        "purples",
                        biomass_r.id,
                        supermodel.reactions,
                        Nletter,
                    )
                    colname_rea = defineNodeColor(
                        colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter
                    )
                    colname_pro = defineNodeColor(
                        colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter
                    )
                    cn_rea = defineEdgeColor(
                        colorBrewer, "purples", rea, biomass_r.reactants, Nletter
                    )
                    cn_pro = defineEdgeColor(
                        colorBrewer, "purples", pro, biomass_r.products, Nletter
                    )
                    g.add_node(colname_r[1], shape="box", color=colname_r[0])
                    g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                    g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                    g.add_edge(
                        colname_rea[1],
                        colname_r[1],
                        color=cn_rea[0],
                        font_color="black",
                    )
                    g.add_edge(
                        colname_r[1],
                        colname_pro[1],
                        color=cn_pro[0],
                        font_color="black",
                    )

    pyvis_graph = Network(
        width="{}px".format(wid),
        height="{}px".format(hei),
        directed=directed,
        notebook=False,
    )
    pyvis_graph.from_nx(g)
    pyvis_graph.show("../Output/" + name + ".html")
    return g


def drawMetSynthesis(
    supermodel,
    met_not_int: list,
    colorBrewer: dict,
    met_of_int: str,
    paths,
    sources: list,
    level: str,
    name_dir: str,
    directed=True,
    wid=2000,
    hei=1000,
    Nletter=1,
):
    g = nx.DiGraph()
    r_done = []
    if type(paths) == list:
        for path in paths:
            for r_id in path:
                if r_id not in r_done:
                    r = supermodel.reactions.converted.get(r_id)
                    for rea in r.reactants.get("core1"):
                        tmp_rea = (
                            rea.id.removesuffix("_c")
                            .removesuffix("_e")
                            .removesuffix("_p")
                        )
                        for pro in r.products.get("core1"):
                            tmp_pro = (
                                pro.id.removesuffix("_c")
                                .removesuffix("_e")
                                .removesuffix("_p")
                            )
                            if (tmp_rea not in met_not_int) & (
                                tmp_pro not in met_not_int
                            ):
                                colname_r = defineNodeColor(
                                    colorBrewer, "1stPath", r, Nletter
                                )
                                cn_rea = defineEdgeColor(
                                    colorBrewer,
                                    "metabolites",
                                    rea,
                                    r,
                                    r.reactants,
                                    Nletter,
                                )
                                cn_pro = defineEdgeColor(
                                    colorBrewer,
                                    "metabolites",
                                    pro,
                                    r,
                                    r.products,
                                    Nletter,
                                )
                                cn_rea_2 = defineEdgeColor(
                                    colorBrewer, "1stPath", rea, r, r.reactants, Nletter
                                )
                                cn_pro_2 = defineEdgeColor(
                                    colorBrewer, "1stPath", pro, r, r.products, Nletter
                                )
                                if (rea.id == met_of_int) or (rea.id in sources):
                                    colname_rea = defineNodeColor(
                                        colorBrewer, "interest", rea, Nletter
                                    )
                                else:
                                    colname_rea = defineNodeColor(
                                        colorBrewer, "metabolites", rea, Nletter
                                    )
                                if (pro.id == met_of_int) or (pro.id in sources):
                                    colname_pro = defineNodeColor(
                                        colorBrewer, "interest", pro, Nletter
                                    )
                                else:
                                    colname_pro = defineNodeColor(
                                        colorBrewer, "metabolites", pro, Nletter
                                    )
                                g.add_node(
                                    colname_r[1], shape="box", color=colname_r[0]
                                )
                                g.add_node(
                                    colname_rea[1], shape="o", color=colname_rea[0]
                                )
                                g.add_node(
                                    colname_pro[1], shape="o", color=colname_pro[0]
                                )
                                if r.lower_bound.get(level)[0] >= 0:
                                    g.add_edge(
                                        colname_rea[1],
                                        colname_r[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_pro[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                if r.upper_bound.get(level)[0] <= 0:
                                    g.add_edge(
                                        colname_pro[1],
                                        colname_r[1],
                                        color=cn_pro[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_rea[1],
                                        color=cn_rea_2[0],
                                        font_color="black",
                                    )
                                else:
                                    g.add_edge(
                                        colname_rea[1],
                                        colname_r[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_pro[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_pro[1],
                                        colname_r[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_rea[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                    r_done.append(r_id)
        pyvis_graph = Network(
            width="{}px".format(wid),
            height="{}px".format(hei),
            directed=directed,
            notebook=False,
        )
        pyvis_graph.from_nx(g)
        pyvis_graph.show("../Output/" + name_dir + "/" + met_of_int + ".html")
    else:
        print(f"No paths for {met_of_int}, because {paths}")
    return g


def drawMetSynthesis1model(
    g,
    col_r,
    supermodel,
    met_not_int: list,
    colorBrewer: dict,
    met_of_int: str,
    paths,
    sources: list,
    level: str,
    directed=True,
    wid=2000,
    hei=1000,
    Nletter=1,
):
    r_done = []
    if type(paths) == list:
        for path in paths:
            for r_id in path:
                if r_id not in r_done:
                    r = supermodel.reactions.converted.get(r_id)
                    for rea in r.reactants.get("core1"):
                        tmp_rea = (
                            rea.id.removesuffix("_c")
                            .removesuffix("_e")
                            .removesuffix("_p")
                        )
                        for pro in r.products.get("core1"):
                            tmp_pro = (
                                pro.id.removesuffix("_c")
                                .removesuffix("_e")
                                .removesuffix("_p")
                            )
                            if (tmp_rea not in met_not_int) & (
                                tmp_pro not in met_not_int
                            ):
                                colname_r = defineNodeColor(
                                    colorBrewer, col_r, r, Nletter
                                )
                                cn_rea = defineEdgeColor(
                                    colorBrewer,
                                    "metabolites",
                                    rea,
                                    r,
                                    r.reactants,
                                    Nletter,
                                )
                                cn_pro = defineEdgeColor(
                                    colorBrewer,
                                    "metabolites",
                                    pro,
                                    r,
                                    r.products,
                                    Nletter,
                                )
                                cn_rea_2 = defineEdgeColor(
                                    colorBrewer, col_r, rea, r, r.reactants, Nletter
                                )
                                cn_pro_2 = defineEdgeColor(
                                    colorBrewer, col_r, pro, r, r.products, Nletter
                                )
                                if (rea.id == met_of_int) or (rea.id in sources):
                                    colname_rea = defineNodeColor(
                                        colorBrewer, "interest", rea, Nletter
                                    )
                                else:
                                    colname_rea = defineNodeColor(
                                        colorBrewer, "metabolites", rea, Nletter
                                    )
                                if (pro.id == met_of_int) or (pro.id in sources):
                                    colname_pro = defineNodeColor(
                                        colorBrewer, "interest", pro, Nletter
                                    )
                                else:
                                    colname_pro = defineNodeColor(
                                        colorBrewer, "metabolites", pro, Nletter
                                    )
                                g.add_node(
                                    colname_r[1], shape="box", color=colname_r[0]
                                )
                                g.add_node(
                                    colname_rea[1], shape="o", color=colname_rea[0]
                                )
                                g.add_node(
                                    colname_pro[1], shape="o", color=colname_pro[0]
                                )
                                if r.lower_bound.get(level)[0] >= 0:
                                    g.add_edge(
                                        colname_rea[1],
                                        colname_r[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_pro[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                if r.upper_bound.get(level)[0] <= 0:
                                    g.add_edge(
                                        colname_pro[1],
                                        colname_r[1],
                                        color=cn_pro[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_rea[1],
                                        color=cn_rea_2[0],
                                        font_color="black",
                                    )
                                else:
                                    g.add_edge(
                                        colname_rea[1],
                                        colname_r[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_pro[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_pro[1],
                                        colname_r[1],
                                        color=cn_pro_2[0],
                                        font_color="black",
                                    )
                                    g.add_edge(
                                        colname_r[1],
                                        colname_rea[1],
                                        color=cn_rea[0],
                                        font_color="black",
                                    )
                    r_done.append(r_id)
    return g
