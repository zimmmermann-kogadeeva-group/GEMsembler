import operator
import networkx as nx
import pyautogui
from pyvis.network import Network
import math
from creation import SetofNewReactions, SetofNewMetabolites, NewObject
import seaborn as sns
from general import findKeysByValue


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
        "notFound": sns.light_palette(colors[7], sources_number + up).as_hex()[up:]
    }
    return out_colors


def defineNodeColor0(colordata: dict, pallitra: str, id_mr: str, objects: SetofNewReactions or SetofNewMetabolites,
                    Nletter=1) -> [str, str]:
    for y in dir(objects):
        if y.startswith("Yes_") or y.startswith("core"):
            if id_mr in getattr(objects, y).keys():
                name = y.split("_No_")[0].replace("Yes_", "")
                if y.startswith("core"):
                    n = int(y.replace("core", ""))
                    col = colordata.get(pallitra)[n - 1]
                else:
                    col = colordata.get(pallitra)[int(len(name) / Nletter) - 1]
                return ([col, id_mr + "\n" + name])


def defineEdgeColor0(colordata: dict, pallitra: str, mr: NewObject, connections: dict, Nletter=1) -> [str, str]:
    for y in connections.keys():
        if y.startswith("Yes_") or y.startswith("core"):
            if mr in connections.get(y):
                name = y.split("_No_")[0].replace("Yes_", "")
                if y.startswith("core"):
                    n = int(y.replace("core", ""))
                    col = colordata.get(pallitra)[n - 1]
                else:
                    col = colordata.get(pallitra)[int(len(name) / Nletter) - 1]
                return ([col, name])


def defineNodeColor(colordata: dict, pallitra: str, object: NewObject, Nletter=1) -> [str, str]:
    tmp_models = findKeysByValue(object.sources, 1, operator.ge)
    n = len(tmp_models)
    col = colordata.get(pallitra)[n - 1]
    name = ""
    for tmp in tmp_models:
        name = name + tmp[:Nletter]
    return ([col, object.id + "\n" + name])


def defineEdgeColor(colordata: dict, pallitra: str, mr: NewObject, r:NewObject, connections: dict, Nletter=1) -> [str, str]:
    tmp_models = findKeysByValue(r.sources, 1, operator.ge)
    n = 0
    name = ""
    for tmp in tmp_models:
        if mr in connections.get(tmp):
            n = n + 1
            name = tmp[:Nletter]
    col = colordata.get(pallitra)[n - 1]
    return ([col, name])


def drawOnePathway(supermodel, pathway, met_not_int, colorBrewer, name, aminoacids=None, directed=False,
                   surrounding=False, Nletter=1):
    g = nx.DiGraph()
    for r_id in pathway["reactions"]:
        r = supermodel.reactions.converted.get(r_id)
        for rea in r.reactants.get("union1"):
            tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
            for pro in r.products.get("union1"):
                tmp_pro = pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                    colname_r = defineNodeColor(colorBrewer, "1stPath", r_id, supermodel.reactions, Nletter)
                    if rea.id in pathway["metabolites"]:
                        colname_rea = defineNodeColor(colorBrewer, "metabolies", rea.id, supermodel.metabolites,
                                                      Nletter)
                        cn_rea = defineEdgeColor(colorBrewer, "1stPath", rea, r.reactants, Nletter)
                    else:
                        colname_rea = defineNodeColor(colorBrewer, "metConnect", rea.id, supermodel.metabolites,
                                                      Nletter)
                        cn_rea = defineEdgeColor(colorBrewer, "connect", rea, r.reactants, Nletter)

                    if pro.id in pathway["metabolites"]:
                        colname_pro = defineNodeColor(colorBrewer, "metabolies", pro.id, supermodel.metabolites,
                                                      Nletter)
                        cn_pro = defineEdgeColor(colorBrewer, "1stPath", pro, r.products, Nletter)
                    else:
                        colname_pro = defineNodeColor(colorBrewer, "metConnect", pro.id, supermodel.metabolites,
                                                      Nletter)
                        cn_pro = defineEdgeColor(colorBrewer, "connect", pro, r.products, Nletter)
                    g.add_node(colname_r[1], shape="box", color=colname_r[0])
                    g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                    g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                    g.add_edge(colname_rea[1], colname_r[1], color=cn_rea[0], font_color="black")
                    g.add_edge(colname_r[1], colname_pro[1], color=cn_pro[0], font_color="black")

    if surrounding:
        for m_id in pathway["metabolites"]:
            m = supermodel.metabolites.converted.get(m_id)
            for ra in m.reactions.get("union1"):
                if ra.id not in pathway["reactions"]:
                    for rea1 in ra.reactants.get("union1"):
                        tmp_rea1 = rea1.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                        for pro1 in ra.products.get("union1"):
                            tmp_pro1 = pro1.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                            if (tmp_rea1 not in met_not_int) & (tmp_pro1 not in met_not_int):
                                if (tmp_rea1 in aminoacids) or (tmp_pro1 in aminoacids):
                                    colname_rn = defineNodeColor(colorBrewer, "connect", ra.id, supermodel.reactions,
                                                                 Nletter)
                                    colname_rean = defineNodeColor(colorBrewer, "metConnect", rea1.id,
                                                                   supermodel.metabolites, Nletter)
                                    colname_pron = defineNodeColor(colorBrewer, "metConnect", pro1.id,
                                                                   supermodel.metabolites, Nletter)
                                    cn_rean = defineEdgeColor(colorBrewer, "connect", rea1, ra.reactants, Nletter)
                                    cn_pron = defineEdgeColor(colorBrewer, "connect", pro1, ra.products, Nletter)
                                    g.add_node(colname_rn[1], shape="box", color=colname_rn[0])
                                    if rea1.id not in pathway["metabolites"]:
                                        g.add_node(colname_rean[1], shape="o", color=colname_rean[0])
                                    if pro1.id not in pathway["metabolites"]:
                                        g.add_node(colname_pron[1], shape="o", color=colname_pron[0])
                                    g.add_edge(colname_rean[1], colname_rn[1], color=cn_rean[0], font_color="black")
                                    g.add_edge(colname_rn[1], colname_pron[1], color=cn_pron[0], font_color="black")

    wid, hei = pyautogui.size()
    pyvis_graph = Network(width='{}px'.format(wid), height='{}px'.format(hei), directed=directed, notebook=False)
    pyvis_graph.from_nx(g)
    pyvis_graph.show("../Output/" + name + ".html")
    return g


def drawTwoPathways(supermodel, pathway0, pathway1, met_not_int, colorBrewer, name, directed=False, Nletter=1):
    g = nx.DiGraph()
    for r_id in pathway1["reactions"]:
        r = supermodel.reactions.converted.get(r_id)
        for rea in r.reactants.get("union1"):
            tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
            for pro in r.products.get("union1"):
                tmp_pro = pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                    colname_r = defineNodeColor(colorBrewer, "greens", r_id, supermodel.reactions, Nletter)
                    if rea.id in pathway1["metabolites"]:
                        colname_rea = defineNodeColor(colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter)
                        cn_rea = defineEdgeColor(colorBrewer, "reds", rea, r.reactants, Nletter)
                    else:
                        colname_rea = defineNodeColor(colorBrewer, "blues", rea.id, supermodel.metabolites, Nletter)
                        cn_rea = defineEdgeColor(colorBrewer, "blues", rea, r.reactants, Nletter)

                    if pro.id in pathway1["metabolites"]:
                        colname_pro = defineNodeColor(colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter)
                        cn_pro = defineEdgeColor(colorBrewer, "reds", pro, r.products, Nletter)
                    else:
                        colname_pro = defineNodeColor(colorBrewer, "blues", pro.id, supermodel.metabolites, Nletter)
                        cn_pro = defineEdgeColor(colorBrewer, "blues", pro, r.products, Nletter)
                    g.add_node(colname_r[1], shape="box", color=colname_r[0])
                    g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                    g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                    g.add_edge(colname_rea[1], colname_r[1], color=cn_rea[0], font_color="black")
                    g.add_edge(colname_r[1], colname_pro[1], color=cn_pro[0], font_color="black")

    for r_id in pathway0["reactions"]:
        r = supermodel.reactions.converted.get(r_id)
        for rea in r.reactants.get("union1"):
            tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
            for pro in r.products.get("union1"):
                tmp_pro = pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                    colname_r = defineNodeColor(colorBrewer, "purples", r_id, supermodel.reactions, Nletter)
                    if rea.id in pathway0["metabolites"]:
                        colname_rea = defineNodeColor(colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter)
                        cn_rea = defineEdgeColor(colorBrewer, "reds", rea, r.reactants, Nletter)
                    else:
                        colname_rea = defineNodeColor(colorBrewer, "blues", rea.id, supermodel.metabolites, Nletter)
                        cn_rea = defineEdgeColor(colorBrewer, "blues", rea, r.reactants, Nletter)

                    if pro.id in pathway0["metabolites"]:
                        colname_pro = defineNodeColor(colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter)
                        cn_pro = defineEdgeColor(colorBrewer, "reds", pro, r.products, Nletter)
                    else:
                        colname_pro = defineNodeColor(colorBrewer, "blues", pro.id, supermodel.metabolites, Nletter)
                        cn_pro = defineEdgeColor(colorBrewer, "blues", pro, r.products, Nletter)
                    g.add_node(colname_r[1], shape="box", color=colname_r[0])
                    g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                    g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                    g.add_edge(colname_rea[1], colname_r[1], color=cn_rea[0], font_color="black")
                    g.add_edge(colname_r[1], colname_pro[1], color=cn_pro[0], font_color="black")

    wid, hei = pyautogui.size()
    pyvis_graph = Network(width='{}px'.format(wid), height='{}px'.format(hei), directed=directed, notebook=False)
    pyvis_graph.from_nx(g)
    pyvis_graph.show("../Output/" + name + ".html")
    return g


def drawPathways(supermodel, pathway, plot_type, met_not_int, colorBrewer, name, directed=False, surrounding=None,
                 second_pathway=None, Nletter=1, union_size=1):
    path_met = []
    for value in pathway.values():
        for v in value:
            path_met.append(v[0])
            path_met.append(v[1])
    path_met = list(set(path_met))
    g = nx.DiGraph()
    for r_id in pathway.keys():
        if r_id in supermodel.reactions.converted.keys():
            r = supermodel.reactions.converted.get(r_id)
            for rea in r.reactants.get("union" + str(union_size)):
                tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                for pro in r.products.get("union" + str(union_size)):
                    tmp_pro = pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")


def drawTCA(supermodel, pathway, met_not_int, colorBrewer, name, aminoacids=None, surrounding=False, Nletter=1):
    path_met = []
    for value in pathway.values():
        for v in value:
            path_met.append(v[0])
            path_met.append(v[1])
    path_met = list(set(path_met))
    g = nx.DiGraph()
    for r_id in pathway.keys():
        if r_id in supermodel.reactions.converted.keys():
            r = supermodel.reactions.converted.get(r_id)
            for rea in r.reactants.get("union1"):
                tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                for pro in r.products.get("union1"):
                    tmp_pro = pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                    if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                        colname_r = defineNodeColor(colorBrewer, "purples", r_id, supermodel.reactions, Nletter)
                        if rea.id in path_met:
                            colname_rea = defineNodeColor(colorBrewer, "reds", rea.id, supermodel.metabolites,
                                                          Nletter)
                            cn_rea = defineEdgeColor(colorBrewer, "reds", rea, r.reactants, Nletter)
                        else:
                            colname_rea = defineNodeColor(colorBrewer, "blues", rea.id, supermodel.metabolites,
                                                          Nletter)
                            cn_rea = defineEdgeColor(colorBrewer, "blues", rea, r.reactants, Nletter)

                        if pro.id in path_met:
                            colname_pro = defineNodeColor(colorBrewer, "reds", pro.id, supermodel.metabolites,
                                                          Nletter)
                            cn_pro = defineEdgeColor(colorBrewer, "reds", pro, r.products, Nletter)
                        else:
                            colname_pro = defineNodeColor(colorBrewer, "blues", pro.id, supermodel.metabolites,
                                                          Nletter)
                            cn_pro = defineEdgeColor(colorBrewer, "blues", pro, r.products, Nletter)
                        g.add_node(colname_r[1], shape="box", color=colname_r[0])
                        g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                        g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                        g.add_edge(colname_rea[1], colname_r[1], color=cn_rea[0], font_color="black")
                        g.add_edge(colname_r[1], colname_pro[1], color=cn_pro[0], font_color="black")
        else:
            g.add_node(r_id, shape="box", color=colorBrewer.get("greens")[-1])
            for pair in pathway.get(r_id):
                if pair[0] in supermodel.metabolites.converted.keys():
                    colname_rea = defineNodeColor(colorBrewer, "reds", pair[0], supermodel.metabolites, Nletter)
                else:
                    colname_rea = [colorBrewer.get("greens")[-1], pair[0]]
                if pair[1] in supermodel.metabolites.converted.keys():
                    colname_pro = defineNodeColor(colorBrewer, "reds", pair[1], supermodel.metabolites, Nletter)
                else:
                    colname_pro = [colorBrewer.get("greens")[-1], pair[1]]
                g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                g.add_edge(colname_rea[1], r_id, color=colorBrewer.get("greens")[-1], font_color="black")
                g.add_edge(r_id, colname_pro[1], color=colorBrewer.get("greens")[-1], font_color="black")
    wid, hei = pyautogui.size()
    pyvis_graph = Network(width='{}px'.format(wid), height='{}px'.format(hei), directed=False, notebook=False)
    pyvis_graph.from_nx(g)
    pyvis_graph.show("../Output/" + name + ".html")
    return g


def drawCore(supermodel, met_not_int, colorBrewer, name, union=False, directed=False, wid=2000, hei=1000, Nletter=1,
             union_size=1, core_size=None):
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
                tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                for pro in r.products.get(attribute_connect):
                    tmp_pro = pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                    if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                        colname_r = defineNodeColor(colorBrewer, "purples", r.id, supermodel.reactions, Nletter)
                        colname_rea = defineNodeColor(colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter)
                        colname_pro = defineNodeColor(colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter)
                        cn_rea = defineEdgeColor(colorBrewer, "purples", rea, r.reactants, Nletter)
                        cn_pro = defineEdgeColor(colorBrewer, "purples", pro, r.products, Nletter)
                        g.add_node(colname_r[1], shape="box", color=colname_r[0])
                        g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                        g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                        g.add_edge(colname_rea[1], colname_r[1], color=cn_rea[0], font_color="black")
                        g.add_edge(colname_r[1], colname_pro[1], color=cn_pro[0], font_color="black")

    pyvis_graph = Network(width='{}px'.format(wid), height='{}px'.format(hei), directed=directed, notebook=False)
    pyvis_graph.from_nx(g)
    pyvis_graph.show("../Output/" + name + ".html")
    return g


def drawBiomass(supermodel, name, only_difference=False, not_converted=False, colorBrewer=None, directed=True, wid=2000,
                hei=1000, Nletter=1, union_size=1, core_size=None):
    if core_size:
        core = "core" + str(core_size)
    else:
        core = "core" + str(len(supermodel.sources))
    if only_difference:
        g = nx.DiGraph()
        biomass_r = supermodel.reactions.converted.get("Biomass")
        for rea in biomass_r.reactants.get("union" + str(union_size)):
            for pro in biomass_r.products.get("union" + str(union_size)):
                colname_r = defineNodeColor(colorBrewer, "purples", biomass_r.id, supermodel.reactions, Nletter)
                g.add_node(colname_r[1], shape="box", color=colname_r[0])
                if rea not in biomass_r.reactants.get(core):
                    colname_rea = defineNodeColor(colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter)
                    cn_rea = defineEdgeColor(colorBrewer, "purples", rea, biomass_r.reactants, Nletter)
                    g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                    g.add_edge(colname_rea[1], colname_r[1], label=cn_rea[1], color=cn_rea[0], font_color="black")
                if pro not in biomass_r.products.get(core):
                    colname_pro = defineNodeColor(colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter)
                    cn_pro = defineEdgeColor(colorBrewer, "purples", pro, biomass_r.products, Nletter)
                    g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                    g.add_edge(colname_r[1], colname_pro[1], label=cn_pro[1], color=cn_pro[0], font_color="black")
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
                            reactants.update({rea.id: f"{rea.id}\n{rea.annotation.get(typ)}\n{typ[0]}"})
                    else:
                        if rea.id in rea.annotation.get(typ):
                            reactants.update({rea.id: f"{reactants.get(rea.id)}\n{rea.id}\n{typ[0]}"})
                        else:
                            reactants.update(
                                {rea.id: f"{reactants.get(rea.id)}\n{rea.id}\n{rea.annotation.get(typ)}\n{typ[0]}"})
                for pro in biomass_r.products.get(typ):
                    if pro.id not in products.keys():
                        if pro.id in pro.annotation.get(typ):
                            products.update({pro.id: f"{pro.id}\n{typ[0]}"})
                        else:
                            reactants.update({pro.id: f"{pro.id}\n{pro.annotation.get(typ)}\n{typ[0]}"})
                    else:
                        if pro.id in pro.annotation.get(typ):
                            products.update({pro.id: f"{products.get(pro.id)}\n{pro.id}\n{typ[0]}"})
                        else:
                            products.update(
                                {pro.id: f"{products.get(pro.id)}\n{pro.id}\n{pro.annotation.get(typ)}\n{typ[0]}"})
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
                    colname_r = defineNodeColor(colorBrewer, "purples", biomass_r.id, supermodel.reactions, Nletter)
                    colname_rea = defineNodeColor(colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter)
                    colname_pro = defineNodeColor(colorBrewer, "reds", pro.id, supermodel.metabolites, Nletter)
                    cn_rea = defineEdgeColor(colorBrewer, "purples", rea, biomass_r.reactants, Nletter)
                    cn_pro = defineEdgeColor(colorBrewer, "purples", pro, biomass_r.products, Nletter)
                    g.add_node(colname_r[1], shape="box", color=colname_r[0])
                    g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                    g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                    g.add_edge(colname_rea[1], colname_r[1], color=cn_rea[0], font_color="black")
                    g.add_edge(colname_r[1], colname_pro[1], color=cn_pro[0], font_color="black")

    pyvis_graph = Network(width='{}px'.format(wid), height='{}px'.format(hei), directed=directed, notebook=False)
    pyvis_graph.from_nx(g)
    pyvis_graph.show("../Output/" + name + ".html")
    return g


def drawMetSynthesis(supermodel, met_not_int: list, colorBrewer: dict, met_of_int: str, paths, sources: list, level: str,
                     name_dir: str, directed=True, wid=2000, hei=1000, Nletter=1):
    g = nx.DiGraph()
    r_done = []
    if type(paths) == list:
        for path in paths:
            for r_id in path:
                if r_id not in r_done:
                    r = supermodel.reactions.converted.get(r_id)
                    for rea in r.reactants.get("core1"):
                        tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                        for pro in r.products.get("core1"):
                            tmp_pro = pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                            if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                                colname_r = defineNodeColor(colorBrewer, "1stPath", r, Nletter)
                                cn_rea = defineEdgeColor(colorBrewer, "metabolites", rea, r, r.reactants, Nletter)
                                cn_pro = defineEdgeColor(colorBrewer, "metabolites", pro, r, r.products, Nletter)
                                cn_rea_2 = defineEdgeColor(colorBrewer, "1stPath", rea, r, r.reactants, Nletter)
                                cn_pro_2 = defineEdgeColor(colorBrewer, "1stPath", pro, r, r.products, Nletter)
                                if (rea.id == met_of_int) or (rea.id in sources):
                                    colname_rea = defineNodeColor(colorBrewer, "interest", rea, Nletter)
                                else:
                                    colname_rea = defineNodeColor(colorBrewer, "metabolites", rea, Nletter)
                                if (pro.id == met_of_int) or (pro.id in sources):
                                    colname_pro = defineNodeColor(colorBrewer, "interest", pro, Nletter)
                                else:
                                    colname_pro = defineNodeColor(colorBrewer, "metabolites", pro, Nletter)
                                g.add_node(colname_r[1], shape="box", color=colname_r[0])
                                g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                                g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                                if r.lower_bound.get(level)[0] >= 0:
                                    g.add_edge(colname_rea[1], colname_r[1], color=cn_rea[0], font_color="black")
                                    g.add_edge(colname_r[1], colname_pro[1], color=cn_pro_2[0], font_color="black")
                                if r.upper_bound.get(level)[0] <= 0:
                                    g.add_edge(colname_pro[1], colname_r[1], color=cn_pro[0], font_color="black")
                                    g.add_edge(colname_r[1], colname_rea[1], color=cn_rea_2[0], font_color="black")
                                else:
                                    g.add_edge(colname_rea[1], colname_r[1], color=cn_rea[0], font_color="black")
                                    g.add_edge(colname_r[1], colname_pro[1], color=cn_pro_2[0], font_color="black")
                                    g.add_edge(colname_pro[1], colname_r[1], color=cn_pro_2[0], font_color="black")
                                    g.add_edge(colname_r[1], colname_rea[1], color=cn_rea[0], font_color="black")
                    r_done.append(r_id)
        pyvis_graph = Network(width='{}px'.format(wid), height='{}px'.format(hei), directed=directed, notebook=False)
        pyvis_graph.from_nx(g)
        pyvis_graph.show("../Output/" + name_dir + "/" + met_of_int + ".html")
    else:
        print(f"No paths for {met_of_int}, because {paths}")
    return g

def drawMetSynthesis1model(g, col_r, supermodel, met_not_int: list, colorBrewer: dict, met_of_int: str, paths, sources: list, level: str,
                     directed=True, wid=2000, hei=1000, Nletter=1):
    r_done = []
    if type(paths) == list:
        for path in paths:
            for r_id in path:
                if r_id not in r_done:
                    r = supermodel.reactions.converted.get(r_id)
                    for rea in r.reactants.get("core1"):
                        tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                        for pro in r.products.get("core1"):
                            tmp_pro = pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                            if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                                colname_r = defineNodeColor(colorBrewer, col_r, r, Nletter)
                                cn_rea = defineEdgeColor(colorBrewer, "metabolites", rea, r, r.reactants, Nletter)
                                cn_pro = defineEdgeColor(colorBrewer, "metabolites", pro, r, r.products, Nletter)
                                cn_rea_2 = defineEdgeColor(colorBrewer, col_r, rea, r, r.reactants, Nletter)
                                cn_pro_2 = defineEdgeColor(colorBrewer, col_r, pro, r, r.products, Nletter)
                                if (rea.id == met_of_int) or (rea.id in sources):
                                    colname_rea = defineNodeColor(colorBrewer, "interest", rea, Nletter)
                                else:
                                    colname_rea = defineNodeColor(colorBrewer, "metabolites", rea, Nletter)
                                if (pro.id == met_of_int) or (pro.id in sources):
                                    colname_pro = defineNodeColor(colorBrewer, "interest", pro, Nletter)
                                else:
                                    colname_pro = defineNodeColor(colorBrewer, "metabolites", pro, Nletter)
                                g.add_node(colname_r[1], shape="box", color=colname_r[0])
                                g.add_node(colname_rea[1], shape="o", color=colname_rea[0])
                                g.add_node(colname_pro[1], shape="o", color=colname_pro[0])
                                if r.lower_bound.get(level)[0] >= 0:
                                    g.add_edge(colname_rea[1], colname_r[1], color=cn_rea[0], font_color="black")
                                    g.add_edge(colname_r[1], colname_pro[1], color=cn_pro_2[0], font_color="black")
                                if r.upper_bound.get(level)[0] <= 0:
                                    g.add_edge(colname_pro[1], colname_r[1], color=cn_pro[0], font_color="black")
                                    g.add_edge(colname_r[1], colname_rea[1], color=cn_rea_2[0], font_color="black")
                                else:
                                    g.add_edge(colname_rea[1], colname_r[1], color=cn_rea[0], font_color="black")
                                    g.add_edge(colname_r[1], colname_pro[1], color=cn_pro_2[0], font_color="black")
                                    g.add_edge(colname_pro[1], colname_r[1], color=cn_pro_2[0], font_color="black")
                                    g.add_edge(colname_r[1], colname_rea[1], color=cn_rea[0], font_color="black")
                    r_done.append(r_id)
    return g
