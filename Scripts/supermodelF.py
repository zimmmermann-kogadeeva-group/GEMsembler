import itertools
import operator
from collections import Counter
from creation import SuperModel, NewObject, SetofNewReactions, SetofNewMetabolites
from general import findKeysByValue
import networkx as nx
from pyvis.network import Network
import pyautogui


def getCoreConnections(connections: dict, core_size: int, sources: [str]) -> [str]:
    """ Getting connections (reactants/products/ for reaction or reactions for metabolites) that are present in more
     then core_size sources = original models. """
    all_connections = []
    for s in sources:
        connection = connections.get(s)
        for oneconnect in connection:
            all_connections.append(oneconnect)
    counted_connection = Counter(all_connections)
    selected_connection = []
    for key, value in counted_connection.items():
        if value >= core_size:
            selected_connection.append(key)
    return selected_connection


def getCore(supermodel: SuperModel, core_size=None, union_size=1):
    """ Getting supermodel core: intersection of at least core_size amount of sources (by default, intersection of all
     sources). Getting supermodel union of all sources. """
    if not core_size:
        core_size = len(supermodel.sources)
    coreN = "core" + str(core_size)
    unionN = "union" + str(union_size)
    setattr(supermodel.metabolites, coreN, {})
    setattr(supermodel.reactions, coreN, {})
    for met in supermodel.metabolites.converted.values():
        tmp_models = findKeysByValue(met.sources, 1, operator.ge)
        core_r = getCoreConnections(met.reactions, core_size, supermodel.sources)
        union_r = getCoreConnections(met.reactions, union_size, supermodel.sources)
        met.reactions.update({coreN: core_r, unionN: union_r})
        if len(tmp_models) >= core_size:
            getattr(supermodel.metabolites, coreN).update({met.id: met})
    for react in supermodel.reactions.converted.values():
        tmp_models = findKeysByValue(react.sources, 1, operator.ge)
        core_reactants = getCoreConnections(react.reactants, core_size, supermodel.sources)
        core_products = getCoreConnections(react.products, core_size, supermodel.sources)
        u_reactants = getCoreConnections(react.reactants, union_size, supermodel.sources)
        u_products = getCoreConnections(react.products, union_size, supermodel.sources)
        react.reactants.update({coreN: core_reactants, unionN: u_reactants})
        react.products.update({coreN: core_products, unionN: u_products})
        if len(tmp_models) >= core_size:
            getattr(supermodel.reactions, coreN).update({react.id: react})


def getDifConnections(connections: dict, sourceIn: [str], sourceNotIn: [str]) -> [NewObject]:
    """ Getting connections (reactants/products/ for reaction or reactions for metabolites) that are present in 
     "sourceIn" list of sources = original models and not present in "sourceNotIn" list of sources = original models. """
    connectionIn = set(connections.get(sourceIn[0]))
    for sIn in sourceIn:
        connectionIn = connectionIn & set(connections.get(sIn))
    connectionNotIn = []
    for sNotIn in sourceNotIn:
        connectionNotIn = connectionNotIn + connections.get(sNotIn)
    difConnection = set(connectionIn) - set(connectionNotIn)
    return list(difConnection)


def getDifference(supermodel: SuperModel, sourceIn: [str], sourceNotIn: [str], Nletter=1):
    """ Getting metabolites and reactions that are present in "sourceIn" list of sources = original models
    and not present in "sourceNotIn" list of sources = original models. """
    name = "Yes_"
    for sI in sourceIn:
        name = name + sI[:Nletter]
    name = name + "_No_"
    for sNI in sourceNotIn:
        name = name + sNI[:Nletter]
    setattr(supermodel.metabolites, name, {})
    setattr(supermodel.reactions, name, {})
    for met in supermodel.metabolites.converted.values():
        dif_r = getDifConnections(met.reactions, sourceIn, sourceNotIn)
        met.reactions.update({name: dif_r})
        sm_present = findKeysByValue(met.sources, 1, operator.ge)
        sm_absent = findKeysByValue(met.sources, 0, operator.eq)
        if ((set(sourceIn) == set(sm_present)) & (set(sourceNotIn) == set(sm_absent))):
            getattr(supermodel.metabolites, name).update({met.id: met})
    for react in supermodel.reactions.converted.values():
        dif_reactants = getDifConnections(react.reactants, sourceIn, sourceNotIn)
        dif_products = getDifConnections(react.products, sourceIn, sourceNotIn)
        react.reactants.update({name: dif_reactants})
        react.products.update({name: dif_products})
        sr_present = findKeysByValue(react.sources, 1, operator.ge)
        sr_absent = findKeysByValue(react.sources, 0, operator.eq)
        if ((set(sourceIn) <= set(sr_present)) & (set(sourceNotIn) <= set(sr_absent))):
            getattr(supermodel.reactions, name).update({react.id: react})


def getVennSegments(supermodel: SuperModel, Nletter=1):
    """ Getting metabolites and reactions networks for each Venn segment in Venn diagram. """
    combinations = []
    for i in range(1, len(supermodel.sources)):
        combinations.extend(itertools.combinations(supermodel.sources, i))
    for combo in combinations:
        sYes = list(combo)
        sNo = list(set(supermodel.sources) - set(combo))
        getDifference(supermodel, sYes, sNo, Nletter)


def defineNodeColor(colordata: dict, pallitra: str, id_mr: str, objects: SetofNewReactions or SetofNewMetabolites,
                    Nletter=1) -> [str, str]:
    for y in dir(objects):
        if y.startswith("Yes_") or y.startswith("core"):
            if id_mr in getattr(objects, y).keys():
                name = y.split("_No_")[0].replace("Yes_", "")
                if y.startswith("core"):
                    col = colordata.get(pallitra)[-1]
                else:
                    col = colordata.get(pallitra)[int(len(name) / Nletter) - 1]
                return ([col, id_mr + "\n" + name])


def defineEdgeColor(colordata: dict, pallitra: str, mr: NewObject, connections: dict, Nletter=1) -> [str, str]:
    for y in connections.keys():
        if y.startswith("Yes_") or y.startswith("core"):
            if mr in connections.get(y):
                name = y.split("_No_")[0].replace("Yes_", "")
                if y.startswith("core"):
                    col = colordata.get(pallitra)[-1]
                else:
                    col = colordata.get(pallitra)[int(len(name) / Nletter) - 1]
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
                    colname_r = defineNodeColor(colorBrewer, "purples", r_id, supermodel.reactions, Nletter)
                    if rea.id in pathway["metabolites"]:
                        colname_rea = defineNodeColor(colorBrewer, "reds", rea.id, supermodel.metabolites, Nletter)
                        cn_rea = defineEdgeColor(colorBrewer, "reds", rea, r.reactants, Nletter)
                    else:
                        colname_rea = defineNodeColor(colorBrewer, "blues", rea.id, supermodel.metabolites, Nletter)
                        cn_rea = defineEdgeColor(colorBrewer, "blues", rea, r.reactants, Nletter)

                    if pro.id in pathway["metabolites"]:
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
                                    colname_rn = defineNodeColor(colorBrewer, "greens", ra.id, supermodel.reactions,
                                                                 Nletter)
                                    colname_rean = defineNodeColor(colorBrewer, "blues", rea1.id,
                                                                   supermodel.metabolites, Nletter)
                                    colname_pron = defineNodeColor(colorBrewer, "blues", pro1.id,
                                                                   supermodel.metabolites, Nletter)
                                    cn_rean = defineEdgeColor(colorBrewer, "blues", rea1, ra.reactants, Nletter)
                                    cn_pron = defineEdgeColor(colorBrewer, "blues", pro1, ra.products, Nletter)
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


def drawBiomass(supermodel, name, only_difference=False, not_converted=False, colorBrewer = None, directed=True, wid=2000,
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
