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


def drawPathway(supermodel, pathway, met_not_int, colorBrewer, name, aminoacids=None, directed=False, surrounding=False, Nletter=1):
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
    pyvis_graph.show("../Output/"+name+".html")
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
    pyvis_graph.show("../Output/"+name+".html")
    return g


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
    pyvis_graph.show("../Output/"+name+".html")
    return g


def drawCore(supermodel, met_not_int, colorBrewer, name, directed=False, wid=2000, hei=1000):
    g = nx.DiGraph()
    for r in supermodel.reactions.core4.values():
        for rea in r.reactants.get("core4"):
            tmp_rea = rea.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
            for pro in r.products.get("core4"):
                tmp_pro = pro.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
                if (tmp_rea not in met_not_int) & (tmp_pro not in met_not_int):
                    g.add_node(r.id, shape="box", color=colorBrewer.get("purples")[-1])
                    g.add_node(rea.id, shape="o", color=colorBrewer.get("reds")[-1])
                    g.add_node(pro.id, shape="o", color=colorBrewer.get("reds")[-1])
                    g.add_edge(rea.id, r.id, color=colorBrewer.get("purples")[-1], font_color="black")
                    g.add_edge(r.id, pro.id, color=colorBrewer.get("purples")[-1], font_color="black")

    pyvis_graph = Network(width='{}px'.format(wid), height='{}px'.format(hei), directed=directed, notebook=False)
    pyvis_graph.from_nx(g)
    pyvis_graph.show("../Output/"+name+".html")
    return g