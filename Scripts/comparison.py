import itertools
import operator
from collections import Counter
import numpy
from scipy.stats import mode
from creation import SuperModel, NewObject, SetofNewReactions, SetofNewMetabolites
from general import findKeysByValue


def getCoreConnections(connections: dict, core_size: int, sources: [str]) -> [str]:
    """ Getting connections (reactants/products/ for reaction or reactions for metabolites or genes for reactions and vv)
     that are present in more then core_size sources = original models. """
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


def getCoreUpperBounds(bounds: dict, core_size: int, sources: [str]) -> [str]:
    """ Getting upper bounds via uniting all possible intersections of core_size number sources """
    if len(sources) < core_size:
        return []
    combinations = list(itertools.combinations(sources, core_size))
    union_bound = 0
    for combination in combinations:
        intersection_bound = 1000
        for c in combination:
            if intersection_bound > bounds.get(c)[0]:
                intersection_bound = bounds.get(c)[0]
        if union_bound < intersection_bound:
            union_bound = intersection_bound
    return [union_bound]


def getCoreLowerBounds(bounds: dict, core_size: int, sources: [str]) -> [str]:
    """ Getting lower bounds via uniting all possible intersections of core_size number sources """
    if len(sources) < core_size:
        return []
    combinations = list(itertools.combinations(sources, core_size))
    union_bound = 0
    for combination in combinations:
        intersection_bound = -1000
        for c in combination:
            if intersection_bound < bounds.get(c)[0]:
                intersection_bound = bounds.get(c)[0]
        if union_bound > intersection_bound:
            union_bound = intersection_bound
    return [union_bound]


def getCoreCoefficients(metabolites: dict, reactants: dict, products: dict, core_name: str, core_size: int,
                        sources: [str]) -> dict:
    """ Getting core coefficients for metabolites via average of all possible modes of core_size number sources """
    core_metabolites = {}
    if len(sources) >= core_size:
        combinations = list(itertools.combinations(sources, core_size))
        for rea in reactants.get(core_name):
            k_mean = []
            for combination in combinations:
                k_mod = []
                for c in combination:
                    if metabolites.get(c).get(rea):
                        k_mod.append(metabolites.get(c).get(rea))
                if k_mod:
                    k_mean.append(mode(k_mod, keepdims=False)[0])
            k = numpy.mean(k_mean)
            core_metabolites.update({rea: k})
        for pro in products.get(core_name):
            k_mean = []
            for combination in combinations:
                k_mod = []
                for c in combination:
                    if metabolites.get(c).get(pro):
                        k_mod.append(metabolites.get(c).get(pro))
                if k_mod:
                    k_mean.append(mode(k_mod, keepdims=False)[0])
            k = numpy.mean(k_mean)
            core_metabolites.update({pro: k})
    return core_metabolites


def getCore(supermodel: SuperModel, core_size):
    """ Getting supermodel core: intersection of at least core_size amount of sources (by default, intersection of all
     sources). Getting supermodel union of all sources. """
    coreN = "core" + str(core_size)
    for met in supermodel.metabolites.converted.values():
        tmp_models = findKeysByValue(met.sources, 1, operator.ge)
        core_r = getCoreConnections(met.reactions, core_size, supermodel.sources)
        met.reactions.update({coreN: core_r})
        if len(tmp_models) >= core_size:
            if core_size == 1:
                supermodel.metabolites.converted.update({met.id: met})
            else:
                getattr(supermodel.metabolites, coreN).update({met.id: met})
    for gene in supermodel.genes.converted.values():
        tmp_models = findKeysByValue(gene.sources, 1, operator.ge)
        core_rg = getCoreConnections(gene.reactions, core_size, supermodel.sources)
        gene.reactions.update({coreN: core_rg})
        if len(tmp_models) >= core_size:
            if core_size == 1:
                supermodel.genes.converted.update({gene.id: gene})
            else:
                getattr(supermodel.genes, coreN).update({gene.id: gene})
    for react in supermodel.reactions.converted.values():
        tmp_models = findKeysByValue(react.sources, 1, operator.ge)
        core_reactants = getCoreConnections(react.reactants, core_size, supermodel.sources)
        core_products = getCoreConnections(react.products, core_size, supermodel.sources)
        core_genes = getCoreConnections(react.genes, core_size, supermodel.sources)
        core_lower_bound = getCoreLowerBounds(react.lower_bound, core_size, tmp_models)
        core_upper_bound = getCoreUpperBounds(react.upper_bound, core_size, tmp_models)
        react.reactants.update({coreN: core_reactants})
        react.products.update({coreN: core_products})
        react.genes.update({coreN: core_genes})
        react.lower_bound.update({coreN: core_lower_bound})
        react.upper_bound.update({coreN: core_upper_bound})
        core_metabolites = getCoreCoefficients(react.metabolites, react.reactants, react.products, coreN, core_size,
                                               tmp_models)
        react.metabolites.update({coreN: core_metabolites})
        if len(tmp_models) >= core_size:
            if core_size == 1:
                supermodel.reactions.converted.update({react.id: react})
            else:
                getattr(supermodel.reactions, coreN).update({react.id: react})


def getDifConnections(connections: dict, sourceIn: [str], sourceNotIn: [str]) -> [NewObject]:
    """ Getting connections (reactants/products/ for reaction or reactions for metabolites or genes for reactions and vv)
     that are present in "sourceIn" list of sources = original models and not present in "sourceNotIn"
     list of sources = original models. """
    connectionIn = set(connections.get(sourceIn[0]))
    for sIn in sourceIn:
        connectionIn = connectionIn & set(connections.get(sIn))
    connectionNotIn = []
    for sNotIn in sourceNotIn:
        connectionNotIn = connectionNotIn + connections.get(sNotIn)
    difConnection = set(connectionIn) - set(connectionNotIn)
    return list(difConnection)


def getSomeBound(bounds: dict, bounds_type: str, sourceIn: [str]):
    """ Getting intersection of lower/upper bounds from sourceIn """
    if bounds_type == "lower":
        bound = -1000
    if bounds_type == "upper":
        bound = 1000
    for sI in sourceIn:
        if not bounds.get(sI):
            return []
        else:
            if bounds_type == "lower":
                if bound < bounds.get(sI)[0]:
                    bound = bounds.get(sI)[0]
            if bounds_type == "upper":
                if bound > bounds.get(sI)[0]:
                    bound = bounds.get(sI)[0]
    return [bound]


def getSomeCoefficients(metabolites: dict, reactants: dict, products: dict, name: str, sourceIn: [str]):
    """ Getting coefficients mode of metabolites from sourceIn """
    coefficients = {}
    if (reactants.get(name)) or (products.get(name)):
        for rea in reactants.get(name):
            k_mod = []
            for sI in sourceIn:
                if metabolites.get(sI).get(rea):
                    k_mod.append(metabolites.get(sI).get(rea))
            coefficients.update({rea: mode(k_mod, keepdims=False)[0]})
        for pro in products.get(name):
            k_mod = []
            for sI in sourceIn:
                if metabolites.get(sI).get(pro):
                    k_mod.append(metabolites.get(sI).get(pro))
            coefficients.update({pro: mode(k_mod, keepdims=False)[0]})
    return coefficients


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
    setattr(supermodel.genes, name, {})
    for met in supermodel.metabolites.converted.values():
        dif_r = getDifConnections(met.reactions, sourceIn, sourceNotIn)
        met.reactions.update({name: dif_r})
        sm_present = findKeysByValue(met.sources, 1, operator.ge)
        sm_absent = findKeysByValue(met.sources, 0, operator.eq)
        if ((set(sourceIn) == set(sm_present)) & (set(sourceNotIn) == set(sm_absent))):
            getattr(supermodel.metabolites, name).update({met.id: met})
    for gene in supermodel.genes.converted.values():
        dif_rg = getDifConnections(gene.reactions, sourceIn, sourceNotIn)
        gene.reactions.update({name: dif_rg})
        sg_present = findKeysByValue(gene.sources, 1, operator.ge)
        sg_absent = findKeysByValue(gene.sources, 0, operator.eq)
        if ((set(sourceIn) == set(sg_present)) & (set(sourceNotIn) == set(sg_absent))):
            getattr(supermodel.genes, name).update({gene.id: gene})
    for react in supermodel.reactions.converted.values():
        dif_reactants = getDifConnections(react.reactants, sourceIn, sourceNotIn)
        dif_products = getDifConnections(react.products, sourceIn, sourceNotIn)
        dif_genes = getDifConnections(react.genes, sourceIn, sourceNotIn)
        sI_lower_bound = getSomeBound(react.lower_bound, "lower", sourceIn)
        sI_upper_bound = getSomeBound(react.upper_bound, "upper", sourceIn)
        react.reactants.update({name: dif_reactants})
        react.products.update({name: dif_products})
        react.genes.update({name: dif_genes})
        react.lower_bound.update({name: sI_lower_bound})
        react.upper_bound.update({name: sI_upper_bound})
        sI_metabolites = getSomeCoefficients(react.metabolites, react.reactants, react.products, name, sourceIn)
        react.metabolites.update({name: sI_metabolites})
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
        sYes = sorted(list(combo))
        sNo = sorted((list(set(supermodel.sources) - set(combo))))
        getDifference(supermodel, sYes, sNo, Nletter)


def swapReactantsAndProducts(r: NewObject, sources_present: list, sources_to_swap: list, Nletter):
    for tmp in sources_present:
        if tmp[:Nletter] in sources_to_swap:
            a = r.reactants.get(tmp)
            b = r.products.get(tmp)
            r.reactants[tmp] = b
            r.products[tmp] = a
            aa = r.lower_bound.get(tmp)[0] * -1
            bb = r.upper_bound.get(tmp)[0] * -1
            r.lower_bound[tmp] = [bb]
            r.upper_bound[tmp] = [aa]
            for met, koef in r.metabolites.get(tmp).items():
                r.metabolites.get(tmp)[met] = koef * -1


def getSwitchedMetabolites(supermodel: SuperModel, Nletter=1):
    for d in dir(supermodel.reactions):
        if (d.startswith("Yes")) | (d.startswith("core")):
            for r in getattr(supermodel.reactions, d).values():
                tmp_models = findKeysByValue(r.sources, 1, operator.ge)
                ex = False
                for tmp in tmp_models:
                    if (not r.reactants.get(tmp)) | (not r.products.get(tmp)):
                        ex = True
                if not ex:
                    if (not r.reactants.get(d)) | (not r.products.get(d)):
                        max_consist = 0
                        consist = []
                        all_consist = []
                        for key, value in r.reactants.items():
                            if (value != []) & (key.startswith("Yes_")):
                                tmp_consist = key.removeprefix("Yes_").split("_No_")[0]
                                all_consist.append(tmp_consist)
                                if max_consist < len(tmp_consist) / Nletter:
                                    max_consist = len(tmp_consist) / Nletter
                                    consist = [tmp_consist]
                                elif max_consist == len(tmp_consist) / Nletter:
                                    consist.append(tmp_consist)
                        if len(consist) == 1:
                            # "Case 1: majority"
                            start_yes_models = [consist[0][i:i + Nletter] for i in range(0, len(consist[0]), Nletter)]
                            start_no_models = [tmp[:Nletter] for tmp in tmp_models if
                                               tmp[:Nletter] not in start_yes_models]
                            swapReactantsAndProducts(r, tmp_models, start_no_models, Nletter)
                        elif len(consist) == 2:
                            lb1 = 0
                            lb2 = 0
                            for tmp in tmp_models:
                                if tmp[:Nletter] in consist[0]:
                                    if r.lower_bound.get(tmp)[0] < lb1:
                                        lb1 = r.lower_bound.get(tmp)[0]
                                if tmp[:Nletter] in consist[1]:
                                    if r.lower_bound.get(tmp)[0] < lb2:
                                        lb2 = r.lower_bound.get(tmp)[0]
                            swap = None
                            if (lb1 >= 0) & (lb2 < 0):
                                swap = consist[1]
                            if (lb1 < 0) & (lb2 >= 0):
                                swap = consist[0]
                            if swap:
                                # "Case 2: boundary"
                                swap = [swap[i:i + Nletter] for i in range(0, len(swap), Nletter)]
                                swapReactantsAndProducts(r, tmp_models, swap, Nletter)
                            else:
                                # "Case 3: Nothing sort"
                                sel = sorted(tmp_models)[0]
                                not_sel = []
                                for tmp in sorted(tmp_models)[1:]:
                                    if not (set(r.reactants.get(tmp)) & set(r.reactants.get(sel))):
                                        not_sel.append(tmp[:Nletter])
                                swapReactantsAndProducts(r, tmp_models, not_sel, Nletter)
                        else:
                            # "Case 3: Nothing sort"
                            sel = sorted(tmp_models)[0]
                            not_sel = []
                            for tmp in sorted(tmp_models)[1:]:
                                if not (set(r.reactants.get(tmp)) & set(r.reactants.get(sel))):
                                    not_sel.append(tmp[:Nletter])
                            swapReactantsAndProducts(r, tmp_models, not_sel, Nletter)


def runComparison(supermodel: SuperModel, run_all=True, core_size=None, sYes=None, sNo=None, union_size=1, Nletter=1):
    if run_all:
        for source in supermodel.sources:
            setattr(supermodel.metabolites, source, {})
            setattr(supermodel.reactions, source, {})
            setattr(supermodel.genes, source, {})
        for met in supermodel.metabolites.converted.values():
            tmp_models = findKeysByValue(met.sources, 1, operator.ge)
            for tmp_source in tmp_models:
                getattr(supermodel.metabolites, tmp_source).update({met.id: met})
        for r in supermodel.reactions.converted.values():
            tmp_models = findKeysByValue(r.sources, 1, operator.ge)
            for tmp_source in tmp_models:
                getattr(supermodel.reactions, tmp_source).update({r.id: r})
        for g in supermodel.genes.converted.values():
            tmp_models = findKeysByValue(g.sources, 1, operator.ge)
            for tmp_source in tmp_models:
                getattr(supermodel.genes, tmp_source).update({g.id: g})
        if not core_size:
            core_size = len(supermodel.sources)
        coreN = "core" + str(core_size)
        setattr(supermodel.metabolites, coreN, {})
        setattr(supermodel.reactions, coreN, {})
        setattr(supermodel.genes, coreN, {})
        getCore(supermodel, core_size)
        getCore(supermodel, union_size)
        getVennSegments(supermodel, Nletter)
        getSwitchedMetabolites(supermodel, Nletter)
        getCore(supermodel, core_size)
        getCore(supermodel, union_size)
        getVennSegments(supermodel, Nletter)
    else:
        if core_size:
            coreN = "core" + str(core_size)
            setattr(supermodel.metabolites, coreN, {})
            setattr(supermodel.reactions, coreN, {})
            setattr(supermodel.genes, coreN, {})
            getCore(supermodel, core_size)
        if sYes and sNo:
            getDifference(supermodel, sYes, sNo, Nletter)
