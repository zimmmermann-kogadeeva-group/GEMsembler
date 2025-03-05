import itertools
import operator
import sys
from collections import Counter

import numpy
from scipy.stats import mode

from .general import findKeysByValue


def getCoreConnections(
    connections: dict, core_size: int, compare_operator: operator, sources: [str]
) -> [str]:
    """ Getting connections (reactants/products/ for reaction or reactions for
    metabolites or genes for reactions and vv)
     that are present in more than core_size sources = original models. """
    all_connections = []
    for s in sources:
        connection = connections.get(s)
        for oneconnect in connection:
            all_connections.append(oneconnect)
    counted_connection = Counter(all_connections)
    selected_connection = []
    for key, value in counted_connection.items():
        if compare_operator(value, core_size):
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


def getCoreCoefficients(
    metabolites: dict,
    reactants: dict,
    products: dict,
    core_name: str,
    core_size: int,
    sources: [str],
) -> dict:
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


def getCoreGPR(
    gprs: dict,
    core_size: int,
    compare_operator: operator,
    sources: [str],
    and_as_solid: bool,
) -> [str]:
    """ Getting logical (or) parts of gene_reaction_rules (...and...)or(...)
    for reaction that are present in more than core_size sources = original models.
    While whether consider genes in (...and...) as solid thing is controlled
    by binary variable"""
    if not and_as_solid:
        s_comb = []
        if compare_operator == operator.ge:
            for i in range(core_size, len(sources) + 1):
                s_comb = s_comb + list(itertools.combinations(sources, i))
        elif compare_operator == operator.eq:
            s_comb = s_comb + list(itertools.combinations(sources, core_size))
        ands_selected = []
        original_ands = []
        for comb in s_comb:
            all_gpr_ands = []
            for s in comb:
                s_gpr_ands = []
                if gprs.get(s):
                    gpr = gprs.get(s)[0]
                    for gene_and in gpr.split(" or "):
                        s_gpr_ands.append(
                            gene_and.replace("(", "").replace(")", "").split(" and ")
                        )
                else:
                    s_gpr_ands.append([])
                all_gpr_ands.append(s_gpr_ands)
            ands_comb = list(itertools.product(*all_gpr_ands))
            for ands in ands_comb:
                intersect = set(ands[0])
                for a in ands:
                    intersect = set(intersect) & set(a)
                    original_ands.append(sorted(a))
                if intersect:
                    ands_selected.append(list(intersect))
        ands_selected_uniq = [
            sorted(list(x)) for x in set(tuple(x) for x in ands_selected)
        ]
        artificial_ands = [y for y in ands_selected_uniq if y not in original_ands]
        ands_selected_simple = []
        for u in ands_selected_uniq:
            inside = False
            if u in artificial_ands:
                for o in ands_selected_uniq:
                    if (u != o) and (set(u) & set(o) == set(u)):
                        inside = True
            if not inside:
                ands_selected_simple.append(u)
        selected_gpr_ands = []
        for gpr_ands in ands_selected_simple:
            if len(gpr_ands) > 1 and len(ands_selected_simple) > 1:
                selected_gpr_ands.append("(" + " and ".join(sorted(gpr_ands)) + ")")
            else:
                selected_gpr_ands.append(" and ".join(sorted(gpr_ands)))
        selected_gpr = " or ".join(sorted(selected_gpr_ands))
    else:
        all_gpr_ands = []
        for s in sources:
            if gprs.get(s):
                gpr = gprs.get(s)[0]
                for gene_and in gpr.split(" or "):
                    all_gpr_ands.append(gene_and)
        counted_gpr_ands = Counter(all_gpr_ands)
        selected_gpr_ands = []
        for key, value in counted_gpr_ands.items():
            if compare_operator(value, core_size):
                selected_gpr_ands.append(key)
        selected_gpr = " or ".join(sorted(selected_gpr_ands))
    if selected_gpr:
        return [selected_gpr]
    else:
        return []


def getCore(
    supermodel, core_size: int, compare_operator: operator, and_as_solid: bool,
):
    """ Getting supermodel core: intersection of at least core_size amount of sources (by default, intersection of all
     sources). Getting supermodel union of all sources. """
    if compare_operator == operator.ge:
        coreN = "core" + str(core_size)
    elif compare_operator == operator.eq:
        coreN = "In" + str(core_size)
    else:
        raise ValueError(
            "Comparison operator is not supported. Has to be operator.ge (>=) or operator.eq (==)"
        )
    for met in supermodel.metabolites.assembly.values():
        core_r = getCoreConnections(
            met.reactions, core_size, compare_operator, supermodel.sources
        )
        met.reactions["comparison"].update({coreN: core_r})
        if compare_operator(met.in_models["models_amount"], core_size):
            supermodel.metabolites.comparison[coreN].update({met.id: met})
    for gene in supermodel.genes.assembly.values():
        core_rg = getCoreConnections(
            gene.reactions, core_size, compare_operator, supermodel.sources
        )
        gene.reactions["comparison"].update({coreN: core_rg})
        if compare_operator(gene.in_models["models_amount"], core_size):
            supermodel.genes.comparison[coreN].update({gene.id: gene})
    for react in supermodel.reactions.assembly.values():
        core_reactants = getCoreConnections(
            react.reactants, core_size, compare_operator, supermodel.sources
        )
        core_products = getCoreConnections(
            react.products, core_size, compare_operator, supermodel.sources
        )
        core_genes = getCoreConnections(
            react.genes, core_size, compare_operator, supermodel.sources
        )
        core_gpr = getCoreGPR(
            react.gene_reaction_rule,
            core_size,
            compare_operator,
            supermodel.sources,
            and_as_solid,
        )
        core_lower_bound = getCoreLowerBounds(
            react.lower_bound, core_size, react.in_models["models_list"]
        )
        core_upper_bound = getCoreUpperBounds(
            react.upper_bound, core_size, react.in_models["models_list"]
        )
        react.reactants["comparison"].update({coreN: core_reactants})
        react.products["comparison"].update({coreN: core_products})
        react.genes["comparison"].update({coreN: core_genes})
        react.gene_reaction_rule["comparison"].update({coreN: core_gpr})
        react.lower_bound["comparison"].update({coreN: core_lower_bound})
        react.upper_bound["comparison"].update({coreN: core_upper_bound})
        core_metabolites = getCoreCoefficients(
            react.metabolites,
            react.reactants["comparison"],
            react.products["comparison"],
            coreN,
            core_size,
            react.in_models["models_list"],
        )
        react.metabolites["comparison"].update({coreN: core_metabolites})
        if compare_operator(react.in_models["models_amount"], core_size):
            supermodel.reactions.comparison[coreN].update({react.id: react})
    return coreN


def getDifConnections(connections: dict, sourceIn: [str], sourceNotIn: [str]):
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


def getSomeCoefficients(
    metabolites: dict, reactants: dict, products: dict, name: str, sourceIn: [str]
):
    """ Getting coefficients mode of metabolites from sourceIn """
    coefficients = {}
    if (reactants["comparison"].get(name)) or (products["comparison"].get(name)):
        for rea in reactants["comparison"].get(name):
            k_mod = []
            for sI in sourceIn:
                if metabolites.get(sI).get(rea):
                    k_mod.append(metabolites.get(sI).get(rea))
            coefficients.update({rea: mode(k_mod, keepdims=False)[0]})
        for pro in products["comparison"].get(name):
            k_mod = []
            for sI in sourceIn:
                if metabolites.get(sI).get(pro):
                    k_mod.append(metabolites.get(sI).get(pro))
            coefficients.update({pro: mode(k_mod, keepdims=False)[0]})
    return coefficients


def getDifGPR(gprs: dict, sourceIn: [str], sourceNotIn: [str], and_as_solid: bool):
    """ Getting logical (or) parts of gene_reaction_rules (...and...)or(...) that are present in "sourceIn"
    list of sources = original models and not present in "sourceNotIn" list of sources = original models.
    While whether consider genes in (...and...) as solid thing is controlled by binary variable. """
    if not and_as_solid:
        ands_selected = []
        original_ands = []
        In_gpr_ands = []
        for s in sourceIn:
            s_gpr_ands = []
            if gprs.get(s):
                gpr = gprs.get(s)[0]
                for gene_and in gpr.split(" or "):
                    s_gpr_ands.append(
                        gene_and.replace("(", "").replace(")", "").split(" and ")
                    )
            else:
                s_gpr_ands.append([])
            In_gpr_ands.append(s_gpr_ands)
        notIn_gpr_ands = []
        for ss in sourceNotIn:
            if gprs.get(ss):
                gpr = gprs.get(ss)[0]
                for gene_and in gpr.split(" or "):
                    notIn_gpr_ands.append(
                        gene_and.replace("(", "").replace(")", "").split(" and ")
                    )
        ands_comb = list(itertools.product(*In_gpr_ands))
        for ands in ands_comb:
            intersect = set(ands[0])
            for a in ands:
                intersect = set(intersect) & set(a)
                original_ands.append(sorted(a))
            for ng in notIn_gpr_ands:
                intersect = set(intersect) - set(ng)
            if intersect:
                ands_selected.append(list(intersect))
        ands_selected_uniq = [
            sorted(list(x)) for x in set(tuple(x) for x in ands_selected)
        ]
        artificial_ands = [y for y in ands_selected_uniq if y not in original_ands]
        ands_selected_simple = []
        for u in ands_selected_uniq:
            inside = False
            if u in artificial_ands:
                for o in ands_selected_uniq:
                    if (u != o) and (set(u) & set(o) == set(u)):
                        inside = True
            if not inside:
                ands_selected_simple.append(u)
        selected_gpr_ands = []
        for gpr_ands in ands_selected_simple:
            if len(gpr_ands) > 1 and len(ands_selected_simple) > 1:
                selected_gpr_ands.append("(" + " and ".join(sorted(gpr_ands)) + ")")
            else:
                selected_gpr_ands.append(" and ".join(sorted(gpr_ands)))
        dif_gpr = " or ".join(sorted(selected_gpr_ands))
    else:
        if gprs.get(sourceIn[0]):
            gpr_and_In = set(gprs.get(sourceIn[0])[0].split(" or "))
            for sIn in sourceIn:
                if gprs.get(sIn):
                    gpr_and_In = gpr_and_In & set(gprs.get(sIn)[0].split(" or "))
            gpr_and_NotIn = []
            for sNotIn in sourceNotIn:
                if gprs.get(sNotIn):
                    gpr_and_NotIn = gpr_and_NotIn + gprs.get(sNotIn)[0].split(" or ")
            dif_gpr = set(gpr_and_In) - set(gpr_and_NotIn)
        else:
            dif_gpr = []
        dif_gpr = " or ".join(list(sorted(dif_gpr)))
    if dif_gpr:
        return [dif_gpr]
    else:
        return []


def getDifference(
    supermodel, sourceIn: [str], sourceNotIn: [str], and_as_solid: bool, nletter: int,
):
    """ Getting metabolites and reactions that are present in "sourceIn"
    list of sources = original models
    and not present in "sourceNotIn" list of sources = original models. """
    if sourceIn:
        name = "Yes_"
        for sI in sorted(sourceIn):
            name = name + sI[:nletter]
        if sourceNotIn:
            name = name + "_No_"
            for sNI in sorted(sourceNotIn):
                name = name + sNI[:nletter]
    else:
        name = "No_"
        for sNI in sourceNotIn:
            name = name + sNI[:nletter]
    for met in supermodel.metabolites.assembly.values():
        sm_present = findKeysByValue(met.sources, 1, operator.ge)
        sm_absent = findKeysByValue(met.sources, 0, operator.eq)
        if sourceIn:
            dif_r = getDifConnections(met.reactions, sourceIn, sourceNotIn)
        else:
            dif_r = getDifConnections(met.reactions, sm_present, sourceNotIn)
        met.reactions["comparison"].update({name: dif_r})
        if (set(sourceIn) <= set(sm_present)) & (set(sourceNotIn) <= set(sm_absent)):
            supermodel.metabolites.comparison[name].update({met.id: met})
    for gene in supermodel.genes.assembly.values():
        sg_present = findKeysByValue(gene.sources, 1, operator.ge)
        sg_absent = findKeysByValue(gene.sources, 0, operator.eq)
        if sourceIn:
            dif_rg = getDifConnections(gene.reactions, sourceIn, sourceNotIn)
        else:
            dif_rg = getDifConnections(gene.reactions, sg_present, sourceNotIn)
        gene.reactions["comparison"].update({name: dif_rg})
        if (set(sourceIn) <= set(sg_present)) & (set(sourceNotIn) <= set(sg_absent)):
            supermodel.genes.comparison[name].update({gene.id: gene})
    for react in supermodel.reactions.assembly.values():
        sr_present = findKeysByValue(react.sources, 1, operator.ge)
        sr_absent = findKeysByValue(react.sources, 0, operator.eq)
        if sourceIn:
            dif_reactants = getDifConnections(react.reactants, sourceIn, sourceNotIn)
            dif_products = getDifConnections(react.products, sourceIn, sourceNotIn)
            dif_genes = getDifConnections(react.genes, sourceIn, sourceNotIn)
            dif_gpr = getDifGPR(
                react.gene_reaction_rule, sourceIn, sourceNotIn, and_as_solid
            )
            sI_lower_bound = getSomeBound(react.lower_bound, "lower", sourceIn)
            sI_upper_bound = getSomeBound(react.upper_bound, "upper", sourceIn)
        else:
            dif_reactants = getDifConnections(react.reactants, sr_present, sourceNotIn)
            dif_products = getDifConnections(react.products, sr_present, sourceNotIn)
            dif_genes = getDifConnections(react.genes, sr_present, sourceNotIn)
            dif_gpr = getDifGPR(
                react.gene_reaction_rule, sr_present, sourceNotIn, and_as_solid
            )
            sI_lower_bound = getSomeBound(react.lower_bound, "lower", sr_present)
            sI_upper_bound = getSomeBound(react.upper_bound, "upper", sr_present)
        react.reactants["comparison"].update({name: dif_reactants})
        react.products["comparison"].update({name: dif_products})
        react.genes["comparison"].update({name: dif_genes})
        react.gene_reaction_rule["comparison"].update({name: dif_gpr})
        react.lower_bound["comparison"].update({name: sI_lower_bound})
        react.upper_bound["comparison"].update({name: sI_upper_bound})
        if sourceIn:
            sI_metabolites = getSomeCoefficients(
                react.metabolites, react.reactants, react.products, name, sourceIn
            )
        else:
            sI_metabolites = getSomeCoefficients(
                react.metabolites, react.reactants, react.products, name, sr_present
            )
        react.metabolites["comparison"].update({name: sI_metabolites})
        if (set(sourceIn) <= set(sr_present)) & (set(sourceNotIn) <= set(sr_absent)):
            supermodel.reactions.comparison[name].update({react.id: react})
    if name not in supermodel.metabolites.comparison.keys():
        supermodel.metabolites.comparison.update({name: {}})
    if name not in supermodel.reactions.comparison.keys():
        supermodel.reactions.comparison.update({name: {}})
    if name not in supermodel.genes.comparison.keys():
        supermodel.genes.comparison.update({name: {}})
    return name
