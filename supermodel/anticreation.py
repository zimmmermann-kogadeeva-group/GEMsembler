from cobra import Model, Reaction, Metabolite
import tempfile
from pprint import pprint
from cobra.io import write_sbml_model, validate_sbml_model
from creation import SuperModel, NewObject


def gapfillTransportR(cobra_model: Model, supermodel: SuperModel):
    for exchange in cobra_model.exchanges:
        met_e = list(exchange.reactants)[0].id
        met_c = list(exchange.reactants)[0].id[:-1]+"c"
        cobra_transport = False
        for r in cobra_model.metabolites.get_by_id(met_e).reactions:
            r_react = [m.id for m in r.reactants]
            r_pro = [m.id for m in r.products]
            if ((met_e in r_react) & (met_c in r_pro)) | ((met_e in r_pro) & (met_c in r_react)):
                cobra_transport = True
        if not cobra_transport:
            if met_c in supermodel.metabolites.converted.keys():
                transport_r = []
                for r_super in supermodel.metabolites.converted.get(met_e).reactions.get("core1"):
                    rs_react = [m.id for m in r_super.reactants.get("core1")]
                    rs_pro = [m.id for m in r_super.products.get("core1"
                                                                 "")]
                    if ((met_e in rs_react) & (met_c in rs_pro)) | ((met_e in rs_pro) & (met_c in rs_react)):
                        transport_r.append(r_super)
                for tr in transport_r:
                    tr_r = Reaction(tr.id)
                    if tr.name:
                        tr_r.name = tr.name
                    else:
                        tr_r.name = ""
                    out_subsystem = ""
                    for source in supermodel.sources:
                        if tr.subsystem.get(source):
                            tmp = tr.subsystem.get(source)[0]
                        else:
                            tmp = "NaN"
                        out_subsystem = out_subsystem + "#" + source + "#" + tmp
                    tr_r.subsystem = out_subsystem
                    tr_r.lower_bound = tr.lower_bound.get("core1")[0]
                    tr_r.upper_bound = tr.upper_bound.get("core1")[0]
                    for met, k in tr.metabolites.get("core1").items():
                        tr_met = Metabolite(met.id, name=met.name, compartment=met.id[-1])
                        tr_r.add_metabolites({tr_met: k})
                    cobra_model.add_reactions([tr_r])


def getModelOfInterest(supermodel: SuperModel, interest_level: str, gene_interest_level=None, extend_zero_bounds=True,
                       gapfill_transport=True, write_sbml=True, name=None,
                       reactions_include: [NewObject] = None,
                       reactions_exclude: [NewObject] = None):
    """ Creating COBRA model from supermodel based on specific level of interest for example core or union.
    Additionaly, some reactions"""
    if not gene_interest_level:
        gene_interest_level = interest_level
    outmodel = Model(interest_level)
    if interest_level == "core1":
        in_reactions = getattr(supermodel.reactions, "converted").values()
    else:
        in_reactions = getattr(supermodel.reactions, interest_level).values()
    if reactions_include:
        in_reactions = list(set(in_reactions) | set(reactions_include))
    else:
        reactions_include = []
    if reactions_exclude:
        in_reactions = list(set(in_reactions) - set(reactions_exclude))
    for r in in_reactions:
        if r in reactions_include:
            interest_level = "core1"
        out_reaction = Reaction(r.id)
        if r.name:
            out_reaction.name = r.name
        else:
            out_reaction.name = ""
        out_subsystem = ""
        for source in supermodel.sources:
            if r.subsystem.get(source):
                tmp = r.subsystem.get(source)[0]
            else:
                tmp = "NaN"
            out_subsystem = out_subsystem + "#" + source + "#" + tmp
        out_reaction.subsystem = out_subsystem
        if extend_zero_bounds and (r.upper_bound.get(interest_level)[0] - r.lower_bound.get(interest_level)[0] == 0):
            out_reaction.lower_bound = -1000.
            out_reaction.upper_bound = 1000.
        else:
            out_reaction.lower_bound = r.lower_bound.get(interest_level)[0]
            out_reaction.upper_bound = r.upper_bound.get(interest_level)[0]
        for met, k in r.metabolites.get(interest_level).items():
            out_met = Metabolite(met.id, name=met.name, compartment=met.id[-1])
            out_reaction.add_metabolites({out_met: k})
        if r.gene_reaction_rule.get(gene_interest_level):
            out_reaction.gene_reaction_rule = r.gene_reaction_rule.get(gene_interest_level)[0]
        else:
            out_reaction.gene_reaction_rule = ""
        outmodel.add_reactions([out_reaction])
    if supermodel.reactions.converted.get("Biomass") not in in_reactions:
        bio_r = Reaction("Biomass")
        bio_r.subsystem = "Growth"
        bio_r.lower_bound = supermodel.reactions.converted.get("Biomass").lower_bound.get(interest_level)[0]
        bio_r.upper_bound = supermodel.reactions.converted.get("Biomass").upper_bound.get(interest_level)[0]
        for met, k in supermodel.reactions.converted.get("Biomass").metabolites.get(interest_level).items():
            out_met = Metabolite(met.id, name=met.name, compartment=met.id[-1])
            bio_r.add_metabolites({out_met: k})
        outmodel.add_reactions([bio_r])
    outmodel.objective = "Biomass"
    if gapfill_transport:
        gapfillTransportR(outmodel, supermodel)
    if write_sbml:
        write_sbml_model(outmodel, "../Output/" + name)
    report = validate_sbml_model(filename="../Output/" + name)
    pprint(report)
    return outmodel
