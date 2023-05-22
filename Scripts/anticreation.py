from cobra import Model, Reaction, Metabolite
import tempfile
from pprint import pprint
from cobra.io import write_sbml_model, validate_sbml_model
from creation import SuperModel, NewObject


def getModelOfInterest(supermodel: SuperModel, interest_level: str, write_sbml=True, name=None,
                       reactions_include: [NewObject] = None,
                       reactions_exclude: [NewObject] = None):
    """ Creating COBRA model from supermodel based on specific level of interest for example core or union.
    Additionaly, some reactions"""
    outmodel = Model(interest_level)
    if interest_level == "union1":
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
            interest_level = "union1"
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
        out_reaction.lower_bound = r.lower_bound.get(interest_level)[0]
        out_reaction.upper_bound = r.upper_bound.get(interest_level)[0]
        for met, k in r.metabolites.get(interest_level).items():
            out_met = Metabolite(met.id, name=met.name, compartment=met.id[-1])
            out_reaction.add_metabolites({out_met: k})
        outmodel.add_reactions([out_reaction])
    outmodel.objective = "Biomass"
    if write_sbml:
        write_sbml_model(outmodel, "../Output/" + name)
    report = validate_sbml_model(filename="../Output/" + name)
    pprint(report)
    return outmodel
