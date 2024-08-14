import warnings
from copy import deepcopy
from pathlib import Path
from pprint import pprint

from cobra import Metabolite, Model, Reaction
from cobra.io import validate_sbml_model, write_sbml_model

from .creation import NewElement, SuperModel


def gapfill_transport_r(cobra_model: Model, supermodel: SuperModel):
    for exchange in cobra_model.exchanges:
        met_e = list(exchange.reactants)[0].id
        met_c = list(exchange.reactants)[0].id[:-1] + "c"
        cobra_transport = False
        for r in cobra_model.metabolites.get_by_id(met_e).reactions:
            r_react = [m.id for m in r.reactants]
            r_pro = [m.id for m in r.products]
            if ((met_e in r_react) & (met_c in r_pro)) | (
                (met_e in r_pro) & (met_c in r_react)
            ):
                cobra_transport = True
        if not cobra_transport:
            if met_c in supermodel.metabolites.assembly.keys():
                transport_r = []
                for r_super in supermodel.metabolites.assembly.get(met_e).reactions.get(
                    "assembly"
                ):
                    rs_react = [m.id for m in r_super.reactants.get("assembly")]
                    rs_pro = [m.id for m in r_super.products.get("assembly")]
                    if ((met_e in rs_react) & (met_c in rs_pro)) | (
                        (met_e in rs_pro) & (met_c in rs_react)
                    ):
                        transport_r.append(r_super)
                for tr in transport_r:
                    tr_r = Reaction(tr.id)
                    if tr.name:
                        tr_r.name = tr.name
                    else:
                        tr_r.name = ""
                    out_subsystem = ""
                    for source in tr.in_models["models_list"]:
                        out_subsystem = (
                            out_subsystem
                            + "#"
                            + source
                            + "#"
                            + tr.subsystem.get(source)[0]
                        )
                    tr_r.subsystem = out_subsystem
                    tr_r.lower_bound = tr.lower_bound.get("assembly")[0]
                    tr_r.upper_bound = tr.upper_bound.get("assembly")[0]
                    for met, k in tr.metabolites.get("assembly").items():
                        tr_met = Metabolite(
                            met.id, name=met.name, compartment=met.id[-1]
                        )
                        tr_r.add_metabolites({tr_met: k})
                    cobra_model.add_reactions([tr_r])


def get_model_of_interest(
    supermodel: SuperModel,
    interest_level: str,
    output_name=None,
    gene_interest_level=None,
    biomass_interest_level=None,
    simple_biomass_products=True,
    extend_zero_bounds=True,
    gapfill_transport=True,
    reactions_include: [NewElement] = None,
    reactions_exclude: [NewElement] = None,
):
    """ Creating COBRA model from supermodel based on specific level of interest for example core or union.
    Additionaly, some reactions"""
    if not gene_interest_level:
        gene_interest_level = interest_level
    if not biomass_interest_level:
        biomass_interest_level = interest_level
    outmodel = Model(interest_level)
    if interest_level in supermodel.sources + ["assembly"]:
        in_reactions = getattr(supermodel.reactions, interest_level).values()
    elif interest_level in supermodel.reactions.comparison.keys():
        in_reactions = supermodel.reactions.comparison[interest_level].values()
    else:
        raise ValueError(
            f"Interest level {interest_level} is not determined yet. "
            f"Please run corresponding supermodel comparison first."
        )
    if reactions_include:
        in_reactions = list(set(in_reactions) | set(reactions_include))
    else:
        reactions_include = []
    if reactions_exclude:
        in_reactions = list(set(in_reactions) - set(reactions_exclude))
    if supermodel.reactions.assembly.get("Biomass") in in_reactions:
        in_reactions = list(
            set(list(in_reactions)) - {supermodel.reactions.assembly.get("Biomass")}
        )
    for r in in_reactions:
        if r in reactions_include:
            interest_level_r = "assembly"
            gene_interest_level_r = "assembly"
        else:
            interest_level_r = interest_level
            gene_interest_level_r = gene_interest_level
        if interest_level_r in supermodel.sources + ["assembly"]:
            r_upper_bound = r.upper_bound
            r_lower_bound = r.lower_bound
            r_metabolites = r.metabolites
        else:
            r_upper_bound = r.upper_bound["comparison"]
            r_lower_bound = r.lower_bound["comparison"]
            r_metabolites = r.metabolites["comparison"]
        if gene_interest_level_r in supermodel.sources + ["assembly"]:
            r_gene_reaction_rule = r.gene_reaction_rule
        else:
            r_gene_reaction_rule = r.gene_reaction_rule["comparison"]
        out_reaction = Reaction(r.id)
        if r.name:
            out_reaction.name = r.name
        else:
            out_reaction.name = ""
        out_subsystem = ""
        for source in r.in_models["models_list"]:
            out_subsystem = (
                out_subsystem + "#" + source + "#" + r.subsystem.get(source)[0]
            )
        out_reaction.subsystem = out_subsystem
        if extend_zero_bounds and (
            r_upper_bound.get(interest_level_r)[0]
            - r_lower_bound.get(interest_level_r)[0]
            == 0
        ):
            out_reaction.lower_bound = -1000.0
            out_reaction.upper_bound = 1000.0
        else:
            out_reaction.lower_bound = r_lower_bound.get(interest_level_r)[0]
            out_reaction.upper_bound = r_upper_bound.get(interest_level_r)[0]
        for met, k in r_metabolites.get(interest_level_r).items():
            out_met = Metabolite(
                met.id, name=met.name, compartment=met.compartments["assembly"][0]
            )
            out_reaction.add_metabolites({out_met: k})
        if r_gene_reaction_rule.get(gene_interest_level_r):
            out_reaction.gene_reaction_rule = r_gene_reaction_rule.get(
                gene_interest_level_r
            )[0]
        else:
            out_reaction.gene_reaction_rule = ""
        outmodel.add_reactions([out_reaction])
    biomass_reaction = Reaction("Biomass")
    if biomass_interest_level in supermodel.sources + ["assembly"]:
        biomass_reaction.upper_bound = supermodel.reactions.assembly[
            "Biomass"
        ].upper_bound.get(biomass_interest_level)[0]
        biomass_reaction.lower_bound = supermodel.reactions.assembly[
            "Biomass"
        ].lower_bound.get(biomass_interest_level)[0]
        bio_r_metabolites = supermodel.reactions.assembly["Biomass"].metabolites
    else:
        biomass_reaction.upper_bound = (
            supermodel.reactions.assembly["Biomass"]
            .upper_bound["comparison"]
            .get(biomass_interest_level)[0]
        )
        biomass_reaction.lower_bound = (
            supermodel.reactions.assembly["Biomass"]
            .lower_bound["comparison"]
            .get(biomass_interest_level)[0]
        )
        bio_r_metabolites = supermodel.reactions.assembly["Biomass"].metabolites[
            "comparison"
        ]
    if simple_biomass_products:
        biomass_prod_id = ["adp_c", "h_c", "pi_c", "ppi_c"]
        for met, k in bio_r_metabolites.get(biomass_interest_level).items():
            if met.id in biomass_prod_id:
                biomass_prod_id.remove(met.id)
                bio_met = Metabolite(
                    met.id, name=met.name, compartment=met.compartments["assembly"][0]
                )
                biomass_reaction.add_metabolites({bio_met: k})
            if k < 0:
                bio_met = Metabolite(
                    met.id, name=met.name, compartment=met.compartments["assembly"][0]
                )
                biomass_reaction.add_metabolites({bio_met: k})
        if biomass_prod_id:
            warnings.warn(
                f"Some expected biomass products are not found: {' '.join(biomass_prod_id)}"
            )
    else:
        for met, k in bio_r_metabolites.get(biomass_interest_level).items():
            bio_met = Metabolite(
                met.id, name=met.name, compartment=met.compartments["assembly"][0]
            )
            biomass_reaction.add_metabolites({bio_met: k})
    outmodel.add_reactions([biomass_reaction])
    outmodel.objective = "Biomass"
    if gapfill_transport:
        gapfill_transport_r(outmodel, supermodel)
    if output_name is not None:
        write_sbml_model(outmodel, output_name)
        report = validate_sbml_model(filename=output_name)
        pprint(report)
    return outmodel


def get_models_with_all_confidence_levels(
    supermodel: SuperModel,
    output_dir=None,
    gene_interest_level=None,
    biomass_interest_level=None,
    simple_biomass_products=True,
    extend_zero_bounds=True,
    gapfill_transport=True,
    reactions_include: [NewElement] = None,
    reactions_exclude: [NewElement] = None,
):
    output_models = {}
    confidence_levels = deepcopy(supermodel.sources)
    for i in range(len(supermodel.sources), 1, -1):
        confidence_levels.append("core" + str(i))
    confidence_levels.append("assembly")
    for level in confidence_levels:
        if output_dir is not None:
            output_dir_lev = Path(output_dir) / (level + ".xml")
        else:
            output_dir_lev = None
        output_models.update(
            {
                level: get_model_of_interest(
                    supermodel,
                    level,
                    output_name=output_dir_lev,
                    gene_interest_level=gene_interest_level,
                    biomass_interest_level=biomass_interest_level,
                    simple_biomass_products=simple_biomass_products,
                    extend_zero_bounds=extend_zero_bounds,
                    gapfill_transport=gapfill_transport,
                    reactions_include=reactions_include,
                    reactions_exclude=reactions_exclude,
                )
            }
        )
    return output_models
