from copy import deepcopy

import cobra
import pandas as pd


def remove_b_type_exchange(
    model: cobra.core.model.Model,
) -> cobra.core.model.Model:
    """
    Removes boundary metabolites / exchange reactions with _b at the end if
    another (_e) compartment exist. Replacing _b with normal compartment if id
    with normal compartment doesn't exist separately.

    Parameters
    ----------
    model : cobra.core.model.Model
        model from which metabolites and reactions in `b` comparment will be
        removed. Input model will be unchanged.

    Returns
    -------
    cobra.core.model.Model
        Model without metabolites/reactions in `b` compartments
    """
    model = deepcopy(model)

    # Rename metabolites ending with "_b" and if new id of a given metabolite
    # is already in the list of metabolites, then remove it
    met_to_remove = []
    for met in model.metabolites:
        if met.id.endswith("_b"):
            new_id = met.id.removesuffix("_b") + "_" + met.compartment
            if new_id in model.metabolites:
                met_to_remove.append(met)
            else:
                met.id = new_id

    # Rename reactions ending with "_b" and if new id of a given reaction
    # is already in the list of reactions, then remove it
    r_to_remove = []
    for r in model.reactions:
        if r.id.endswith("_b"):
            new_id = r.id.removesuffix("_b") + "_" + list(r.compartments)[0]
            if new_id in model.reactions:
                r_to_remove.append(r)
            else:
                r.id = new_id

    # Remove appropriate metabolites and reactions
    model.remove_metabolites(met_to_remove)
    model.remove_reactions(r_to_remove)

    # Repair the model
    model.repair()

    return model


def get_rows(reaction: cobra.core.Reaction, do_not_care_about_directions=True):
    """
    Returns a list of row(s) to build a table of reaction information (id,
    products, reactants, and gene reaction rule). Two reactions are returned
    when flux can be positive or negative i.e. flux lower bound is negative and
    upper bound is positive.
    """

    # Get all reactants and products and convert the resulting list to a string
    reactants = " ".join(sorted([x.id for x in reaction.reactants]))
    products = " ".join(sorted([x.id for x in reaction.products]))

    # Temporary solution not directions into account, so it is the same as structural
    # TODO: change to option bellow with directionality in bigg network dict
    if do_not_care_about_directions:
        return [
            [reaction.id, reactants, products, reaction.gene_reaction_rule],
            [reaction.id, products, reactants, reaction.gene_reaction_rule],
        ]

    # Depending on the flux lower and upper bounds return products and reactants
    # in specific order.
    # TODO: what about reactions that have lower_bound and upper_bound both equal to
    # zero like reaction 1368 in ./example/BU_gapseq.xml.gz model
    if (reaction.lower_bound == 0) and (reaction.upper_bound >= 0):
        return [[reaction.id, reactants, products, reaction.gene_reaction_rule]]
    elif (reaction.lower_bound < 0) and (reaction.upper_bound == 0):
        return [[reaction.id, products, reactants, reaction.gene_reaction_rule]]
    elif (reaction.lower_bound < 0) and (reaction.upper_bound > 0):
        return [
            [reaction.id, reactants, products, reaction.gene_reaction_rule],
            [reaction.id, products, reactants, reaction.gene_reaction_rule],
        ]


def get_duplicated_reactions(model):
    """
    Parameters
    ----------
    model : cobra.core.model.Model
        a cobra model

    Returns
    -------
    DataFrame : pandas.DataFrame
        table of duplicated reactions based on reactants and products
    """
    # Create a table of information on reactions. Double list comprehension is
    # needed as sometimes two reactions are returned by `get_rows()` function
    data = pd.DataFrame(
        [row for reaction in model.reactions for row in get_rows(reaction)],
        columns=["ID", "Reactants", "Products", "GPR"],
    )

    return (
        data[data.duplicated(subset=["Reactants", "Products"], keep=False)]
        .sort_values(["Reactants", "Products"])
        .reset_index(drop=True)
    )
