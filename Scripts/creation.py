import pandas as pd


class NewObject():
    """ New object class - one metabolite or reaction for supermodel. """
    def __init__(self, new_id: str, old_id: str, compartments: [str], source: str, possible_sources: [str],
                 possible_conversion: dict):
        self.id = new_id
        self.compartments = {}
        self.possible_conversion = {}
        self.sources = {}
        self.annotation = {}
        for ps in possible_sources:
            if ps == source:
                self.compartments.update({ps: compartments})
                self.possible_conversion.update({ps: [possible_conversion]})
                self.sources.update({ps: 1})
                self.annotation.update({ps: [old_id]})
            else:
                self.compartments.update({ps: []})
                self.possible_conversion.update({ps: []})
                self.sources.update({ps: 0})
                self.annotation.update({ps: []})

    def updateNewObject(self, id_to_update: str, compart_to_update: [str], conv_to_updata: dict, source: str):
        self.sources.update({source: self.sources.get(source) + 1})
        self.annotation.get(source).append(id_to_update)
        self.possible_conversion.get(source).append(conv_to_updata)
        self.compartments.update({source: self.compartments.get(source) + compart_to_update})


class SetofNewObjects():
    """ Setting dictionaries for all metabolites or reactions:
    selected for supermodel - self.converted and not selected - self.notconverted. """
    def __init__(self):
        self.converted = {}
        self.notconverted = {}

    def addNewObjs(self, selected: dict, converted: dict, where_to_add, sources):
        for source in sources:
            if source in list(selected.keys()):
                objects = selected.get(source)
                for key in objects.keys():
                    if not objects[key][1]:
                        new_id = key
                        comp = objects[key][0]
                        if source in list(converted.keys()):
                            conversion = converted.get(source).get(key)[1]
                        else:
                            conversion = {}
                        if new_id in where_to_add.keys():
                            where_to_add.get(new_id).updateNewObject(key, comp, conversion, source)
                        else:
                            new = NewObject(new_id, key, comp, source, sources, conversion)
                            where_to_add.update({new_id: new})
                    elif objects[key][1] != ["Biomass"]:
                        for new_id in objects[key][1]:
                            comp = objects[key][0]
                            if source in list(converted.keys()):
                                conversion = converted.get(source).get(key)[1]
                            else:
                                conversion = {}
                            if new_id in where_to_add.keys():
                                where_to_add.get(new_id).updateNewObject(key, comp, conversion, source)
                            else:
                                new = NewObject(new_id, key, comp, source, sources, conversion)
                                where_to_add.update({new_id: new})

    def makeSetofNew(self, selected: dict, not_selected: dict, converted: dict, sources, additional=None):
        self.addNewObjs(selected, converted, self.converted, sources)
        self.addNewObjs(not_selected, converted, self.notconverted, sources) #TODO connect not_converted for realy not converted only with old id
        if additional:
            self.addNewObjs(additional, converted, self.converted, sources)

    def makeForwardBackward(self, all_models: dict, selected: dict, obj_type: "metabolites" or "reactions",
                            additional=None):
        """ Creating dictionaries linking metabolites/reactions:
        NewObject in supermodel with old original ID and OldObject in original models with new ID in supermodel """
        goOldNew = {}
        goNewOld = {}
        types = list(selected.keys())
        for typ in types:
            goOldNew.update({typ: {}})
            for key, value in selected.get(typ).items():
                new_obj = [self.converted.get(value[1][0])]
                if additional:
                    if typ in list(additional.keys()):
                        if key in list(additional.get(typ).keys()):
                            new_obj.append(self.converted.get(additional.get(typ).get(key)[1][0]))
                goOldNew.get(typ).update({key: new_obj})
        for k, v in self.converted.items():
            goNewOld.update({k: {}})
            for t in types:
                goNewOld.get(k).update({t: {}})
                if v.annotation.get(t):
                    if obj_type == "metabolites":
                        old_obj = [all_models.get(t).metabolites.get_by_id(i) for i in v.annotation.get(t)]
                    elif obj_type == "reactions":
                        old_obj = [all_models.get(t).reactions.get_by_id(v.annotation.get(t)[0])]
                    if old_obj: goNewOld.get(k)[t] = old_obj
        return goOldNew, goNewOld


class SetofNewMetabolites(SetofNewObjects):
    """ Metabolites class that add name and blank reaction attribute to metabolite """
    def getName(self, database_info: pd.core.frame.DataFrame):
        for obj in self.converted.values():
            id_noc = obj.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
            name = database_info[database_info["universal_bigg_id"] == id_noc]["name"].values[0]
            obj.name = name
            obj.reactions = {k: [] for k in obj.sources.keys()}
            obj.formula = {k: [] for k in obj.sources.keys()}
        for ncobj in self.notconverted.values():
            name = "Not converted"
            ncobj.name = name
            ncobj.reactions = {k: [] for k in ncobj.sources.keys()}
            ncobj.formula = {k: [] for k in ncobj.sources.keys()}


class SetofNewReactions(SetofNewObjects):
    """ Reactions class that add name, reaction equation and blank reactants/products attributes to reaction """
    def getNameANDEquation(self, database_info: pd.core.frame.DataFrame):
        for obj in self.converted.values():
            id_noc = obj.id.replace("sink_", "DM_")
            name = database_info[database_info["bigg_id"] == id_noc]["name"]
            if not name.empty:
                name = name.values[0]
            equation = database_info[database_info["bigg_id"] == id_noc]["reaction_string"]
            if not equation.empty:
                equation = equation.values[0]
            else:
                equation = None
            obj.name = name
            obj.reaction = equation
            obj.reactants = {k: [] for k in obj.sources.keys()}
            obj.products = {k: [] for k in obj.sources.keys()}
            obj.metabolites = {k: [] for k in obj.sources.keys()}
            obj.lower_bound = {k: [] for k in obj.sources.keys()}
            obj.upper_bound = {k: [] for k in obj.sources.keys()}
            obj.subsystem = {k: [] for k in obj.sources.keys()}

        for ncobj in self.notconverted.values():
            name = "Not converted"
            equation = None
            ncobj.name = name
            ncobj.reaction = equation
            ncobj.reactants = {k: [] for k in ncobj.sources.keys()}
            ncobj.products = {k: [] for k in ncobj.sources.keys()}
            ncobj.metabolites = {k: [] for k in ncobj.sources.keys()}
            ncobj.lower_bound = {k: [] for k in ncobj.sources.keys()}
            ncobj.upper_bound = {k: [] for k in ncobj.sources.keys()}
            ncobj.subsystem = {k: [] for k in ncobj.sources.keys()}



class SuperModel(): #TODO add transport reactions for periplasmic metabolites for models without periplasmic compartments
    """ Supermodel class with metabolites and reactions. Sources - names of original models used to create supermodel.
    Creating connections between metabolites and reaction via dictionaries with sources as keys and links to
    reactants/products/reactions as values.  """
    def __init__(self, metabolites, reactions, types: [str]):
        self.metabolites = metabolites
        self.reactions = reactions
        self.sources = types

    def findReactions(self, metabolite: NewObject, m_goNewOld: dict, r_goOldNew: dict, types: [str],
                      periplasmic_r: dict, periplasmic_m: dict):
        for typ in types:
            old_mets = m_goNewOld.get(metabolite.id).get(typ)
            if old_mets:
                new_r = []
                for old_met in old_mets:
                    if typ in list(periplasmic_m.keys()):
                        if old_met.id in list(periplasmic_m.get(typ).keys()):
                            for reaction in old_met.reactions:
                                if r_goOldNew.get(typ).get(reaction.id):
                                    if (metabolite.id.endswith("_p")) & (
                                            reaction.id in list(periplasmic_r.get(typ).keys())):
                                        new_r.append(r_goOldNew.get(typ).get(reaction.id)[0])
                                    elif (not metabolite.id.endswith("_p")) & (
                                            reaction.id not in list(periplasmic_r.get(typ).keys())):
                                        new_r.append(r_goOldNew.get(typ).get(reaction.id)[0])
                        else:
                            for r in old_met.reactions:
                                if r_goOldNew.get(typ).get(r.id):
                                    new_r.append(r_goOldNew.get(typ).get(r.id)[0])
                    else:
                        new_r = []
                        for r in old_met.reactions:
                            if r_goOldNew.get(typ).get(r.id):
                                new_r.append(r_goOldNew.get(typ).get(r.id)[0])
                if new_r: metabolite.reactions[typ] = list(set(new_r))

    def findMetabolites(self, reaction: NewObject, r_goNewOld: dict, m_goOldNew: dict, types: [str],
                        periplasmic_r: dict):
        for typ in types:
            old_react = r_goNewOld.get(reaction.id).get(typ)
            if old_react:
                old_react_reactants = old_react[0].reactants
                old_react_products = old_react[0].products
                if typ in periplasmic_r.keys():
                    if old_react[0].id in periplasmic_r.get(typ).keys():
                        for reactant in old_react_reactants:
                            new_reactants = m_goOldNew.get(typ).get(reactant.id)
                            if new_reactants:
                                if len(new_reactants) == 1:
                                    reaction.reactants.get(typ).append(new_reactants[0])
                                elif len(new_reactants) > 1:
                                    for new_reactant in new_reactants:
                                        if (new_reactant.id.endswith("_p")) & (
                                                reactant.id in periplasmic_r.get(typ).get(old_react[0].id).keys()):
                                            reaction.reactants.get(typ).append(new_reactant)
                                        if (not new_reactant.id.endswith("_p")) & (reactant.id not in periplasmic_r.get(typ).get(old_react[0].id).keys()):
                                            reaction.reactants.get(typ).append(new_reactant)
                        for product in old_react_products:
                            new_products = m_goOldNew.get(typ).get(product.id)
                            if new_products:
                                if len(new_products) == 1:
                                    reaction.products.get(typ).append(new_products[0])
                                elif len(new_products) > 1:
                                    for new_product in new_products:
                                        if (new_product.id.endswith("_p")) & (
                                                product.id in periplasmic_r.get(typ).get(old_react[0].id).keys()):
                                            reaction.products.get(typ).append(new_product)
                                        if (not new_product.id.endswith("_p")) & (product.id not in periplasmic_r.get(typ).get(old_react[0].id).keys()):
                                            reaction.products.get(typ).append(new_product)
                    else:
                        for reactant in old_react_reactants:
                            new_reactants = m_goOldNew.get(typ).get(reactant.id)
                            if new_reactants:
                                if len(new_reactants) == 1:
                                    reaction.reactants.get(typ).append(new_reactants[0])
                                elif len(new_reactants) > 1:
                                    for new_reactant in new_reactants:
                                        if not new_reactant.id.endswith("_p"):
                                            reaction.reactants.get(typ).append(new_reactant)
                        for product in old_react_products:
                            new_products = m_goOldNew.get(typ).get(product.id)
                            if new_products:
                                if len(new_products) == 1:
                                    reaction.products.get(typ).append(new_products[0])
                                elif len(new_products) > 1:
                                    for new_product in new_products:
                                        if not new_product.id.endswith("_p"):
                                            reaction.products.get(typ).append(new_product)
                else:
                    for reactant in old_react_reactants:
                        new_reactants = m_goOldNew.get(typ).get(reactant.id)
                        if new_reactants: reaction.reactants.get(typ).append(new_reactants[0])
                    for product in old_react_products:
                        new_products = m_goOldNew.get(typ).get(product.id)
                        if new_products: reaction.products.get(typ).append(new_products[0])

    def findConnections(self, m_goNewOld: dict, m_goOldNew: dict, r_goNewOld: dict, r_goOldNew: dict, types: [str], periplasmic_r: dict, periplasmic_m: dict):
        for met in self.metabolites.converted.values():
            self.findReactions(met, m_goNewOld, r_goOldNew, types, periplasmic_r, periplasmic_m)
        for r in self.reactions.converted.values():
            self.findMetabolites(r, r_goNewOld, m_goOldNew, types, periplasmic_r)


    def addBiomass(self, types: [str], m_goOldNew: dict, all_models: dict, final_r_not_sel: dict, final_m_not_sel: dict):
        new_biomass = None
        for typ in types:
            for r in all_models.get(typ).reactions:
                if len(r.reactants) > 24:
                    if not new_biomass:
                        new_biomass = NewObject("Biomass", r.id, final_r_not_sel.get(typ).get(r.id)[0], typ, types, {})
                        new_biomass.reactants = {typ: []}
                        new_biomass.products = {typ: []}
                        nc_biomass = NewObject("Biomass", r.id, final_r_not_sel.get(typ).get(r.id)[0], typ, types, {})
                        nc_biomass.reactants = {typ: []}
                        nc_biomass.products = {typ: []}

                    else:
                        new_biomass.updateNewObject(r.id, final_r_not_sel.get(typ).get(r.id)[0], {}, typ)
                        new_biomass.reactants.update({typ: []})
                        new_biomass.products.update({typ: []})
                        nc_biomass.updateNewObject(r.id, final_r_not_sel.get(typ).get(r.id)[0], {}, typ)
                        nc_biomass.reactants.update({typ: []})
                        nc_biomass.products.update({typ: []})
                    biomass_react = [mr.id for mr in r.reactants]
                    for reactant in biomass_react:
                        new_reactants = m_goOldNew.get(typ).get(reactant)
                        if new_reactants:
                            new_biomass.reactants.get(typ).append(new_reactants[0])
                            new_reactants[0].reactions.get(typ).append(new_biomass)
                        else:
                            if not final_m_not_sel.get(typ).get(reactant)[1]:
                                nc_biomass.reactants.get(typ).append(self.metabolites.notconverted.get(reactant))
                                self.metabolites.notconverted.get(reactant).reactions.get(typ).append(nc_biomass)
                            else:
                                for bigg_reactant in final_m_not_sel.get(typ).get(reactant)[1]:
                                    nc_biomass.reactants.get(typ).append(self.metabolites.notconverted.get(bigg_reactant))
                                    self.metabolites.notconverted.get(bigg_reactant).reactions.get(typ).append(nc_biomass)
                    biomass_pro = [mp.id for mp in r.products]
                    for product in biomass_pro:
                        new_products = m_goOldNew.get(typ).get(product)
                        if new_products:
                            new_biomass.products.get(typ).append(new_products[0])
                            new_products[0].reactions.get(typ).append(new_biomass)
                        else:
                            if not final_m_not_sel.get(typ).get(product)[1]:
                                nc_biomass.products.get(typ).append(self.metabolites.notconverted.get(product))
                                self.metabolites.notconverted.get(product).reactions.get(typ).append(nc_biomass)
                            else:
                                for bigg_product in final_m_not_sel.get(typ).get(product)[1]:
                                    nc_biomass.products.get(typ).append(self.metabolites.notconverted.get(bigg_product))
                                    self.metabolites.notconverted.get(bigg_product).reactions.get(typ).append(nc_biomass)
        self.reactions.converted.update({"Biomass": new_biomass})
        self.reactions.notconverted.update({"Biomass": nc_biomass})


def runSupermodelCreation(model_type, final_m, final_m_not_sel, final_r, final_r_not_sel, all_models, bigg_all_m,
                          bigg_all_r, additional_periplasmic_m, periplasmic_r):
    """ Creating supermodel with metabolites and reactions. """
    metabolites = SetofNewMetabolites()
    metabolites.makeSetofNew(final_m, final_m_not_sel, bigg_all_m, model_type, additional_periplasmic_m)
    metabolites.getName(bigg_all_m)
    reactions = SetofNewReactions()
    reactions.makeSetofNew(final_r, final_r_not_sel, bigg_all_r, model_type)
    reactions.getNameANDEquation(bigg_all_r)
    m_goOldNew, m_goNewOld = metabolites.makeForwardBackward(all_models, final_m, "metabolites",
                                                             additional_periplasmic_m)
    r_goOldNew, r_goNewOld = reactions.makeForwardBackward(all_models, final_r, "reactions")
    supermodel = SuperModel(metabolites, reactions, model_type)
    supermodel.findConnections(m_goNewOld, m_goOldNew, r_goNewOld, r_goOldNew, model_type, periplasmic_r, additional_periplasmic_m)
    supermodel.addBiomass(model_type, m_goOldNew, all_models, final_r_not_sel, final_m_not_sel)
    return supermodel