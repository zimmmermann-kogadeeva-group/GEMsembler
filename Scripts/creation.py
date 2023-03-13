import pandas as pd


class NewObject():
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
    def __init__(self):
        self.converted = {}
        self.notconverted = {}

    def addNewObjs(self, selected: dict, converted: dict, where_to_add, sources):
        for source in sources:
            if source in list(selected.keys()):
                objects = selected.get(source)
                for key in objects.keys():
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
        self.addNewObjs(not_selected, converted, self.notconverted, sources)
        if additional:
            self.addNewObjs(additional, converted, self.converted, sources)

    def makeForwardBackward(self, all_models: dict, selected: dict, obj_type: "metabolites" or "reactions",
                            additional=None):
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

    def getName(self, database_info: pd.core.frame.DataFrame):
        for obj in self.converted.values():
            id_noc = obj.id.removesuffix("_c").removesuffix("_e").removesuffix("_p")
            name = database_info[database_info["universal_bigg_id"] == id_noc]["name"].values[0]
            obj.name = name


class SetofNewReactions(SetofNewObjects):
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


class SuperModel():
    def __init__(self, metabolites, reactions, types: [str]):
        self.metabolites = metabolites
        self.reactions = reactions
        self.sources = types

    def findReactions(self, metabolite: NewObject, m_goNewOld: dict, r_goOldNew: dict, types: [str],
                      periplasmic_r: dict, periplasmic_m: dict):
        metabolite.reactions = {}
        for typ in types:
            metabolite.reactions.update({typ: []})
            old_mets = m_goNewOld.get(metabolite.id).get(typ)
            if old_mets:
                for old_met in old_mets:
                    # if metabolite.id == "h_p":
                        # print(old_met.id)
                    if typ in list(periplasmic_m.keys()):
                        if old_met.id in list(periplasmic_m.get(typ).keys()):
                            new_r = []
                            # if metabolite.id == "h_p":
                            #     print(old_met.id)
                            for reaction in old_met.reactions:
                                # if metabolite.id == "h_p":
                                #     if reaction.id in ['rxn09163_c0', 'rxn08815_c0', 'rxn09149_c0', 'rxn12549_c0']:
                                #         print(reaction.id)
                                #         print(r_goOldNew.get(typ).get(reaction.id)[0].id)
                                if r_goOldNew.get(typ).get(reaction.id):
                                    if (metabolite.id.endswith("_p")) & (
                                            reaction.id in list(periplasmic_r.get(typ).keys())):
                                        # if metabolite.id == "h_p":
                                        #     if reaction.id in ['rxn09163_c0', 'rxn08815_c0', 'rxn09149_c0',
                                        #                        'rxn12549_c0']:
                                        #         print(reaction.id)
                                        #         print("enter_p")
                                        new_r.append(r_goOldNew.get(typ).get(reaction.id)[0])
                                        # if metabolite.id == "h_p":
                                        #     print("enter")
                                        #     print(reaction.id)
                                        #     print(r_goOldNew.get(typ).get(reaction.id)[0].id)
                                        #     print(len(new_r))
                                        #     print(new_r)
                                    elif (not metabolite.id.endswith("_p")) & (
                                            reaction.id not in list(periplasmic_r.get(typ).keys())):
                                        new_r.append(r_goOldNew.get(typ).get(reaction.id)[0])
                        else:
                            new_r = []
                            for r in old_met.reactions:
                                if r_goOldNew.get(typ).get(r.id):
                                    new_r.append(r_goOldNew.get(typ).get(r.id)[0])
                    else:
                        new_r = []
                        for r in old_met.reactions:
                            if r_goOldNew.get(typ).get(r.id):
                                new_r.append(r_goOldNew.get(typ).get(r.id)[0])
                    # if metabolite.id == "h_p":
                    #     print(new_r)
                    #     print(len(new_r))
                    if new_r: metabolite.reactions[typ] = new_r

    def findMetabolites(self, reaction: NewObject, r_goNewOld: dict, m_goOldNew: dict, types: [str],
                        periplasmic_r: dict):
        reaction.reactants = {}
        reaction.products = {}
        for typ in types:
            reaction.reactants.update({typ: []})
            reaction.products.update({typ: []})
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


def runSupermodelCreation(model_type, final_m, final_m_not_sel, final_r, final_r_not_sel, all_models, bigg_all_m,
                          bigg_all_r, additional_periplasmic_m, periplasmic_r):
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
    return supermodel