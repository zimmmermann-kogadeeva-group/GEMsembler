import pandas as pd


class NewObject():
    def __init__(self, new_id: str, old_id: str, compartments: [str], source: str, possible_sources: [str], possible_conversion: dict):
        self.id = new_id
        self.compartments = compartments
        self.possible_conversion = possible_conversion
        self.sources = {}
        self.annotation = {}
        for ps in possible_sources:
            if ps == source:
                self.sources.update({ps: 1})
                self.annotation.update({ps: [old_id]})
            else:
                self.sources.update({ps: 0})
                self.annotation.update({ps: []})

    def updateNewObject(self, id_to_update: str, source: str):
        self.sources.update({source: self.sources.get(source) + 1})
        self.annotation.get(source).append(id_to_update)


class SetofNewObjects():
    def __init__(self):
        self.converted = {}
        self.notconverted = {}

    def addNewObjs(self, selected: dict, converted: dict, where_to_add):
        sources = list(selected.keys())
        for source in sources:
            objects = selected.get(source)
            for key in objects.keys():
                for new_id in objects[key][1]:
                    if new_id in where_to_add.keys():
                        where_to_add.get(new_id).updateNewObject(key, source)
                    else:
                        comp = objects[key][0]
                        if source in list(converted.keys()):
                            conversion = converted.get(source).get(key)[1]
                        else:
                            conversion = {}
                        new = NewObject(new_id, key, comp, source, sources, conversion)
                        where_to_add.update({new_id: new})

    def makeSetofNew(self, selected: dict, not_selected: dict, converted: dict):
        self.addNewObjs(selected, converted, self.converted)
        self.addNewObjs(not_selected, converted, self.notconverted)

    def makeForwardBackward(self, all_models: dict, selected: dict, obj_type: "metabolites" or "reactions"):
        goOldNew = {}
        goNewOld = {}
        types = list(selected.keys())
        for typ in types:
            goOldNew.update({typ: {}})
            for key, value in selected.get(typ).items():
                if obj_type == "metabolites":
                    old_obj = all_models.get(typ).metabolites.get_by_id(key)
                elif obj_type == "reactions":
                    old_obj = all_models.get(typ).reactions.get_by_id(key)
                new_obj = self.converted.get(value[1][0])
                goOldNew.get(typ).update({key: new_obj})
                if goNewOld.get(value[1][0]) == None:
                    goNewOld.update({value[1][0]: {typ: old_obj}})
                else:
                    goNewOld.get(value[1][0]).update({typ: old_obj})
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

    def findReactions(self, metabolite: NewObject, m_goNewOld: dict, r_goOldNew: dict, types: [str]):
        metabolite.reactions = {}
        for typ in types:
            metabolite.reactions.update({typ: []})
            old_met = m_goNewOld.get(metabolite.id).get(typ)
            if old_met is not None:
                old_met_react = old_met.reactions
                for r in old_met_react:
                    new_r = r_goOldNew.get(typ).get(r.id)
                    if new_r is not None:
                        metabolite.reactions.get(typ).append(new_r)

    def findMetabolites(self, reaction: NewObject, r_goNewOld: dict, m_goOldNew: dict, types: [str]):
        reaction.reactants = {}
        reaction.products = {}
        for typ in types:
            reaction.reactants.update({typ: []})
            reaction.products.update({typ: []})
            old_react = r_goNewOld.get(reaction.id).get(typ)
            if old_react is not None:
                old_react_reactants = old_react.reactants
                for reactant in old_react_reactants:
                    new_reactant = m_goOldNew.get(typ).get(reactant.id)
                    if new_reactant is not None:
                        reaction.reactants.get(typ).append(new_reactant)
                old_react_products = old_react.products
                for product in old_react_products:
                    new_product = m_goOldNew.get(typ).get(product.id)
                    if new_product is not None:
                        reaction.products.get(typ).append(new_product)

    def findConnections(self, m_goNewOld: dict, m_goOldNew: dict, r_goNewOld: dict, r_goOldNew: dict, types: [str]):
        for met in self.metabolites.converted.values():
            self.findReactions(met, m_goNewOld, r_goOldNew, types)
        for r in self.reactions.converted.values():
            self.findMetabolites(r, r_goNewOld, m_goOldNew, types)



def runSupermodelCreation(model_type, final_m, final_m_not_sel, final_r, final_r_not_sel, all_models, bigg_all_m, bigg_all_r):
    metabolites = SetofNewMetabolites()
    metabolites.makeSetofNew(final_m, final_m_not_sel, bigg_all_m)
    metabolites.getName(bigg_all_m)
    reactions = SetofNewReactions()
    reactions.makeSetofNew(final_r, final_r_not_sel, bigg_all_r)
    reactions.getNameANDEquation(bigg_all_r)
    m_goOldNew, m_goNewOld = metabolites.makeForwardBackward(all_models, final_m, "metabolites")
    r_goOldNew, r_goNewOld = reactions.makeForwardBackward(all_models, final_r, "reactions")
    supermodel = SuperModel(metabolites, reactions, model_type)
    supermodel.findConnections(m_goNewOld, m_goOldNew, r_goNewOld, r_goOldNew, model_type)
    return supermodel