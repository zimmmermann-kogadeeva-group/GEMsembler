import cobra

from .structural import getReaction


class PeriplasmicM(object):
    def __init__(self, bigg_o: str, bigg_p: str, tr_r_c: dict, tr_r_e: dict):
        self.bigg_orig_comp = bigg_o
        self.bigg_p = bigg_p
        if tr_r_c:
            self.transport_r_c = list(tr_r_c.keys())[0]
        else:
            self.transport_r_c = "not_exist_transport_p_c"
        if tr_r_e:
            self.transport_r_e = list(tr_r_e.keys())[0]
        else:
            self.transport_r_e = "not_exist_transport_p_e"
        self.periplasmic_r = 1
        self.all_r = None
        self.replace = None

    def add_periplasmic_r(self):
        self.periplasmic_r = self.periplasmic_r + 1

    def replace_status(self, all_r_for_met: int):
        self.all_r = all_r_for_met
        if all_r_for_met > self.periplasmic_r:
            self.replace = False
        else:
            self.replace = True


class PeriplasmicR(object):
    def __init__(self, bigg_r: str):
        self.bigg_r_p = bigg_r
        self.metabolites_changed_p = {}

    def add_changed_met(self, orig_met: str, bigg_p: str):
        self.metabolites_changed_p.update({orig_met: bigg_p})


def getSuggestionPeriplasmic(
    structural_r_sel: dict,
    structural_r_second: dict,
    met_sel: dict,
    model: cobra.core.model.Model,
    bigg_network: dict,
):
    """ Getting dictionary of metabolites that gave reaction equation if their compartment is changed to periplasmic and
     dictionary with corresponding reactions. Also, getting ids for transport reactions for metabolite from original
     compartment to periplasmic. Checking which metabolites are changed to periplasmic in all their reactions (replace)
     and wich in part of their reactions (not_replace - split). """
    met_periplasmic = {}
    react_periplasmic = {}
    bigg_met_comp_sel = [
        sel.highest_consistent[0]
        for sel in met_sel.values()
        if (sel.to_one_id == True and sel.from_one_id == True)
    ]
    for orig_id, struct in structural_r_second.items():
        if (
            structural_r_sel[orig_id].to_one_id == True
            and structural_r_sel[orig_id].from_one_id == True
            and struct.comment.startswith("Found_via_adding_periplasmic_compartment")
        ):
            react_periplasmic.update({orig_id: PeriplasmicR(struct.structural[0])})
            for i in range(len(struct.suggestions["orig_m"])):
                react_periplasmic.get(orig_id).add_changed_met(
                    struct.suggestions["orig_m"][i], struct.suggestions["p_m"][i]
                )
                if struct.suggestions["orig_m"][i] not in met_periplasmic.keys():
                    # transport connection with cellular metabolite
                    if struct.suggestions["p_m"][i][:-1] + "c" in bigg_met_comp_sel:
                        tr_p_c = getReaction(
                            [struct.suggestions["p_m"][i][:-1] + "c"],
                            [struct.suggestions["p_m"][i]],
                            bigg_network,
                            "transport_r_for_periplasmic_c",
                        )
                    else:
                        tr_p_c = {"No_c_met": "no_c_transport_possible"}
                    # transport connection with extracellular metabolite
                    if struct.suggestions["p_m"][i][:-1] + "e" in bigg_met_comp_sel:
                        tr_p_e = getReaction(
                            [struct.suggestions["p_m"][i][:-1] + "e"],
                            [struct.suggestions["p_m"][i]],
                            bigg_network,
                            "transport_r_for_periplasmic_e",
                        )
                    else:
                        tr_p_e = {"No_e_met": "no_e_transport_possible"}
                    met_periplasmic.update(
                        {
                            struct.suggestions["orig_m"][i]: PeriplasmicM(
                                struct.suggestions["b_m"][i],
                                struct.suggestions["p_m"][i],
                                tr_p_c,
                                tr_p_e,
                            )
                        }
                    )
                else:
                    met_periplasmic.get(
                        struct.suggestions["orig_m"][i]
                    ).add_periplasmic_r()
    for orig_m, periplasmic_m in met_periplasmic.items():
        periplasmic_m.replace_status(len(model.metabolites.get_by_id(orig_m).reactions))
    return met_periplasmic, react_periplasmic
