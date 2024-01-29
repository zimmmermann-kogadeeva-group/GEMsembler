from collections import defaultdict
import pandas as pd
from cobra.flux_analysis import pfba


def run_growth_full_flux_analysis(
    models_to_analyse: dict,
    medium: list,
    biomass_r_id=None,
    m_compartment_suffix=None,
    bp_compartment_suffix=None,
):
    if biomass_r_id is None:
        biomass_r_id = "Biomass"
    if m_compartment_suffix is None:
        m_compartment_suffix = "_e"
    if bp_compartment_suffix is None:
        bp_compartment_suffix = "_c"
    biomass_p = []
    for m in models_to_analyse.values():
        biomass_p = list(
            set(
                biomass_p
                + [react.id for react in m.reactions.get_by_id(biomass_r_id).reactants]
            )
        )
    out_bp_production = {"biomass_precursors": ["overall_growth"] + biomass_p}
    flux_res_out = defaultdict(dict)
    path_pfba_out = defaultdict(dict)
    for k, model in models_to_analyse.items():
        med_mod = {}
        for e in model.exchanges:
            if (
                list(e.metabolites.keys())[0].id.removesuffix(m_compartment_suffix)
                in medium
            ):
                med_mod.update({e.id: 1000})
            else:
                med_mod.update({e.id: 0})
        old_medium = model.medium
        model.medium = med_mod
        res = model.optimize()
        flux_res_out[k] = {"overall_fba": res}
        model_data = [res.objective_value]
        all_met = [mm.id for mm in model.metabolites]
        all_r = [rr.id for rr in model.reactions]
        for bp in biomass_p:
            if bp not in all_met:
                model_data.append("nan")
            else:
                demand_added = False
                if ("DM_" + bp) not in all_r:
                    model.add_boundary(model.metabolites.get_by_id(bp), type="demand")
                    demand_added = True
                model.objective = model.demands.get_by_id("DM_" + bp)
                res_bp = model.optimize()
                flux_res_out[k].update(
                    {bp.removesuffix(bp_compartment_suffix) + "_fba": res_bp}
                )
                if res_bp.objective_value > 0:
                    model_data.append(1)
                    pfba_res = pfba(model)
                    flux_res_out[k].update(
                        {bp.removesuffix(bp_compartment_suffix) + "_pfba": pfba_res}
                    )
                    reactions = list(
                        pfba_res.to_frame()[
                            (pfba_res.to_frame()["fluxes"] > 0.001)
                            | (pfba_res.to_frame()["fluxes"] < -0.001)
                        ].index
                    )
                    path_pfba_out[k].update(
                        {
                            bp.removesuffix(bp_compartment_suffix)
                            + "_path_pfba": reactions
                        }
                    )
                else:
                    model_data.append(0)
                if demand_added:
                    model.remove_reactions(["DM_" + bp])
        model.medium = old_medium
        out_bp_production.update({k: model_data})
    out_bp_production_tab = pd.DataFrame(out_bp_production)
    return out_bp_production_tab, flux_res_out, path_pfba_out
