import operator
import re
import warnings
from collections import defaultdict
import pandas as pd
from cobra.flux_analysis import pfba
from .creation import SuperModel
from .comparison import getCoreGPR
import seaborn as sns
from .drawing import draw_one_known_pathway, draw_pfba_results


def write_biomass_precursors_fba_production(
    out_bp_production_tab,
    write_output_to_folder: str,
    yticklabels=True,
    cmap="mako",
    metric="jaccard",
    method="single",
    center=0.3,
    linewidths=0.003,
    **kwargs,
):
    """Function to plot heatmap of produced or not produced metabolites.
    Parameters: table with production, output fold and
    attributes of seaborn clustermap function."""
    out_bp_production_tab.to_csv(
        f"{write_output_to_folder}/all_biomass_precursors_production.tsv",
        sep="\t",
        index=False,
    )
    out_core_bp_production_tab = out_bp_production_tab[
        (out_bp_production_tab != "nan").all(axis=1)
    ]
    out_core_bp_production_tab.index = out_core_bp_production_tab["biomass_precursors"]
    num = out_core_bp_production_tab[
        out_core_bp_production_tab.columns.difference(["biomass_precursors"])
    ]
    num_no_grow = num.drop(index="overall_growth")
    num_no_grow = num_no_grow.astype(
        {mod_name: "float" for mod_name in num_no_grow.columns}
    )
    num_no_grow.to_csv(
        f"{write_output_to_folder}/core_biomass_precursors_production.tsv", sep="\t",
    )
    bp_heatmap = sns.clustermap(
        num_no_grow,
        yticklabels=yticklabels,
        metric=metric,
        method=method,
        cmap=cmap,
        center=center,
        linewidths=linewidths,
        **kwargs,
    )
    bp_heatmap.savefig(
        f"{write_output_to_folder}/core_biomass_precursors_production.png"
    )


def table_reactions_confidence(
    supermodel: SuperModel,
    output_name: str,
    pathway_r=None,
    yes_range=1,
    no_range=1,
    genes=True,
    and_as_solid=False,
):
    output = {"Reaction": pathway_r, "R_confidence": [], "Status": [], "R_models": []}
    if genes:
        output.update({"GPR_confidence": [], "GPR_core": [], "GPR_assembly": []})
    if pathway_r is None:
        pathway_r = list(supermodel.reactions.assembly.keys())
    for r_id in pathway_r:
        if r_id in supermodel.reactions.assembly.keys():
            r = supermodel.reactions.assembly[r_id]
            output["R_confidence"].append(f"Core{r.in_models['models_amount']}")
            if r.in_models["models_amount"] >= len(supermodel.sources) - yes_range:
                output["Status"].append("yes")
            elif r.in_models["models_amount"] <= no_range:
                output["Status"].append("no")
            else:
                output["Status"].append("q")
            output["R_models"].append(" ".join(sorted(r.in_models["models_list"])))
            if genes:
                if not r.gene_reaction_rule["assembly"]:
                    output["GPR_confidence"].append("Core0")
                    output["GPR_core"].append("")
                    output["GPR_assembly"].append("")
                else:
                    for i in range(r.in_models["models_amount"], 0, -1):
                        gpr_core = getCoreGPR(
                            r.gene_reaction_rule,
                            i,
                            operator.ge,
                            r.in_models["models_list"],
                            and_as_solid,
                        )
                        if gpr_core:
                            output["GPR_confidence"].append(f"Core{i}")
                            output["GPR_core"].append(gpr_core[0])
                            output["GPR_assembly"].append(
                                r.gene_reaction_rule["assembly"][0]
                            )
                            break
        else:
            output["R_confidence"].append("Core0")
            output["Status"].append("no")
            output["R_models"].append("")
            if genes:
                output["GPR_confidence"].append("Core0")
                output["GPR_core"].append("")
                output["GPR_assembly"].append("")
    output_tb = pd.DataFrame(output)
    output_tb.to_csv(output_name, sep="\t", index=False)
    return output_tb


def draw_glycolysis(
    supermodel: SuperModel,
    output_name: str,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=None,
    hei=None,
    write_table=True,
    yes_range=1,
    no_range=1,
):
    glycolysis = {
        "HEX1": [("glc__D_c", "g6p_c")],
        "PGI": [("g6p_c", "f6p_c")],
        "PFK": [("f6p_c", "fdp_c")],
        "FBA": [("fdp_c", "g3p_c")],
        "GAPD": [("g3p_c", "13dpg_c")],
        "PGK": [("13dpg_c", "3pg_c")],
        "PGM": [("3pg_c", "2pg_c")],
        "ENO": [("2pg_c", "pep_c")],
        "PYK": [("pep_c", "pyr_c")],
    }
    g = draw_one_known_pathway(
        supermodel,
        glycolysis,
        output_name,
        additional_met,
        genes,
        and_as_solid,
        directed,
        met_not_int,
        n_letter,
        wid,
        hei,
    )
    if write_table:
        t_name = re.sub(".html$", ".tsv", output_name)
        t = table_reactions_confidence(
            supermodel,
            list(glycolysis.keys()),
            t_name,
            yes_range,
            no_range,
            genes,
            and_as_solid,
        )
        return g, t
    else:
        return g


def draw_tca(
    supermodel: SuperModel,
    output_name: str,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=None,
    hei=None,
    write_table=True,
    yes_range=1,
    no_range=1,
):
    tca = {
        "PYK": [("pep_c", "pyr_c")],
        "PPC": [("pep_c", "oaa_c")],
        "PPCK": [("pep_c", "oaa_c")],
        "PEPCK_re": [("pep_c", "oaa_c")],
        "PC": [("pyr_c", "oaa_c")],
        "PDH": [("pyr_c", "accoa_c")],
        "PFL": [("pyr_c", "accoa_c")],
        "CS": [("accoa_c", "cit_c"), ("oaa_c", "cit_c")],
        "ACONT": [("cit_c", "icit_c")],
        "ACONTa": [("cit_c", "acon_C_c")],
        "ACONTb": [("acon_C_c", "icit_c")],
        "ICL": [("icit_c", "succ_c"), ("icit_c", "glx_c")],
        "ICDHyr": [("icit_c", "akg_c")],
        "ICDHx": [("icit_c", "akg_c")],
        "ICITRED": [("icit_c", "osuc_c")],
        "OSUCCL": [("osuc_c", "akg_c")],
        "AKGDH": [("akg_c", "succoa_c")],
        "OOR2r": [("akg_c", "succoa_c"), ("fdxo_42_c", "fdxr_42_c")],
        "AKGDa": [("akg_c", "sdhlam_c"), ("lpam_c", "sdhlam_c")],
        "AKGDb": [("sdhlam_c", "succoa_c"), ("sdhlam_c", "dhlam_c")],
        "PDHcr": [("dhlam_c", "lpam_c")],
        "SUCOAS": [("succoa_c", "succ_c")],
        "SUCDi": [("succ_c", "fum_c")],
        "FRD7": [("fum_c", "succ_c")],
        "FUM": [("fum_c", "mal__L_c")],
        "MALS": [("glx_c", "mal__L_c")],
        "MDH": [("mal__L_c", "oaa_c")],
        "MDH2": [("mal__L_c", "oaa_c")],
        "MDH3": [("mal__L_c", "oaa_c")],
    }
    g = draw_one_known_pathway(
        supermodel,
        tca,
        output_name,
        additional_met,
        genes,
        and_as_solid,
        directed,
        met_not_int,
        n_letter,
        wid,
        hei,
    )
    if write_table:
        t_name = re.sub(".html$", ".tsv", output_name)
        t = table_reactions_confidence(
            supermodel,
            list(tca.keys()),
            t_name,
            yes_range,
            no_range,
            genes,
            and_as_solid,
        )
        return g, t
    else:
        return g


def draw_pentose_phosphate(
    supermodel: SuperModel,
    output_name: str,
    additional_met=False,
    genes=True,
    and_as_solid=False,
    directed=True,
    met_not_int=None,
    n_letter=None,
    wid=None,
    hei=None,
    write_table=True,
    yes_range=1,
    no_range=1,
):
    pentose_phosphate_pathway = {
        "G6PBDH": [("g6p_B_c", "6pgl_c")],
        "G6PDH2r": [("g6p_c", "6pgl_c")],
        "PGL": [("6pgl_c", "6pgc_c")],
        "GND": [("6pgc_c", "ru5p__D_c")],
        "RPE": [("ru5p__D_c", "xu5p__D_c")],
        "RPI": [("ru5p__D_c", "r5p_c")],
        "TKT1": [("xu5p__D_c", "g3p_c"), ("r5p_c", "s7p_c")],
        "TALA": [("g3p_c", "f6p_c"), ("s7p_c", "e4p_c")],
        "TKT2": [("e4p_c", "f6p_c"), ("xu5p__D_c", "g3p_c")],
        "PGI": [("f6p_c", "g6p_c")],
        "G6PI": [("g6p_c", "g6p_B_c")],
    }
    g = draw_one_known_pathway(
        supermodel,
        pentose_phosphate_pathway,
        output_name,
        additional_met,
        genes,
        and_as_solid,
        directed,
        met_not_int,
        n_letter,
        wid,
        hei,
    )
    if write_table:
        t_name = re.sub(".html$", ".tsv", output_name)
        t = table_reactions_confidence(
            supermodel,
            list(pentose_phosphate_pathway.keys()),
            t_name,
            yes_range,
            no_range,
            genes,
            and_as_solid,
        )
        return g, t
    else:
        return g


def run_growth_full_flux_analysis(
    models_to_analyse: dict,
    medium: list,
    write_output_to_folder=None,
    draw_pfba_for_models=None,
    supermodel_for_draw_pfba=None,
    draw_met_not_int=False,
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
        bp_model = [
            react.id for react in model.reactions.get_by_id(biomass_r_id).reactants
        ]
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
                if res_bp.objective_value > 0.001:
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
                    if bp in bp_model:
                        model_data.append(0)
                    else:
                        model_data.append(0.5)
                if demand_added:
                    model.remove_reactions(["DM_" + bp])
        model.medium = old_medium
        out_bp_production.update({k: model_data})
    out_bp_production_tab = pd.DataFrame(out_bp_production)
    if write_output_to_folder is not None:
        write_biomass_precursors_fba_production(
            out_bp_production_tab, write_output_to_folder
        )
    if supermodel_for_draw_pfba is not None:
        if write_output_to_folder is None:
            raise ValueError(
                "Output folder is not provided. "
                "Please provide output folder if you want to draw pfba pathways."
            )
        if draw_pfba_for_models is None:
            warnings.warn(
                "Models for which to draw pfba of biomass precursors is not provided. "
                "So for each biomass precursor the model with "
                "the highest confidence level will be used."
            )
        draw_pfba_results(
            path_pfba_out,
            supermodel_for_draw_pfba,
            medium,
            write_output_to_folder,
            draw_pfba_for_models,
            draw_met_not_int,
        )

    return out_bp_production_tab, flux_res_out, path_pfba_out
