from pathlib import Path
import tracemalloc
from copy import deepcopy
import dill
from cobra import Model
from cobra.io import read_sbml_model, write_sbml_model
import metquest


def preprocessMedium(tag: str, nutritional_sources: [str], other_medium: [str], cofactors: [str]):
    if not cofactors:
        cofactors = ["co2", "hco3", "pi", "ppi", "ACP", "atp", "adp", "amp", "nad", "nadh", "nadp", "nadph", "coa",
                     "cmp", "cdp", "ctp", "gmp", "gdp", "gtp", "ump", "udp", "utp", "fadh2", "fad", "q8", "q8h2",
                     "mqn8", "mql8", "mqn6", "mql6", "thf"]
    cofactors_c = [c + "_c" for c in cofactors]
    other_medium_c = [m + "_c" for m in other_medium]
    nutritional_sources_c = [n + "_c" for n in nutritional_sources]
    seed_metabolites = nutritional_sources_c + other_medium_c + cofactors_c
    seed_metabolites_tag = set([tag + " " + s for s in seed_metabolites])
    return seed_metabolites_tag, cofactors_c, other_medium_c


def removeCofactorsFromModel(model: Model, medium_wo_cof_c: [str], cofactors_c: [str]):
    r_to_remove = []
    met_c = medium_wo_cof_c + cofactors_c
    for r in model.reactions:
        p = True
        if (r.lower_bound == 0) and (r.upper_bound > 0):
            for react in r.reactants:
                if react.id not in met_c:
                    p = False
        elif (r.lower_bound < 0) and (r.upper_bound == 0):
            for pro in r.products:
                if pro.id not in met_c:
                    p = False
        else:
            pp = True
            for react in r.reactants:
                if react.id not in met_c:
                    pp = False
            if not pp:
                pp = True
                for pro in r.products:
                    if pro.id not in met_c:
                        pp = False
            p = pp
        if p:
            r_to_remove.append(r.id)

    model.remove_reactions(r_to_remove)
    return model


def runMetQuest(path_to_models: str, name: str, seed_metabolites_tag: [str], number_of_xml, max_paths_length):
    graph, name_map = metquest.construct_graph.create_graph(path_to_models, number_of_xml)
    lowboundmet_pic, status_dict_pic, scope_pic = metquest.forward_pass(graph, seed_metabolites_tag)
    # Number of reactions visited in individual network for the given media (seed)
    visited_pic = {k: status_dict_pic[k] for k in list(status_dict_pic) if status_dict_pic[k] == 'V'}
    print(f'Number of reactions visited in {name}: {len(visited_pic)}')
    # starting the monitoring
    tracemalloc.start()
    # function call
    # BE CAREFULLY can crash computer: memory heavy and long. Better do on cluster.
    pathway_table, c_path, scope = metquest.find_pathways(graph, seed_metabolites_tag, max_paths_length)
    # displaying the memory
    print(tracemalloc.get_traced_memory())
    # stopping the library
    tracemalloc.stop()
    output = {"all_synthesized": scope_pic, str(max_paths_length) + "_synthesized": scope,
              str(max_paths_length) + "_linear_paths": pathway_table,
              str(max_paths_length) + "_circular_paths": c_path,
              "graph": graph, "name_map": name_map}
    return output


def pathsForMetabolites(metabolites: [str], metquest_output: dict, tag: str, max_paths_length: int = 40,
                        len_diversity: int = 3):
    met_paths = {}
    for met in metabolites:
        met_paths.update({met: []})
        met_tag = tag + " " + met
        if met_tag not in metquest_output.get("all_synthesized"):
            met_paths.update({met: "can not be synthesized"})
        elif met_tag not in metquest_output.get(str(max_paths_length) + "_synthesized"):
            met_paths.update({met: f"can not be synthesized with max {max_paths_length} length path"})
        else:
            if met_tag in metquest_output.get(str(max_paths_length) + "_linear_paths").keys():
                tmp_paths = metquest_output.get(str(max_paths_length) + "_linear_paths").get(met_tag)
                if 0 in tmp_paths.keys():
                    met_paths.update({met: "No need to synthesize"})
                else:
                    for i in list(tmp_paths.keys())[:min(len(tmp_paths.keys()), len_diversity)]:
                        for path in tmp_paths.get(i):
                            met_paths.get(met).append([metquest_output.get("name_map").get(p) for p in path])
            elif met_tag in metquest_output.get(str(max_paths_length) + "_circular_paths").keys():
                tmp_paths = metquest_output.get(str(max_paths_length) + "_circular_paths").get(met_tag)
                if 0 in tmp_paths.keys():
                    met_paths.update({met: "No need to synthesize"})
                else:
                    for i in list(tmp_paths.keys())[:min(len(tmp_paths.keys()), len_diversity)]:
                        for path in tmp_paths.get(i):
                            met_paths.get(met).append([metquest_output.get("name_map").get(p) for p in path])
            else:
                if not (set(list(metquest_output.get("all_synthesized"))) - set(
                        list(metquest_output.get(str(max_paths_length) + "_synthesized")))):
                    met_paths.update({
                        met: f"Maybe can not be synthesized with max {max_paths_length} length path, because "
                             f"all_synthesized not bigger {str(max_paths_length)}_synthesized"})
                else:
                    met_paths.update({met: "Problematic: can be synthesized but does not have paths"})
    return met_paths


def runPaths(data_dir: str, model_name: str, nutritional_sources: [str], other_medium: [str], cofactors=None,
             write_metquest=True, metabolites: [str] = None, met_name: str = None, max_paths_length=40,
             len_diversity=3, number_of_xml=1):
    Path(data_dir).mkdir(parents=True, exist_ok=True)
    model = read_sbml_model(model_name + ".xml")
    tag = model.id
    seed_metabolites_tag, cofactors_c, other_medium_c = preprocessMedium(tag, nutritional_sources,
                                                                         other_medium, cofactors)
    model_wo_cof = removeCofactorsFromModel(deepcopy(model), other_medium_c, cofactors_c)
    write_sbml_model(model_wo_cof, data_dir + "/" + model_name + "_wo_cofactors.xml")
    metquest = runMetQuest(data_dir, tag, seed_metabolites_tag, number_of_xml, max_paths_length)
    if write_metquest:
        dill.dump(metquest, open(model_name + "_wo_cofactors_metquest_results.pkl", "wb"))
    if metabolites:
        met_paths = pathsForMetabolites(metabolites, metquest, tag, max_paths_length, len_diversity)
        dill.dump(met_paths, open(model_name + "_" + met_name + "_paths.pkl", "wb"))


if __name__ == '__main__':
    """user change here """
    model_name = "BU_union_model"
    dir_path = "test_MetQuest_BU_diff_models"
    nutritional_sources = ["glc__D"]
    other_medium = ["pheme", "b12", "mndn", "h2s", "k", "pi", "na1", "cl", "nh4", "so4", "mg2", "fe2", "fe3", "ca2",
                    "zn2", "mn2", "cu2", "cobalt2", "h2o", "h"]
    model = read_sbml_model(model_name + ".xml")
    union_model = read_sbml_model("BU_union_model.xml")
    met_name = "biomass_precursors"
    metabolites = [m.id for m in union_model.reactions.get_by_id("Biomass").reactants]
    find_paths = True
    get_met_from_path = False
    # end of user parameters

    if find_paths:
        runPaths(dir_path, model_name, nutritional_sources, other_medium, metabolites=metabolites, met_name=met_name, len_diversity=1)
    if get_met_from_path:
        metquest_out = dill.load(open("./" + dir_path + "/" + model_name + "_wo_cofactors_metquest_results.pkl", "rb"))
        met_paths = pathsForMetabolites(metabolites, metquest_out, model.id)
        dill.dump(met_paths, open("./" + dir_path + "/" + model_name + "_" + met_name + "_paths.pkl", "wb"))
