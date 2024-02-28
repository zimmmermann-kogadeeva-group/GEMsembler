import argparse
from os import chdir, getcwd
import pickle
import tempfile
from ast import literal_eval
from contextlib import contextmanager
from copy import deepcopy
from pathlib import Path
from .general import findKeysByValue
import h5py
import metquest
from cobra import Model
from cobra.io import read_sbml_model, write_sbml_model


# Context manager to get around metquest chdir internally
@contextmanager
def tmp_cwd():
    oldpwd = getcwd()
    try:
        yield
    finally:
        chdir(oldpwd)


def write_mq_data(mq_data, filename):
    with h5py.File(filename, "w") as fh:
        fh.create_dataset("all_synthesized", data=list(mq_data["all_synthesized"]))
        fh.create_dataset("max_synthesized", data=list(mq_data["max_synthesized"]))
        fh.create_dataset("max_paths_length", data=list(mq_data["max_paths_length"]))
        fh.create_dataset("tag", data=list(mq_data["tag"]))

        name_map = fh.create_group("name_map")
        for k, v in mq_data["name_map"].items():
            name_map[k] = v

        linear_paths = fh.create_group("max_linear_paths")
        for k, v in mq_data["max_linear_paths"].items():
            linear_paths[k] = str(v)

        circular_paths = fh.create_group("max_circular_paths")
        for k, v in mq_data["max_circular_paths"].items():
            circular_paths[k] = str(v)


def read_mq_data(filename, key, subkey=None, len_diversity=None, map_back=True):
    possible_keys = (
        "all_synthesized",
        "max_synthesized",
        "name_map",
        "max_linear_paths",
        "max_circular_paths",
        "max_paths_length",
        "tag"
    )
    # Raise error if not one of the possible keys
    if key not in possible_keys:
        raise ValueError(f"key has to one of the following: {possible_keys}\nGot {key}")

    with h5py.File(filename, "r") as fh:
        if key.endswith("synthesized") or (key in ["tag", "max_paths_length"]):
            return [x.decode() for x in fh[key]]
        elif key == "name_map":
            return {k: fh[key][k][()].decode() for k in fh[key].keys()}
        else:
            # If subkey is defined (for a given metabolite), then restrict the
            # output to that
            keys = fh[key].keys() if subkey is None else [subkey]
            if len_diversity is None:
                return {k: literal_eval(fh[key][k][()].decode()) for k in keys}
            else:
                output = {}
                for k in keys:
                    if k not in fh[key].keys():
                        output.update({k: None})
                    elif 0 in fh[key][k].keys():
                        output.update({k: "No need to synthesize"})
                    else:
                        len_diversity = min(len(fh[key][k].keys()), len_diversity)
                        path_lengths = list(fh[key][k].keys())[:len_diversity]
                        if map_back:
                            output.update({k: [
                                [fh["name_map"][p][()].decode() for p in path]
                                for path_length in path_lengths
                                for path in literal_eval(fh[key][k][path_length][()].decode())
                            ]})
                        else:
                            output.update({k: [
                                [p for p in path]
                                for path_length in path_lengths
                                for path in literal_eval(fh[key][k].get(path_length)[()].decode())
                            ]})
                return output


def get_all_paths(mq_data_path, metabolites, len_diversity=3):
    met_paths = {}
    comments = {}
    tag = read_mq_data(mq_data_path, "tag")[0]
    max_paths_length = read_mq_data(mq_data_path, "max_paths_length")[0]
    name_map = read_mq_data(mq_data_path, "name_map")
    all_syn = read_mq_data(mq_data_path,"all_synthesized")
    max_syn = read_mq_data(mq_data_path,"max_synthesized")

    for met in metabolites:
        # metquest adds model id to the beginning (sometimes - TODO)
        met_tag = tag + " " + met
        met_name_map = findKeysByValue(name_map, met)
        if met_name_map:
            met_name_map = met_name_map[0]
        else:
            met_name_map = ""

        # if metabolite is not in all_synthesized list, then it cannot be synthesized
        if (met_tag not in all_syn) and (met not in all_syn) and (met_name_map not in all_syn):
            comments[met] = "can not be synthesized"

        # if metabolite is not in X_synthesized then it cannot be synthesized with max path length of X
        elif (met_tag not in max_syn) and (met not in max_syn) and (met_name_map not in max_syn):
            comment = f"can not be synthesized with max {max_paths_length} length path"
            comments[met] = comment
        else:
            # Check linear paths for given metabolite
            if met_tag in max_syn:
                met_found = met_tag
            elif met in max_syn:
                met_found = met
            elif met_name_map in max_syn:
                met_found = met_name_map
            # Get linear and circular paths from metquest output for a given metabolite
            linear_paths = read_mq_data(mq_data_path, "max_linear_paths", met_found,
                                        len_diversity)
            if type(linear_paths[met_found]) == str:
                comments[met] = linear_paths[met_found]
            elif type(linear_paths[met_found]) == list:
                met_paths[met] = linear_paths[met_found]
            elif linear_paths[met_found] is None:
                # Check circular paths for given metabolite
                circular_paths = read_mq_data(mq_data_path, "max_circular_paths",
                                              met_found,
                                              len_diversity)
                if type(circular_paths[met_found]) == str:
                    comments[met] = circular_paths[met_found]
                elif type(circular_paths[met_found]) == list:
                    met_paths[met] = circular_paths[met_found]
                elif circular_paths[met_found] is None:
                    # Otherwise check ...
                    maybe_synthesised = set(all_syn) - set(max_syn)
                    if maybe_synthesised:  # if non-empty
                        comment = "Problematic: can be synthesized but does not have paths"
                    else:
                        comment = (
                            f"Maybe can not be synthesized with max {max_paths_length} length "
                            f"path, because all_synthesized not bigger than max_synthesized"
                        )
                    comments[met] = comment

    return met_paths, comments


def preprocess_medium(
    tag: str, nutritional_sources: [str], other_medium: [str], cofactors: [str],
):
    if cofactors is None:
        cofactors = [
            "co2",
            "hco3",
            "pi",
            "ppi",
            "ACP",
            "atp",
            "adp",
            "amp",
            "nad",
            "nadh",
            "nadp",
            "nadph",
            "coa",
            "cmp",
            "cdp",
            "ctp",
            "gmp",
            "gdp",
            "gtp",
            "ump",
            "udp",
            "utp",
            "fadh2",
            "fad",
            "q8",
            "q8h2",
            "mqn8",
            "mql8",
            "mqn6",
            "mql6",
            "thf",
        ]
    # TODO: Why is `_c` added to each metabolite id?
    cofactors_c = [c + "_c" for c in cofactors]
    other_medium_c = [m + "_c" for m in other_medium]
    nutritional_sources_c = [n + "_c" for n in nutritional_sources]
    seed_metabolites = nutritional_sources_c + other_medium_c + cofactors_c
    seed_metabolites_tag = set([tag + " " + s for s in seed_metabolites])
    return seed_metabolites_tag, cofactors_c, other_medium_c


def remove_cofactors_from_model(
    model: Model, medium_wo_cof_c: [str], cofactors_c: [str]
):
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


def run_metquest(
    path_to_models: str,
    seed_metabolites_tag: [str],
    number_of_xml,
    max_paths_length,
):
    graph, name_map = metquest.construct_graph.create_graph(
        path_to_models, number_of_xml
    )

    lowboundmet_pic, status_dict_pic, scope_pic = metquest.forward_pass(
        graph, seed_metabolites_tag
    )

    # Number of reactions visited in individual network for the given media (seed)
    # TODO: is this needed?
    visited_pic = {
        k: status_dict_pic[k]
        for k in list(status_dict_pic)
        if status_dict_pic[k] == "V"
    }

    # BE CAREFULLY can crash computer: memory heavy and long. Better do on cluster.
    pathway_table, c_path, scope = metquest.find_pathways(
        graph, seed_metabolites_tag, max_paths_length
    )

    output = {
        "all_synthesized": scope_pic,
        "max_paths_length": max_paths_length,
        "max_synthesized": scope,
        "max_linear_paths": pathway_table,
        "max_circular_paths": c_path,
        "name_map": name_map,
    }
    return output, graph


def pathsfinding(
    path_to_model,
    output_dir,
    nutritional_sources: [str],
    other_medium: [str],
    cofactors=None,
    max_paths_length=40,
    number_of_xml=1,
    **kwargs,
):
    # Convert paths to Path objects
    path_to_model = Path(path_to_model)
    model_name = path_to_model.stem
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Open the models
    model = read_sbml_model(path_to_model)
    tag = model.id

    # Alter metabolite ids
    seed_metabolites_tag, cofactors_c, other_medium_c = preprocess_medium(
        tag, nutritional_sources, other_medium, cofactors
    )

    # Remove cofactors from the model
    model_wo_cof = remove_cofactors_from_model(
        deepcopy(model), other_medium_c, cofactors_c
    )

    # Save the model without cofactors
    path_to_model_wo_cof = output_dir / (model_name + "_wo_cofactors.xml")
    write_sbml_model(model_wo_cof, path_to_model_wo_cof)

    # Make a tmp dir in which to run metquest
    with tempfile.TemporaryDirectory(dir=output_dir) as tmp_dir, tmp_cwd():
        # Make symbolic link to model
        tmp_model_path = Path(tmp_dir) / (model_name + ".xml")
        tmp_model_path.symlink_to(f"../{model_name}_wo_cofactors.xml")

        # Run metquest
        mq_output, graph = run_metquest(
            tmp_dir,
            seed_metabolites_tag,
            number_of_xml,
            max_paths_length,
        )

    # Add necessary parameters to output data
    mq_output.update({"tag": tag, "max_paths_length": max_paths_length})

    # Save the output
    write_mq_data(mq_output, output_dir / "metquest.h5")

    with open(output_dir / "graph.pkl", "wb") as fh:
        pickle.dump(graph, fh)
