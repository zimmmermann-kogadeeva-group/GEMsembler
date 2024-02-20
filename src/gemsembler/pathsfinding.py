import argparse
from os import chdir, getcwd
import pickle
import tempfile
from ast import literal_eval
from contextlib import contextmanager
from copy import deepcopy
from pathlib import Path

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

        name_map = fh.create_group("name_map")
        for k, v in mq_data["name_map"].items():
            name_map[k] = v

        linear_paths = fh.create_group("max_linear_paths")
        for k, v in mq_data["max_linear_paths"].items():
            linear_paths[k] = str(v)

        circular_paths = fh.create_group("max_circular_paths")
        for k, v in mq_data["max_circular_paths"].items():
            circular_paths[k] = str(v)


def read_mq_data(filename, key, subkey=None):
    with h5py.File(filename, "r") as fh:
        if key.endswith("synthesized"):
            return [x.decode() for x in fh[key]]
        else:
            # If subkey is defined (for a given metabolite), then restrict the
            # output to that
            keys = fh[key].keys() if subkey is None else [subkey]
            if key == "name_map":
                return {k: fh[key][k][()].decode() for k in keys}
            elif key.endswith("paths"):
                return {k: literal_eval(fh[key][k][()].decode()) for k in keys}


def get_paths(subset, name_map, len_diversity=3):
    # Update len_diversity - len_diversity cannot be greater than number of
    # path lengths
    len_diversity = min(len(subset.keys()), len_diversity)
    path_lengths = list(subset.keys())[:len_diversity]
    # convert ids from metquest back to original using the name_map
    return [
        [name_map.get(p) for p in path]
        for path_length in path_lengths
        for path in subset.get(path_length)
    ]


def get_all_paths(mq_data, metabolites, tag, max_paths_length=40, len_diversity=3):
    met_paths = {}
    comments = {}
    name_map = mq_data.get("name_map")

    for met in metabolites:
        # metquest adds model id to the beginning (sometimes - TODO)
        met_tag = tag + " " + met

        # Get linear and circular paths from metquest output for a given metabolite
        linear_paths = mq_data.get(f"{max_paths_length}_linear_paths").get(met_tag)
        circular_paths = mq_data.get(f"{max_paths_length}_circular_paths").get(met_tag)

        # if metabolite is not in all_synthesized list, then it cannot be synthesized
        if met_tag not in mq_data.get("all_synthesized"):
            comments[met] = "can not be synthesized"

        # if metabolite is not in X_synthesized then it cannot be synthesized with max path length of X
        elif met_tag not in mq_data.get(f"{max_paths_length}_synthesized"):
            comment = f"can not be synthesized with max {max_paths_length} length path"
            comments[met] = comment

        # Check linear paths for given metabolite
        elif linear_paths is not None:
            if 0 in linear_paths:
                comments[met] = "No need to synthesize"
            else:
                met_paths[met] = get_paths(linear_paths, name_map, len_diversity)

        # Check circular paths for given metabolite
        elif circular_paths is not None:
            if 0 in circular_paths:
                comments[met] = "No need to synthesize"
            else:
                met_paths[met] = get_paths(circular_paths, name_map, len_diversity)

        # Otherwise check ...
        else:
            all_syn = set(mq_data.get("all_synthesized"))
            max_syn = set(mq_data.get(f"{max_paths_length}_synthesized"))
            maybe_synthesised = all_syn - max_syn
            if maybe_synthesised:  # if non-empty
                comment = "Problematic: can be synthesized but does not have paths"
            else:
                comment = (
                    f"Maybe can not be synthesized with max {max_paths_length} length "
                    f"path, because all_synthesized not bigger "
                    f"{max_paths_length}_synthesized"
                )
            comments[met] = comment
    return met_paths, comments


def preprocess_medium(
    tag: str, nutritional_sources: [str], other_medium: [str], cofactors: [str] = None
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


def main(
    path_to_model,
    path_to_union_model,
    output_dir,
    nutritional_sources: [str],
    other_medium: [str],
    cofactors=None,
    metabolites: [str] = None,
    max_paths_length=40,
    len_diversity=3,
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
    union_model = read_sbml_model(path_to_union_model)
    tag = model.id

    # If type of metabolites is given, rather than list of metabolites, extract
    # from union_model
    if isinstance(metabolites, str):
        metabolites = [
            m.id for m in union_model.reactions.get_by_id(metabolites).reactants
        ]

    # Alter metabolite ids
    seed_metabolites_tag, cofactors_c, other_medium_c = preprocess_medium(
        tag, nutritional_sources, other_medium
    )

    # Remove cofactors from the model
    model_wo_cof = remove_cofactors_from_model(
        deepcopy(model), other_medium_c, cofactors_c
    )

    # Save the model without cofactors
    path_to_model_wo_cof = output_dir / (model_name + "_wo_cofactors.xml")
    write_sbml_model(model_wo_cof, path_to_model_wo_cof)

    # TODO: save config to output_dir

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

    # Save the output
    write_mq_data(mq_output, output_dir / "metquest.h5")

    with open(output_dir / "graph.pkl", "wb") as fh:
        pickle.dump(graph, fh)
