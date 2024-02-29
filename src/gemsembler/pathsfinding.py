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
        fh.create_dataset("max_paths_length", data=mq_data["max_paths_length"])
        fh.create_dataset("tag", data=mq_data["tag"])

        name_map = fh.create_group("name_map")
        for mq_id, orig_id in mq_data["name_map"].items():
            name_map[mq_id] = orig_id

        linear_paths = fh.create_group("linear_paths")
        for r_id, path_lengths in mq_data["linear_paths"].items():
            linear_paths.create_group(r_id)
            for path_length, path in path_lengths.items():
                linear_paths[r_id][str(path_length)] = str(path)

        circular_paths = fh.create_group("circular_paths")
        for r_id, path_lengths in mq_data["circular_paths"].items():
            circular_paths.create_group(r_id)
            for path_length, path in path_lengths.items():
                circular_paths[r_id][str(path_length)] = str(path)


class MQData(object):
    def __init__(self, filename):
        self.__filename = filename

    @property
    def all_synthesized(self):
        with h5py.File(self.__filename, "r") as fh:
            return [x.decode() for x in fh["all_synthesized"]]

    @property
    def max_synthesized(self):
        with h5py.File(self.__filename, "r") as fh:
            return [x.decode() for x in fh["max_synthesized"]]

    @property
    def name_map(self):
        with h5py.File(self.__filename, "r") as fh:
            return {
                mq_id: orig_id[()].decode() for mq_id, orig_id in fh["name_map"].items()
            }

    @property
    def max_paths_length(self):
        with h5py.File(self.__filename, "r") as fh:
            return fh["max_paths_length"][()]

    @property
    def tag(self):
        with h5py.File(self.__filename, "r") as fh:
            return fh["tag"][()].decode()

    def __getitem__(self, key):
        with h5py.File(self.__filename, "r") as fh:
            if isinstance(key, str):
                return {
                    r_id: {
                        path_length: literal_eval(
                            fh[f"{key}/{r_id}/{path_length}"][()].decode()
                        )
                        for path_length, path in paths.items()
                    }
                    for r_id, paths in fh[key].items()
                }
            elif hasattr(key, "__iter__"):
                if len(key) == 2:
                    return {
                        path_length: literal_eval(path[()].decode())
                        for path_length, path in fh["/".join(key)].items()
                    }
                elif len(key) == 3:
                    return literal_eval(fh[f"/".join(key)][()].decode())
                else:
                    raise ValueError(f"Too many args: {len(key)}")

    def get_keys(self, path_type, r_id=None):
        with h5py.File(self.__filename, "r") as fh:
            if r_id is not None:
                return list(fh[f"{path_type}/{r_id}"].keys())
            else:
                return list(fh[path_type].keys())

    def get_ids(self, path_type):
        return self.get_keys(path_type)

    def get_paths_lengths(self, path_type, r_id):
        return self.get_keys(path_type, r_id)

    def _get_paths(
        self, path_type, r_id, path_lengths=None, len_diversity=None, map_back=True
    ):
        # Get path lengths if not defined, otherwise check the correct type
        if path_lengths is None:
            path_lengths = get_paths_lengths("linear_paths", r_id)
            len_diversity = min(len(path_lengths), len_diversity)
            path_lengths = path_lengths[:len_diversity]
        elif isinstance(path_lengths, int):
            path_lengths = [path_lengths]
        elif isinstance(path_lengths, list):
            pass
        else:
            raise ValueError("path_lengths can only be: None, int or list")

        with h5py.File(filename, "r") as fh:
            if map_back:
                return {
                    k: [
                        [fh["name_map"][p][()].decode() for p in paths]
                        for path_length in path_lengths
                        for paths in literal_eval(
                            fh[f"{path_type}/{r_id}/{path_length}"][()].decode()
                        )
                    ]
                }
            else:
                return {
                    k: [
                        [p for p in paths]
                        for path_length in path_lengths
                        for paths in literal_eval(
                            fh[f"{path_type}/{r_id}/{path_length}"][()].decode()
                        )
                    ]
                }

        return self.__getitem__("linear_paths")

    def get_linear_paths(
        self, r_id, path_lengths=None, len_diversity=None, map_back=True
    ):
        return self._get_paths(
            "linear_paths", r_id, path_lengths, len_diversity, map_back
        )

    def get_circular_paths(
        self, r_id, path_length=None, len_diversity=None, map_back=True
    ):
        return self._get_paths(
            "circular_paths", r_id, path_lengths, len_diversity, map_back
        )


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

    tag = read_mq_data(mq_data_path, "tag")[0]
    max_paths_length = read_mq_data(mq_data_path, "max_paths_length")[0]
    name_map = read_mq_data(mq_data_path, "name_map")
    all_syn = read_mq_data(mq_data_path, "all_synthesized")
    max_syn = read_mq_data(mq_data_path, "max_synthesized")

    for met in metabolites:
        # metquest adds model id to the beginning (sometimes - TODO)
        met_tag = tag + " " + met
        met_name_map = findKeysByValue(name_map, met)
        if met_name_map:
            met_name_map = met_name_map[0]
        else:
            met_name_map = ""

        # if metabolite is not in all_synthesized list, then it cannot be synthesized
        if (
            (met_tag not in all_syn)
            and (met not in all_syn)
            and (met_name_map not in all_syn)
        ):
            comments[met] = "can not be synthesized"

        # if metabolite is not in X_synthesized then it cannot be synthesized with max path length of X
        elif (
            (met_tag not in max_syn)
            and (met not in max_syn)
            and (met_name_map not in max_syn)
        ):
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
            linear_paths = read_mq_data(
                mq_data_path, "max_linear_paths", met_found, len_diversity
            )
            if type(linear_paths[met_found]) == str:
                comments[met] = linear_paths[met_found]
            elif type(linear_paths[met_found]) == list:
                met_paths[met] = linear_paths[met_found]
            elif linear_paths[met_found] is None:
                # Check circular paths for given metabolite
                circular_paths = read_mq_data(
                    mq_data_path, "max_circular_paths", met_found, len_diversity
                )
                if type(circular_paths[met_found]) == str:
                    comments[met] = circular_paths[met_found]
                elif type(circular_paths[met_found]) == list:
                    met_paths[met] = circular_paths[met_found]
                elif circular_paths[met_found] is None:
                    # Otherwise check ...
                    maybe_synthesised = set(all_syn) - set(max_syn)
                    if maybe_synthesised:  # if non-empty
                        comment = (
                            "Problematic: can be synthesized but does not have paths"
                        )
                    else:
                        comment = (
                            f"Maybe can not be synthesized with max {max_paths_length} length "
                            f"path, because all_synthesized not bigger than max_synthesized"
                        )
                    comments[met] = comment

    return met_paths, comments


def preprocess_medium(
    tag: str,
    nutritional_sources: [str],
    other_medium: [str],
    cofactors: [str],
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
        "linear_paths": pathway_table,
        "circular_paths": c_path,
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
