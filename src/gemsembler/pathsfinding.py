import argparse
import pickle
import tempfile
from ast import literal_eval
from contextlib import contextmanager
from copy import deepcopy
from os import chdir, getcwd
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


def write_mq_data(mq_data: dict, filename: Path):
    """
    Function to save output from metquest in hdf5 format.
    HDF5 format was chosen due to large size of the output from metquest, which
    caused crashes when reading the data.

    Parameters
    ----------
    mq_data : dict
        Dictionary containing output from metquest
    filename : str
        path to output file
    """

    # Open file handle with h5py
    with h5py.File(filename, "w") as fh:
        # Save the metadata at root of the hierarchy in the file
        fh.create_dataset("all_synthesized", data=list(mq_data["all_synthesized"]))
        fh.create_dataset("max_synthesized", data=list(mq_data["max_synthesized"]))
        fh.create_dataset("max_paths_length", data=mq_data["max_paths_length"])
        fh.create_dataset("tag", data=mq_data["tag"])

        # Save the dictionary containing the mapping between metquest internal
        # ids and original ids
        name_map = fh.create_group("name_map")
        for mq_id, orig_id in mq_data["name_map"].items():
            name_map[mq_id] = orig_id

        # HDF5 format does not support saving dictionaries, so they need to
        # converted to strings before being saved. To ease memory footprint, we
        # stringify data at path length level, i.e.
        # /path_type/mb_id/path_length/<list_of_paths>
        # Applies for both linear and circular paths

        # Save linear linear paths
        linear_paths = fh.create_group("linear_paths")
        for mb_id, path_lengths in mq_data["linear_paths"].items():
            linear_paths.create_group(mb_id)
            for path_length, path in path_lengths.items():
                linear_paths[mb_id][str(path_length)] = str(path)

        # Save the circular paths
        circular_paths = fh.create_group("circular_paths")
        for mb_id, path_lengths in mq_data["circular_paths"].items():
            circular_paths.create_group(mb_id)
            for path_length, path in path_lengths.items():
                circular_paths[mb_id][str(path_length)] = str(path)


class MQData(object):
    def __init__(self, filename):
        # Save the filename as private variable to prevent the user from
        # accidentally overwriting it. Get and store metadata from the file.
        self.__filename = filename
        with h5py.File(self.__filename, "r") as fh:
            self.__max_paths_length = fh["max_paths_length"][()]
            self.__tag = fh["tag"][()].decode()

            self.__name_map = {
                mq_id: orig_id[()].decode() for mq_id, orig_id in fh["name_map"].items()
            }

            self.__rev_name_map = {
                orig_id: mq_id for mq_id, orig_id in self.__name_map.items()
            }

    # To minimise memory footprint, only load this data when requested by user
    @property
    def all_synthesized(self):
        with h5py.File(self.__filename, "r") as fh:
            return {x.decode() for x in fh["all_synthesized"]}

    # To minimise memory footprint, only load this data when requested by user
    @property
    def max_synthesized(self):
        with h5py.File(self.__filename, "r") as fh:
            return {x.decode() for x in fh["max_synthesized"]}

    # Property to prevent the user from overwriting it
    @property
    def name_map(self):
        return self.__name_map

    # Property to prevent the user from overwriting it
    @property
    def reverse_name_map(self):
        return self.__rev_name_map

    # Property to prevent the user from overwriting it
    @property
    def max_paths_length(self):
        return self.__max_paths_length

    # Property to prevent the user from overwriting it
    @property
    def tag(self):
        return self.__tag

    # Function to get available keys at the given level
    def get_keys(self, path_type, mb_id=None):
        with h5py.File(self.__filename, "r") as fh:
            if mb_id is not None:
                return list(fh[f"{path_type}/{mb_id}"].keys())
            else:
                return list(fh[path_type].keys())

    # Wrapper function around get_keys for better ux
    def get_ids(self, path_type):
        return self.get_keys(path_type)

    # Wrapper function around get_keys for better ux
    def get_paths_lengths(self, path_type, mb_id):
        return self.get_keys(path_type, mb_id)

    # Function to get all paths given path type, metabolite id and path length
    def __get_by_path_length(self, path_type, mb_id, path_length, convert=False):
        with h5py.File(self.__filename) as fh:
            paths = literal_eval(fh[f"{path_type}/{mb_id}/{path_length}"][()].decode())
            if convert:
                paths = [{self.__name_map[_id] for _id in path} for path in paths]
            return paths

    # Function to get all paths for all path lengths given path type and
    # metabolite id
    def __get_by_mb_id(self, path_type, mb_id, convert=False):
        return {
            path_len: self.__get_by_path_length(path_type, mb_id, path_len, convert)
            for path_len in self.get_paths_lengths(path_type, mb_id)
        }

    # Defining __getitem__ method to be able to access the data with
    # mq_data[<path_type>, <mb_id>, <path_length>]
    def __getitem__(self, key):
        # If only path type is given, return all paths for all metabolites
        if isinstance(key, str):
            return {
                mb_id: self.__get_by_mb_id(key, mb_id) for mb_id in self.get_ids(key)
            }
        # If more arguments were specified
        elif hasattr(key, "__iter__"):
            # Return all paths for all path lengths given metabolite id
            if len(key) == 2:
                return self.__get_by_mb_id(*key)
            # Return all paths for given path type, metabolite id and path length
            elif len(key) == 3:
                return self.__get_by_path_length(*key)
            else:
                raise ValueError(f"Too many args: {len(key)}")

    def _subset_paths(
        self, path_type, mb_id, path_lengths=None, len_diversity=None, map_back=True
    ):
        # Get path lengths if not defined, otherwise check the correct type
        if path_lengths is None:
            # path lengths are stored as strings in the file, so needs
            # conversion back to ints and sorting
            path_lengths = sorted(
                [int(x) for x in self.get_paths_lengths(path_type, mb_id)]
            )
        elif isinstance(path_lengths, int):
            path_lengths = [path_lengths]
        elif isinstance(path_lengths, list):
            pass
        else:
            raise ValueError("path_lengths can only be: None, int or list")

        # Subset path lengths based on len_diversity arg if given
        if len_diversity is not None:
            len_diversity = min(len(path_lengths), len_diversity)
            path_lengths = path_lengths[:len_diversity]

        # Return all paths for given path lengths (if specified)
        # In case paths are there: length bigger than 0
        # otherwise return empty dict
        if 0 in path_lengths:
            return {0: []}
        else:
            return {
                path_len: self.__get_by_path_length(
                    path_type, mb_id, path_len, map_back
                )
                for path_len in path_lengths
            }

    # Wrapper function around _subset_paths to allow for getting all paths for
    # all metabolites
    def subset_paths(
        self,
        path_type,
        mb_id=None,
        path_lengths=None,
        len_diversity=None,
        map_back=True,
    ):
        if mb_id is None:
            return {
                self._subset_paths(
                    path_type, mb_id, path_lengths, len_diversity, map_back
                )
                for mb_id in self.get_ids(path_type)
            }
        else:
            return self._subset_paths(
                path_type, mb_id, path_lengths, len_diversity, map_back
            )

    # metquest is not consistent with metabolite ids, hence we get all
    # possibilites
    def _get_alt_mb_ids(self, mb_id):
        mb_id_tag = self.__tag + " " + mb_id
        mb_id_name_map = self.__rev_name_map.get(mb_id, "")
        return mb_id_tag, mb_id_name_map

    # Function to check which of the given metabolites is present in set of ids
    # for given path type. args is a list of possible metabolite ids.
    def _check_mb_id(self, path_type, *args):
        pos_ids = self.get_ids(path_type)
        for mb_id in args:
            if mb_id in pos_ids:
                return mb_id

    def get_all_paths(self, metabolites, len_diversity=3):
        met_paths = {}
        comments = {}
        all_syn = self.all_synthesized
        max_syn = self.max_synthesized

        for met in metabolites:
            # metquest adds model id to the beginning (sometimes - TODO)
            # Get all possible ways metquest can give ids to metabolites
            met_tag, met_name_map = self._get_alt_mb_ids(met)

            # if metabolite is not in all_synthesized list, then it cannot be
            # synthesized
            if (
                (met_tag not in all_syn)
                and (met not in all_syn)
                and (met_name_map not in all_syn)
            ):
                comments[met] = "can not be synthesized"
                met_paths[met] = {}

            # if metabolite is not in X_synthesized then it cannot be
            # synthesized with max path length of X
            elif (
                (met_tag not in max_syn)
                and (met not in max_syn)
                and (met_name_map not in max_syn)
            ):
                comment = f"can not be synthesized with max {self.max_paths_length} length path"
                comments[met] = comment
                met_paths[met] = {}
            else:
                met_found = self._check_mb_id(
                    "linear_paths", met_tag, met, met_name_map
                )
                if met_found is not None:
                    # Check linear paths for given metabolite. Get linear and
                    # circular paths from metquest output for a given metabolite
                    linear_paths = self._subset_paths(
                        "linear_paths", met_found, len_diversity=len_diversity
                    )
                    if linear_paths == {0: []}:
                        comments[met] = "Is in the medium, no need to synthesys"
                        met_paths[met] = {}
                    elif isinstance(linear_paths, str):
                        comments[met] = linear_paths
                    elif isinstance(linear_paths, dict):
                        comments[met] = "Found in linear paths"
                        met_paths[met] = linear_paths
                else:
                    met_found = self._check_mb_id(
                        "circular_paths", met_tag, met, met_name_map
                    )
                    if met_found is not None:
                        # Check circular paths for given metabolite
                        circular_paths = self._subset_paths(
                            "circular_paths", met_found, len_diversity=len_diversity
                        )
                        if circular_paths == {0: []}:
                            comments[met] = "Is in the medium, no need to synthesys"
                            met_paths[met] = {}
                        elif isinstance(circular_paths, str):
                            comments[met] = circular_paths
                        elif isinstance(circular_paths, dict):
                            comments[met] = "Found in circular paths"
                            met_paths[met] = circular_paths
                    else:
                        # Otherwise check ...
                        maybe_synthesised = self.all_synthesized - self.max_synthesized
                        if maybe_synthesised:  # if non-empty
                            comment = (
                                "Problematic: can be synthesized "
                                "but does not have paths"
                            )
                        else:
                            comment = (
                                "Maybe can not be synthesized with max "
                                f"{self.max_paths_length} length path, because "
                                "all_synthesized not bigger than max_synthesized"
                            )
                        comments[met] = comment
                        met_paths[met] = {}

        return met_paths, comments


def check_which_met_id_found(path_data: dict, *args):
    pos_ids = path_data.keys()
    for mb_id in args:
        if mb_id in pos_ids:
            return mb_id


def get_paths_for_met(path_data: dict, met: str, len_diversity: int, name_map: dict):
    len_diversity = min(len(path_data[met]), len_diversity)
    path_lengths = list(path_data[met].keys())[:len_diversity]
    # Return all paths for given path lengths (if specified)
    # In case paths are there: length bigger than 0, otherwise return empty dict
    if 0 in path_lengths:
        return {0: []}
    met_paths = {}
    for path_len in path_lengths:
        paths = path_data[met][path_len]
        new_paths = []
        for path in paths:
            new_paths.append({name_map[r_id] for r_id in path})
        met_paths[path_len] = new_paths
    return met_paths


def get_shortest_paths_data(mq_out_data: dict, metabolites: list, len_diversity=3):
    met_paths = {}
    comments = {}
    all_syn = mq_out_data["all_synthesized"]
    max_syn = mq_out_data["max_synthesized"]

    for met in metabolites:
        # Get all possible ways metquest can give ids to metabolites
        met_tag = mq_out_data["tag"] + " " + met
        reverse_name_map = {
            orig_id: mq_id for mq_id, orig_id in mq_out_data["name_map"].items()
        }
        met_name_map = reverse_name_map.get(met, "")

        # if metabolite is not in all_synthesized list, then it cannot be
        # synthesized
        if (
            (met_tag not in all_syn)
            and (met not in all_syn)
            and (met_name_map not in all_syn)
        ):
            comments[met] = "can not be synthesized"
            met_paths[met] = {}

        # if metabolite is not in X_synthesized then it cannot be
        # synthesized with max path length of X
        elif (
            (met_tag not in max_syn)
            and (met not in max_syn)
            and (met_name_map not in max_syn)
        ):
            comment = f"can not be synthesized with max {mq_out_data['max_paths_length']} length path"
            comments[met] = comment
            met_paths[met] = {}
        else:
            met_found = check_which_met_id_found(
                mq_out_data["linear_paths"], met, met_tag, met_name_map
            )
            if met_found is not None:
                # Check linear paths for given metabolite. Get linear and
                # circular paths from metquest output for a given metabolite
                linear_paths = get_paths_for_met(
                    mq_out_data["linear_paths"],
                    met_found,
                    len_diversity=len_diversity,
                    name_map=mq_out_data["name_map"],
                )
                if linear_paths == {0: []}:
                    comments[met] = "Is in the medium, no need to synthesys"
                    met_paths[met] = {}
                elif isinstance(linear_paths, dict):
                    comments[met] = "Found in linear paths"
                    met_paths[met] = linear_paths
            else:
                c_met_found = check_which_met_id_found(
                    mq_out_data["circular_paths"], met, met_tag, met_name_map
                )
                if c_met_found is not None:
                    # Check circular paths for given metabolite
                    circular_paths = get_paths_for_met(
                        mq_out_data["circular_paths"],
                        c_met_found,
                        len_diversity=len_diversity,
                        name_map=mq_out_data["name_map"],
                    )
                    if circular_paths == {0: []}:
                        comments[met] = "Is in the medium, no need to synthesys"
                        met_paths[met] = {}
                    elif isinstance(circular_paths, dict):
                        comments[met] = "Found in circular paths"
                        met_paths[met] = circular_paths
                else:
                    # Otherwise check ...
                    maybe_synthesised = all_syn - max_syn
                    if maybe_synthesised:  # if non-empty
                        comment = (
                            "Problematic: can be synthesized " "but does not have paths"
                        )
                    else:
                        comment = (
                            "Maybe can not be synthesized with max "
                            f"{mq_out_data['max_paths_length']} length path, because "
                            "all_synthesized not bigger than max_synthesized"
                        )
                    comments[met] = comment
                    met_paths[met] = {}
    return {"paths": met_paths, "comments": comments}


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
    path_to_models: str, seed_metabolites_tag: [str], number_of_xml, max_paths_length,
):
    graph, name_map = metquest.construct_graph.create_graph(
        path_to_models, number_of_xml
    )

    lowboundmet_pic, status_dict_pic, scope_pic = metquest.forward_pass(
        graph, seed_metabolites_tag
    )

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
    len_diversity=3,
    metabolites=None,
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
            tmp_dir, seed_metabolites_tag, number_of_xml, max_paths_length,
        )

    # Add necessary parameters to output data
    mq_output.update({"tag": tag, "max_paths_length": max_paths_length})

    # Save the output
    write_mq_data(mq_output, output_dir / "metquest.h5")

    if metabolites is None:
        metabolites = [m.id for m in model_wo_cof.metabolites]
    shortest_paths = get_shortest_paths_data(mq_output, metabolites, len_diversity)

    with open(output_dir / "shortest_paths.pkl", "wb") as fh:
        pickle.dump(shortest_paths, fh)

    with open(output_dir / "graph.pkl", "wb") as fh:
        pickle.dump(graph, fh)
