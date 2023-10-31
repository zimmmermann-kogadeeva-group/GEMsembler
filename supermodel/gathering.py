import sys
import warnings
import logging
import operator
from collections import Counter, defaultdict
from copy import deepcopy

from cobra.io import read_sbml_model

from .conversion import ConvCarveme, ConvGapseq, ConvModelseed, ConvAgora
from .curation import remove_b_type_exchange, get_duplicated_reactions
from .general import findKeysByValue
from .selection import checkDBConsistency, checkFromOneFromMany
from .structural import runStructuralConversion, runStructuralCheck
from .dbs import get_bigg_network


class LoggerContext:
    def __init__(self, logger_name, show_logs=False):
        self.__show_logs__ = show_logs
        self.__logger__ = logging.getLogger(logger_name)

    def __enter__(self):
        if not self.__show_logs__:
            self.__logger__.setLevel(logging.CRITICAL)

    def __exit__(self, exit_type, exit_value, exit_traceback):
        if not self.__show_logs__:
            self.__logger__.setLevel(logging.NOTSET)


class GatheredModels:
    """
    Class that gathers information and necessary conversion results for all
    models. Input for the class and tool in general is dictionary with all
    models and related information.

    Parameters
    ----------
    dict_of_all_models_with_feature : dict
        Dictionary of the following format:
        {
            model_id:
            {
                'path_to_model': str,
                'model_type': str one of (agora, carveme, gapseq, modelseed) or custom,
                'path_to_genome': str (can be '' or None if convert_genes = False)
            }
        }
        Note: if custom `model_type` create type class in advance,
    assembly : None, optional
    path_final_genome_nt : None, optional
    path_final_genome_aa : None, optional
    custom_model_type : None, optional

    Notes
    -----
    If all 3 parameters bellow are None then gene conversion is not done and
    genomes for model_id are not need.
    """

    def __init__(
        self,
        assembly_id=None,
        path_final_genome_nt=None,
        path_final_genome_aa=None,
        custom_model_type=None,
    ):
        self.__conf__ = {
            "carveme": {
                "remove_b": False,
                "db_name": "bigg",
                "wo_periplasmic": False,
                "conv_strategy": ConvCarveme(),
            },
            "gapseq": {
                "remove_b": False,
                "db_name": "modelseed",
                "wo_periplasmic": False,
                "conv_strategy": ConvGapseq(),
            },
            "modelseed": {
                "remove_b": True,
                "db_name": "modelseed",
                "wo_periplasmic": False,
                "conv_strategy": ConvModelseed(),
            },
            "agora": {
                "remove_b": False,
                "db_name": "weird_bigg",
                "wo_periplasmic": True,
                "conv_strategy": ConvAgora(),
            },
        }
        self.__models__ = {}
        self.first_selected = None

        # Check if assembly and final genome are present.
        # If not, throw a warning.
        if not assembly_id and not path_final_genome_nt and not path_final_genome_aa:
            warnings.warn(
                "\nWarning! No final genome for gene conversion is provided. "
                "Gene conversion will not be performed.\nIf you want to "
                "convert genes, please provide either assembly id or custom "
                "fasta files (nt/aa/both), \nto wich genes must be converted."
            )
            convert_genes = False
        else:
            convert_genes = True

    def _get_same_db_models(self):
        same_db_models = defaultdict(dict)
        for model_id, model_attrs in self.__models__.items():
            model_type = model_attrs["model_type"]
            db_name = self.__conf__.get(model_type).get("db_name")
            same_db_models[db_name][model_id] = model_type
        return same_db_models

    def run(self):
        same_db_models = self._get_same_db_models()

        self.first_selected = checkDBConsistency(
            same_db_models,
            {
                model_id: model_attrs["converted"]
                for model_id, model_attrs in self.__models__.items()
            },
            "highest",
        )

        for s in self.first_selected.values():
            checkFromOneFromMany(s)

        self.structural = {}
        for model_id, sel in self.first_selected.items():
            checkFromOneFromMany(sel)
            if (
                strategies.db_name[
                    dict_of_all_models_with_feature[model_id]["model_type"]
                ]
                == "bigg"
            ):
                self.structural.update(
                    {
                        model_id: {
                            "reactions": runStructuralCheck(
                                self.first_selected[model_id]["reactions"],
                                self.preprocessed_models[model_id],
                                bigg_network,
                            )
                        }
                    }
                )
            else:
                self.structural.update(
                    {
                        model_id: {
                            "reactions": runStructuralConversion(
                                self.first_selected[model_id],
                                self.preprocessed_models[model_id],
                                bigg_network,
                                strategies.wo_periplasmic[
                                    dict_of_all_models_with_feature[model_id][
                                        "model_type"
                                    ]
                                ],
                            )
                        }
                    }
                )

    def set_configuration(
        self,
        model_type: str,
        remove_b: bool,
        db_name: str,
        wo_periplasmic: bool,
        conv_strategy,
        **kwargs,
    ):
        assert isinstance(conv_strategy, ConvBase)

        # TODO: add checks on conf input arg
        self.__conf__[model_type] = {
            "remove_b": remove_b,
            "db_name": db_name,
            "wo_periplasmic": wo_periplasmic,
            "conv_strategy": conv_strategy,
            **kwargs,
        }

    def add_model(
        self,
        model_id: str,
        path_to_model: str,
        model_type: str,
        path_to_genome: str = None,
        show_logs: bool = False,
    ):
        # Run checks on model_id and model_type
        # TODO: check with conversion dictionaries
        assert model_id not in self.__models__, f"model_id {model_id} already used"
        assert model_type in self.__conf__, f"Missing configuration for {model_type}"

        # Read the cobra model
        with LoggerContext("cobra", show_logs):
            model = read_sbml_model(path_to_model)

        self.__models__[model_id] = {
            "original_model": deepcopy(model),
            "path_to_model": path_to_model,
            "model_type": model_type,
            "path_to_genome": path_to_genome,
        }

        # If model_type requires it, remove `_b` extensions
        if self.__conf__.get(model_type).get("remove_b"):
            model = remove_b_type_exchange(model)
        self.__models__[model_id]["model"] = model

        dupl_r, dupl_r_gpr = get_duplicated_reactions(model)
        self.__models__[model_id]["duplicated_r"] = dupl_r

        # TODO: add isinstance
        conv = self.__conf__.get(model_type).get("conv_strategy")
        self.__models__[model_id]["converted"] = (
            self.__conf__.get(model_type).get("conv_strategy").convert_model(model)
        )
