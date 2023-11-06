from collections import Counter, defaultdict
from copy import deepcopy
import logging
import operator
from pathlib import Path
import pickle
import sys
import warnings

from cobra.io import read_sbml_model

from .conversion import ConvCarveme, ConvGapseq, ConvModelseed, ConvAgora, ConvBase
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
                "wo_periplasmic": True,
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
        self.converted_metabolites = defaultdict(dict)
        self.converted_reactions = defaultdict(dict)
        self.first_stage_selected_metabolites = None
        self.first_stage_selected_reactions = None
        self.structural_first_run_reactions = defaultdict(dict)

        # Check if assembly and final genome are present.
        # If not, throw a warning.
        if not assembly_id and not path_final_genome_nt and not path_final_genome_aa:
            warnings.warn(
                "\nWarning! No final genome for gene conversion is provided. "
                "Gene conversion will not be performed.\nIf you want to "
                "convert genes, please provide either assembly id or custom "
                "fasta files (nt/aa/both), \nto which genes must be converted."
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
        # run conversion
        for model_id, model_attrs in self.__models__.items():
            conv = self.__conf__.get(model_attrs["model_type"]).get("conv_strategy")
            self.converted_metabolites[model_id] = {
                m.id: conv.convert_metabolite(m)
                for m in model_attrs["preprocess_model"].metabolites
            }
            self.converted_reactions[model_id] = {
                r.id: conv.convert_reaction(r)
                for r in model_attrs["preprocess_model"].reactions
            }

        # prepare dictionary for models per used database
        same_db_models = self._get_same_db_models()

        # run first stage selection
        self.first_stage_selected_metabolites = checkDBConsistency(
            same_db_models,
            self.converted_metabolites,
            "highest",
        )
        for s in self.first_stage_selected_metabolites.values():
            checkFromOneFromMany(s)

        self.first_stage_selected_reactions = checkDBConsistency(
            same_db_models,
            self.converted_reactions,
            "highest",
        )
        for s in self.first_stage_selected_reactions.values():
            checkFromOneFromMany(s)

        # run first structural conversion
        bigg_network = get_bigg_network()
        for model_id, first_sel in self.first_stage_selected_reactions.items():
            model_type = self.__models__[model_id]["model_type"]
            db_name = self.__conf__.get(model_type).get("db_name")
            if db_name == "bigg":
                self.structural_first_run_reactions[model_id] = runStructuralCheck(
                    first_sel,
                    self.__models__[model_id]["preprocess_model"],
                    bigg_network,
                )
            else:
                self.structural_first_run_reactions[model_id] = runStructuralConversion(
                    first_sel,
                    self.first_stage_selected_metabolites[model_id],
                    self.__models__[model_id]["preprocess_model"],
                    bigg_network,
                    self.__conf__.get(model_type).get("wo_periplasmic"),
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
        cache: bool = True,
        show_logs: bool = False,
    ):
        # Run checks on model_id and model_type
        # TODO: check with conversion dictionaries
        assert model_id not in self.__models__, f"model_id {model_id} already used"
        assert model_type in self.__conf__, f"Missing configuration for {model_type}"

        cache_path = Path(path_to_model).with_suffix(".p")
        if cache_path.exists():
            with open(cache_path, "rb") as fh:
                model = pickle.load(fh)
        else:
            # Read the cobra model
            with LoggerContext("cobra", show_logs):
                model = read_sbml_model(path_to_model)
            # Cache the model to a pickle if option is set
            if cache:
                with open(cache_path, "wb") as fh:
                    pickle.dump(model, fh)

        # Populate the internal data
        self.__models__[model_id] = {
            "original_model": deepcopy(model),
            "path_to_model": path_to_model,
            "model_type": model_type,
            "path_to_genome": path_to_genome,
        }

        # If model_type requires it, remove `_b` extensions
        if self.__conf__.get(model_type).get("remove_b"):
            model = remove_b_type_exchange(model)
        self.__models__[model_id]["preprocess_model"] = model

        dupl_r, dupl_r_gpr = get_duplicated_reactions(model)
        self.__models__[model_id]["duplicated_reactions"] = dupl_r
