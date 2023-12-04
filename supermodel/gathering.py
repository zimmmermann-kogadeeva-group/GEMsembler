import os
from collections import defaultdict
from copy import deepcopy
import logging
from pathlib import Path
import pickle
import sys
from typing import Type
import warnings

from cobra.io import read_sbml_model
from .conversion import ConvCarveme, ConvGapseq, ConvModelseed, ConvAgora, ConvBase
from .curation import remove_b_type_exchange, get_duplicated_reactions
from .periplasmic import getSuggestionPeriplasmic
from .selection import run_selection
from .structural import runStructuralConversion, runSuggestionsMet
from .dbs import get_bigg_network
from .genes import (
    get_final_fasta_with_ncbi_assemble,
    get_genes_gapseq,
    get_genes_not_gapseq,
)


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


def load_sbml_model(path_to_model, cache: bool = True, show_logs: bool = False):
    cache_path = Path(path_to_model).with_suffix(".pkl")
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
    return model


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
        custom_model_type=None,
        clear_db_cache=False,
    ):
        # If specified, clear the cached conversion tables and dictionaries
        if clear_db_cache:
            for p in Path("~/.gemsembler/").expanduser().iterdir():
                p.unlink()

        self.__conf = {
            "agora": {
                "remove_b": False,
                "db_name": "weird_bigg",
                "wo_periplasmic": True,
                "conv_strategy": ConvAgora(),
                "genome_model_strategy": get_genes_not_gapseq,
            },
            "carveme": {
                "remove_b": False,
                "db_name": "bigg",
                "wo_periplasmic": False,
                "conv_strategy": ConvCarveme(),
                "genome_model_strategy": get_genes_not_gapseq,
            },
            "gapseq": {
                "remove_b": False,
                "db_name": "modelseed",
                "wo_periplasmic": False,
                "conv_strategy": ConvGapseq(),
                "genome_model_strategy": get_genes_gapseq,
            },
            "modelseed": {
                "remove_b": True,
                "db_name": "modelseed",
                "wo_periplasmic": True,
                "conv_strategy": ConvModelseed(),
                "genome_model_strategy": get_genes_not_gapseq,
            },
        }

        self.__models = {}
        self.converted_metabolites = defaultdict(dict)
        self.converted_reactions = defaultdict(dict)
        self.first_stage_selected_metabolites = None
        self.first_stage_selected_reactions = None
        self.structural_first_run_reactions = defaultdict(dict)
        self.second_stage_selected_reactions = None
        self.structural_first_run_metabolites = defaultdict(dict)
        self.second_stage_selected_metabolites = None
        self.structural_second_run_reactions = defaultdict(dict)
        self.third_stage_selected_reactions = None
        self.many_to_one_sug = defaultdict(dict)
        self.periplasmic_metabolites = defaultdict(dict)
        self.periplasmic_reactions = defaultdict(dict)

    def __contains__(self, item):
        return item in self.__models

    def __len__(self):
        return len(self.__models)

    def get_conf(self, model_type=None):
        if model_type is None:
            return deepcopy(self.__conf)
        else:
            return deepcopy(self.__conf.get(model_type))

    def get_model_attrs(self, model_id=None):
        if model_id is None:
            return deepcopy(self.__models)
        else:
            return deepcopy(self.__models.get(model_id))

    def _get_same_db_models(self):
        same_db_models = defaultdict(dict)
        for model_id, model_attrs in self.__models.items():
            model_type = model_attrs["model_type"]
            db_name = self.__conf.get(model_type).get("db_name")
            same_db_models[db_name][model_id] = model_type
        return same_db_models

    def run(self):
        # run conversion
        for model_id, model_attrs in self.__models.items():
            conv = self.__conf.get(model_attrs["model_type"]).get("conv_strategy")
            converted_model = conv.convert_model(model_attrs["preprocess_model"])

            self.converted_metabolites[model_id] = converted_model["metabolites"]
            self.converted_reactions[model_id] = converted_model["reactions"]

        # prepare dictionary for models per used database
        same_db_models = self._get_same_db_models()

        # run first stage selection for converted
        self.first_stage_selected_metabolites = run_selection(
            same_db_models, self.converted_metabolites, "highest"
        )
        self.first_stage_selected_reactions = run_selection(
            same_db_models, self.converted_reactions, "highest"
        )

        # run first structural conversion
        bigg_network = get_bigg_network()
        for model_id, first_sel in self.first_stage_selected_reactions.items():
            model_type = self.__models[model_id]["model_type"]
            db_name = self.__conf.get(model_type).get("db_name")
            self.structural_first_run_reactions[model_id] = runStructuralConversion(
                db_name,
                first_sel,
                self.first_stage_selected_metabolites[model_id],
                self.__models[model_id]["preprocess_model"],
                bigg_network,
                False,
            )
        # run second stage selection for first structural reactions
        self.second_stage_selected_reactions = run_selection(
            same_db_models,
            self.structural_first_run_reactions,
            "structural",
            replace_with_consistent=False,
        )

        # get suggestions from structural reactions for metabolites
        for model_id, rs_struct_sel in self.second_stage_selected_reactions.items():
            model_type = self.__models[model_id]["model_type"]
            db_mod = self.__conf.get(model_type).get("db_name")
            (
                self.structural_first_run_metabolites[model_id],
                self.many_to_one_sug[model_id],
            ) = runSuggestionsMet(
                db_mod,
                self.structural_first_run_reactions[model_id],
                rs_struct_sel,
                self.first_stage_selected_metabolites[model_id],
            )
        # run second stage selection for suggestions for metabolites from structural
        self.second_stage_selected_metabolites = run_selection(
            same_db_models, self.structural_first_run_metabolites, "structural"
        )

        # run second structural conversion with suggestions for metabolites
        for model_id, sec_sel in self.second_stage_selected_reactions.items():
            model_type = self.__models[model_id]["model_type"]
            db_name = self.__conf.get(model_type).get("db_name")
            self.structural_second_run_reactions[model_id] = runStructuralConversion(
                db_name,
                sec_sel,
                self.second_stage_selected_metabolites[model_id],
                self.__models[model_id]["preprocess_model"],
                bigg_network,
                self.__conf.get(model_type).get("wo_periplasmic"),
            )
        # run third stage selection for first structural reactions
        self.third_stage_selected_reactions = run_selection(
            same_db_models,
            self.structural_second_run_reactions,
            "structural",
            replace_with_consistent=False,
        )

        # introducing periplasmic compartment for models, that don't have it originally
        for model_id, th_sel in self.third_stage_selected_reactions.items():
            if self.__conf.get(self.__models[model_id]["model_type"]).get(
                "wo_periplasmic"
            ):
                (
                    self.periplasmic_metabolites[model_id],
                    self.periplasmic_reactions[model_id],
                ) = getSuggestionPeriplasmic(
                    th_sel,
                    self.structural_second_run_reactions[model_id],
                    self.second_stage_selected_metabolites[model_id],
                    self.__models[model_id]["preprocess_model"],
                    bigg_network,
                )
            else:
                (
                    self.periplasmic_metabolites[model_id],
                    self.periplasmic_reactions[model_id],
                ) = ({}, {})

    def get_input_dictionaries(self):
        final_r_sel = defaultdict(dict)
        final_r_not_sel = defaultdict(dict)
        final_m_sel = defaultdict(dict)
        final_m_not_sel = defaultdict(dict)
        periplasmic_m = defaultdict(dict)
        for model_id in self.__models.keys():
            for orig_r_id, sel_r in self.third_stage_selected_reactions[
                model_id
            ].items():
                if (sel_r.to_one_id == True and sel_r.from_one_id == True) or (
                    sel_r.to_one_id == True
                    and sel_r.from_one_id == False
                    and orig_r_id
                    in self.__models[model_id]["duplicated_reactions"]["ID"].tolist()
                ):
                    final_r_sel[model_id].update(
                        {orig_r_id: [sel_r.compartments, sel_r.highest_consistent]}
                    )
                else:
                    final_r_not_sel[model_id].update(
                        {orig_r_id: [sel_r.compartments, sel_r.highest_consistent]}
                    )
                if (
                    len(
                        self.__models[model_id]["preprocess_model"]
                        .reactions.get_by_id(orig_r_id)
                        .reactants
                    )
                    > 24
                ):
                    final_r_not_sel[model_id].update(
                        {orig_r_id: [sel_r.compartments, sel_r.highest_consistent]}
                    )
            for orig_m_id, sel_m in self.second_stage_selected_metabolites[
                model_id
            ].items():
                if orig_m_id in self.periplasmic_metabolites[model_id]:
                    if self.periplasmic_metabolites[model_id][orig_m_id].replace:
                        final_m_sel[model_id].update(
                            {
                                orig_m_id: [
                                    ["p"],
                                    self.periplasmic_metabolites[model_id][
                                        orig_m_id
                                    ].bigg_p,
                                ]
                            }
                        )
                    else:
                        final_m_sel[model_id].update(
                            {orig_m_id: [sel_m.compartments, sel_m.highest_consistent]}
                        )
                        periplasmic_m[model_id].update(
                            {
                                orig_m_id: [
                                    ["p"],
                                    self.periplasmic_metabolites[model_id][
                                        orig_m_id
                                    ].bigg_p,
                                ]
                            }
                        )
                else:
                    if sel_m.to_one_id == True and sel_m.from_one_id == True:
                        final_m_sel[model_id].update(
                            {orig_m_id: [sel_m.compartments, sel_m.highest_consistent]}
                        )
                    else:
                        final_m_not_sel[model_id].update(
                            {orig_m_id: [sel_m.compartments, sel_m.highest_consistent]}
                        )
        return final_m_sel, final_m_not_sel, final_r_sel, final_r_not_sel, periplasmic_m

    def assemble_supermodel(
        self,
        output_folder: str,
        assembly_id=None,
        path_final_genome_nt=None,
        path_final_genome_aa=None,
    ):
        # Check if assembly and final genome are present.
        # If not, throw a warning.
        if not assembly_id and not path_final_genome_nt and not path_final_genome_aa:
            warnings.warn(
                "\nWarning! No final genome for gene conversion is provided. "
                "Gene conversion will not be performed.\nIf you want to "
                "convert genes, please provide either assembly id or custom "
                "fasta files (nt/aa/both), \nto which genes must be converted."
            )
        elif assembly_id and (path_final_genome_nt or path_final_genome_aa):
            warnings.warn(
                "\nWarning! Both assembly and user final genome for gene conversion are provided. "
                "Gene conversion will not be performed.\nIf you want to "
                "convert genes, please provide one of both either assembly id or custom "
                "fasta files (nt/aa/both), to which genes must be converted."
            )
        else:
            gene_path = Path(output_folder, "tmp_gene_conversion")
            Path(gene_path).mkdir(exist_ok=True)
            db_path = Path(gene_path, "blast_db")
            Path(db_path).mkdir(exist_ok=True)
            if assembly_id:
                (
                    path_final_genome_nt,
                    path_final_genome_aa,
                ) = get_final_fasta_with_ncbi_assemble(output_folder, assembly_id)
            if path_final_genome_nt is not None:
                print(path_final_genome_nt)
                os.system(
                    f"makeblastdb -in {path_final_genome_nt} -out "
                    f"{Path(db_path, 'nt_db')} -dbtype nucl"
                    f" -title 'nt_db' -parse_seqids"
                )
            if path_final_genome_aa is not None:
                print(path_final_genome_aa)
                os.system(
                    f"makeblastdb -in {path_final_genome_aa} -out "
                    f"{Path(db_path, 'aa_db')} -dbtype"
                    f" prot -title 'aa_db' -parse_seqids"
                )
            for model_id, model_data in self.__models.items():
                out_blast_file = Path(gene_path, model_id + "_blast.tsv")
                model_gene_file, aa_status = self.__conf[model_data["model_type"]][
                    "genome_model_strategy"
                ](
                    gene_path,
                    model_data["path_to_genome"],
                    model_data["preprocess_model"],
                    model_data["model_type"],
                    model_id,
                )
                blast_command = ""
                db_name = ""
                if aa_status and path_final_genome_aa is not None:
                    blast_command = "blastp"
                    db_name = "'aa_db'"
                elif aa_status and path_final_genome_nt is not None:
                    blast_command = "tblastn"
                    db_name = "'nt_db'"
                elif not aa_status and path_final_genome_nt is not None:
                    blast_command = "blastn"
                    db_name = "'nt_db'"
                elif not aa_status and path_final_genome_aa is not None:
                    blast_command = "blastx"
                    db_name = "'aa_db'"
                if blast_command == "" or db_name == "":
                    warnings.warn("\nWarning! Something wrong with aa/nt in files/DB")
                else:
                    os.system(
                        f"{blast_command} -query {model_gene_file} "
                        f"-db {Path(db_path, db_name)} "
                        f"-max_target_seqs 1 -outfmt '6' -out {out_blast_file}"
                    )

    def set_configuration(
        self,
        model_type: str,
        remove_b: bool,
        db_name: str,
        wo_periplasmic: bool,
        conv_strategy,
        genome_model_strategy,
        **kwargs,
    ):
        assert isinstance(conv_strategy, ConvBase)

        # TODO: add checks on conf input arg
        self.__conf[model_type] = {
            "remove_b": remove_b,
            "db_name": db_name,
            "wo_periplasmic": wo_periplasmic,
            "conv_strategy": conv_strategy,
            "genome_model_strategy": genome_model_strategy,
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
        assert model_id not in self.__models, f"model_id {model_id} already used"
        assert model_type in self.__conf, f"Missing configuration for {model_type}"

        model = load_sbml_model(path_to_model, cache, show_logs)

        # Populate the internal data
        self.__models[model_id] = {
            "original_model": deepcopy(model),
            "path_to_model": path_to_model,
            "model_type": model_type,
            "path_to_genome": path_to_genome,
        }

        # If model_type requires it, remove `_b` extensions
        if self.__conf.get(model_type).get("remove_b"):
            model = remove_b_type_exchange(model)

        self.__models[model_id]["preprocess_model"] = model

        # Record duplicated reactions
        self.__models[model_id]["duplicated_reactions"] = get_duplicated_reactions(
            model
        )
