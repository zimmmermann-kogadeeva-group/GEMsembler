import sys
import operator
from collections import Counter
from copy import deepcopy
from .conversion import ConvCarveme, ConvGapseq, ConvModelseed, ConvAgora
from cobra.io import read_sbml_model
from .curation import remove_b_type_exchange, get_duplicated_reactions
from .general import findKeysByValue
from .selection import checkDBConsistency, checkFromOneFromMany
from .structural import runStructuralConversion, runStructuralCheck
from .dbs import get_bigg_network


class StrategiesForModelType:
    """Whole strategy of processing particular model type"""

    def __init__(self, path_to_db=None):
        self.model_type_list = ["carveme", "gapseq", "modelseed", "agora"]
        self.remove_b = {
            "carveme": False,
            "gapseq": False,
            "modelseed": True,
            "agora": False,
        }
        self.db_name = {
            "carveme": "bigg",
            "gapseq": "modelseed",
            "modelseed": "modelseed",
            "agora": "wierd_bigg",
        }
        self.wo_periplasmic = {
            "carveme": False,
            "gapseq": False,
            "modelseed": True,
            "agora": True,
        }
        self.conversion_strategies = {
            "carveme": ConvCarveme(),
            "gapseq": ConvGapseq(),
            "modelseed": ConvModelseed(),
            "agora": ConvAgora(),
        }


class GatheredModels:
    """Class, that gathers information and necessary conversion results for all models. Input for the class and
    tool in general is dictionary with all models and related information.
    This dictionary dict_of_all_models_with_feature:

    {model_id:
    {'path_to_model':str,
    'model_type':str one of (agora, carveme, gapseq, modelseed) or if custom, create type class in advance,
    'path_to_genome': str (can be '' or None if convert_genes = False)}}

    And other parameters:
    if all 3 parameters bellow are None then gene conversion is not done and genomes for model_id are not need
    assembly = None
    path_final_genome_nt = None
    path_final_genome_aa = None

    Optional: custom_model_type """

    def __init__(
        self,
        dict_of_all_models_with_feature: dict,
        assembly_id=None,
        path_final_genome_nt=None,
        path_final_genome_aa=None,
        path_to_db=None,
        custom_model_type=None,
    ):
        model_ids_checking = Counter(list(dict_of_all_models_with_feature.keys()))
        not_uniq_ids = findKeysByValue(model_ids_checking, 1, operator.gt)
        if len(not_uniq_ids) >= 1:
            sys.exit(f"Some model ids are not unique: {' '.join(not_uniq_ids)}")
        if not assembly_id and not path_final_genome_nt and not path_final_genome_aa:
            print(
                "Warning! No final genome for gene conversion is provided. Gene conversion will not be performed.\n"
                "If you want to convert genes, please provide either assembly id or custom fasta files (nt/aa/both), to wich genes must be converted."
            )
            convert_genes = False
        else:
            convert_genes = True
        strategies = StrategiesForModelType(path_to_db)
        bigg_network = get_bigg_network()
        models_same_db = {db: {} for db in strategies.db_name.values()}
        self.original_models = {}
        self.preprocessed_models = {}
        self.duplicated_r = {}
        self.converted = {}
        for k, v in dict_of_all_models_with_feature.items():
            model = read_sbml_model(v["path_to_model"])
            models_same_db[strategies.db_name[v["model_type"]]].update(
                {k: v["model_type"]}
            )
            self.original_models.update({k: model})
            if strategies.remove_b[v["model_type"]]:
                model_b_removed = remove_b_type_exchange(deepcopy(model))
                self.preprocessed_models.update({k: model_b_removed})
            else:
                self.preprocessed_models.update({k: model})
            dupl_r, dupl_r_gpr = get_duplicated_reactions(model)
            self.duplicated_r.update({k: dupl_r})
            self.converted.update(
                {
                    k: strategies.conversion_strategies[v["model_type"]].convert_model(
                        self.preprocessed_models[k]
                    )
                }
            )
        self.first_selected = checkDBConsistency(
            models_same_db, self.converted, "highest"
        )
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
