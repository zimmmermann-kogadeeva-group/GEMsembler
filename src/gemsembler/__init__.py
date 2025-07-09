from importlib.resources import files

from .anticreation import get_model_of_interest, get_models_with_all_confidence_levels
from .creation import read_supermodel_from_json
from .data import BU, LP
from .gathering import GatheredModels, load_sbml_model

__version__ = "0.10.5"

lp_example = [
    dict(
        model_id="agora_LP",
        path_to_model=files(LP) / "LP_agora_model.xml.gz",
        model_type="agora",
        path_to_genome=files(LP) / "LP_genome_agora.fna",
    ),
    dict(
        model_id="carveme_LP",
        path_to_model=files(LP) / "LP_carveme_model.xml.gz",
        model_type="carveme",
        path_to_genome=files(LP) / "LP_genome.faa",
    ),
    dict(
        model_id="gapseq_LP",
        path_to_model=files(LP) / "LP_gapseq_model.xml.gz",
        model_type="gapseq",
        path_to_genome=files(LP) / "LP_genome.fna",
    ),
    dict(
        model_id="modelseed_LP",
        path_to_model=files(LP) / "LP_modelseed_model.xml.gz",
        model_type="modelseed",
        path_to_genome=files(LP) / "LP_genome.faa",
    ),
]

bu_example = [
    dict(
        model_id="carveme_BU",
        path_to_model=files(BU) / "BU_carveme_hom.xml.gz",
        model_type="carveme",
    ),
    dict(
        model_id="agora_BU",
        path_to_model=files(BU) / "BU_agora.xml.gz",
        model_type="agora",
    ),
    dict(
        model_id="modelseed_BU",
        path_to_model=files(BU) / "BU_modelSEED.sbml.gz",
        model_type="modelseed",
    ),
    dict(
        model_id="gapseq_BU",
        path_to_model=files(BU) / "BU_gapseq.xml.gz",
        model_type="gapseq",
    ),
]
