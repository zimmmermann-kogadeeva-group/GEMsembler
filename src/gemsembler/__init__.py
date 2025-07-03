from importlib.resources import files

from . import data
from .anticreation import get_model_of_interest, get_models_with_all_confidence_levels
from .creation import read_supermodel_from_json
from .gathering import GatheredModels, load_sbml_model

__version__ = "0.10.3"

lp_example = [
    dict(
        model_id="agora_LP",
        path_to_model=files(data) / "agora_model.xml.gz",
        model_type="agora",
        path_to_genome=files(data) / "genome_agora.fna",
    ),
    dict(
        model_id="carveme_LP",
        path_to_model=files(data) / "carveme_model.xml.gz",
        model_type="carveme",
        path_to_genome=files(data) / "genome.faa",
    ),
    dict(
        model_id="gapseq_LP",
        path_to_model=files(data) / "gapseq_model.xml.gz",
        model_type="gapseq",
        path_to_genome=files(data) / "genome.fna",
    ),
    dict(
        model_id="modelseed_LP",
        path_to_model=files(data) / "modelseed_model.xml.gz",
        model_type="modelseed",
        path_to_genome=files(data) / "genome.faa",
    ),
]
