import json
from importlib.resources import files

import pandas as pd

from . import data
from .anticreation import get_model_of_interest, get_models_with_all_confidence_levels
from .creation import read_supermodel_from_json
from .gathering import GatheredModels, load_sbml_model

__version__ = "0.10.5"


with open(files(data) / "masses_and_charges.json") as fh:
    df_bigg_extra = (
        pd.DataFrame(
            [
                (bigg_id, *x)
                for bigg_id, info in json.load(fh).items() if info is not None
                for x in zip(info["formula"], info["mass"], info["charges"])
            ],
            columns=["bigg_id", "formula", "mass", "charges"]
        )
        .drop_duplicates(subset="bigg_id", keep="first")
    )

lp_example = [
    dict(
        model_id="agora_LP",
        path_to_model=files(data) / "LP" / "LP_agora_model.xml.gz",
        model_type="agora",
        path_to_genome=files(data) / "LP" / "LP_genome_agora.fna.gz",
    ),
    dict(
        model_id="carveme_LP",
        path_to_model=files(data) / "LP" / "LP_carveme_model.xml.gz",
        model_type="carveme",
        path_to_genome=files(data) / "LP" / "LP_genome.faa.gz",
    ),
    dict(
        model_id="gapseq_LP",
        path_to_model=files(data) / "LP" / "LP_gapseq_model.xml.gz",
        model_type="gapseq",
        path_to_genome=files(data) / "LP" / "LP_genome.fna.gz",
    ),
    dict(
        model_id="modelseed_LP",
        path_to_model=files(data) / "LP" / "LP_modelseed_model.xml.gz",
        model_type="modelseed",
        path_to_genome=files(data) / "LP" / "LP_genome.faa.gz",
    ),
]

bu_example = [
    dict(
        model_id="carveme_BU",
        path_to_model=files(data) / "BU" / "BU_carveme_hom.xml.gz",
        model_type="carveme",
    ),
    dict(
        model_id="agora_BU",
        path_to_model=files(data) / "BU" / "BU_agora.xml.gz",
        model_type="agora",
    ),
    dict(
        model_id="modelseed_BU",
        path_to_model=files(data) / "BU" / "BU_modelSEED.sbml.gz",
        model_type="modelseed",
    ),
    dict(
        model_id="gapseq_BU",
        path_to_model=files(data) / "BU" / "BU_gapseq.xml.gz",
        model_type="gapseq",
    ),
]
