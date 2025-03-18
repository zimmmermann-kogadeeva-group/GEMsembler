from importlib.resources import files
from .anticreation import get_model_of_interest, get_models_with_all_confidence_levels
from .creation import read_supermodel_from_json
from .data import BU, LP
from .gathering import GatheredModels, load_sbml_model

__version__ = "0.9.9"

lp_example = [
    dict(
        model_id="curated_LP",
        path_to_model=files(LP) / "LP_iLP728_revision_data_met_C_c.xml.gz",
        model_type="carveme",
        path_to_genome=files(LP) / "LP_protein_fasta.faa.gz",
    ),
    dict(
        model_id="cauniv_LP",
        path_to_model=files(LP) / "LP_CA1.xml.gz",
        model_type="carveme",
        path_to_genome=files(LP) / "LP_protein_fasta.faa.gz",
    ),
    dict(
        model_id="cagram_LP",
        path_to_model=files(LP) / "LP_CA2.xml.gz",
        model_type="carveme",
        path_to_genome=files(LP) / "LP_protein_fasta.faa.gz",
    ),
    dict(
        model_id="msgram_LP",
        path_to_model=files(LP) / "LP_MS2.sbml.gz",
        model_type="modelseed",
        path_to_genome=files(LP) / "LP_protein_fasta.faa.gz",
    ),
    dict(
        model_id="agora_LP",
        path_to_model=files(LP) / "LP_WCFS1_agora.xml.gz",
        model_type="agora",
        path_to_genome=files(LP) / "LP_WCFS1.fasta.gz",
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
