import pandas as pd
from .dbs import get_dbs
import conversion


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
            "carveme": conversion.ConvCarveme(),
            "gapseq": conversion.ConvGapseq(),
            "modelseed": conversion.ConvModelseed(),
            "agora": conversion.ConvAgora(),
        }
