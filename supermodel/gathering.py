import pandas as pd
from BiGGnetwork import get_dbs
import conversion


class StrategiesForModelType:
    """Whole strategy of processing particular model type"""

    def __init__(self, path_to_db=None):
        db = get_dbs(path_to_db)
        bigg_m = list(set(db.get("bigg_all_m")["universal_bigg_id"]))
        bigg_r = list(set(db.get("bigg_all_r")["universal_bigg_id"]))
        self.model_type_list = ["carveme", "gapseq", "modelseed", "agora"]
        self.remove_b = {"carveme": False, "gapseq": False, "modelseed": True, "agora": False}
        self.db_name = {"carveme": "bigg", "gapseq": "modelseed", "modelseed": "modelseed", "agora": "wierd_bigg"}
        self.wo_periplasmic = {"carveme": False, "gapseq": False, "modelseed": True, "agora": True}
        self.conversion_strategies = {"carveme": conversion.ConversionForCarveMe(db["old_new_bigg_m"],
                                                                                 db["old_new_bigg_r"], bigg_m, bigg_r),
                                      "gapseq": conversion.ConversionForGapseq(db["seed_orig_m"], db["seed_orig_r"],
                                                                db["seed_addit_m"], db["seed_addit_r"], bigg_m, bigg_r),
                                      "modelseed": conversion.ConversionForModelseed(db["seed_orig_m"], db["seed_orig_r"],
                                                                db["seed_addit_m"], db["seed_addit_r"], bigg_m, bigg_r),
                                      "agora": conversion.ConversionForAgora(db["old_new_bigg_m"], db["old_new_bigg_r"],
                                                                db["kegg_bigg_m"], db["kegg_bigg_r"], bigg_m, bigg_r)}
