from importlib.resources import files

from cobra.core.model import Model

from gemsembler import GatheredModels
from gemsembler.conversion import (
    ConvAgora,
    ConvBase,
    ConvCarveme,
    ConvGapseq,
    ConvModelseed,
)
from gemsembler.data import BU


class TestGathering:
    def test_gathered_models(self):
        # Test GatheredModels init method
        g = GatheredModels()

        # Check configuration in GatheredModels object
        conf = g.get_conf()

        assert "agora" in conf
        assert conf.get("agora").get("remove_b") == False
        assert conf.get("agora").get("db_name") == "weird_bigg"
        assert conf.get("agora").get("wo_periplasmic") == True
        assert type(conf.get("agora").get("conv_strategy")) is ConvAgora

        assert "carveme" in conf
        assert conf.get("carveme").get("remove_b") == False
        assert conf.get("carveme").get("db_name") == "bigg"
        assert conf.get("carveme").get("wo_periplasmic") == False
        assert type(conf.get("carveme").get("conv_strategy")) is ConvCarveme

        assert "gapseq" in conf
        assert conf.get("gapseq").get("remove_b") == False
        assert conf.get("gapseq").get("db_name") == "modelseed"
        assert conf.get("gapseq").get("wo_periplasmic") == False
        assert type(conf.get("gapseq").get("conv_strategy")) is ConvGapseq

        assert "modelseed" in conf
        assert conf.get("modelseed").get("remove_b") == True
        assert conf.get("modelseed").get("db_name") == "modelseed"
        assert conf.get("modelseed").get("wo_periplasmic") == True
        assert type(conf.get("modelseed").get("conv_strategy")) is ConvModelseed

        # Check models in GatheredModels object
        assert len(g) == 0

        g.add_model("test_carveme", files(BU) / "BU_carveme_hom.xml.gz", "carveme")
        assert "test_carveme" in g
        model_attrs = g.get_model_attrs("test_carveme")
        assert model_attrs["path_to_model"] == files(BU) / "BU_carveme_hom.xml.gz"
        assert model_attrs["model_type"] == "carveme"
        model = model_attrs["preprocess_model"]
        assert type(model) is Model
        assert len(model.metabolites) == 1244
        assert len(model.reactions) == 1849

        g.add_model("test_agora", files(BU) / "BU_agora.xml.gz", "agora")
        assert "test_agora" in g
        model_attrs = g.get_model_attrs("test_agora")
        assert model_attrs["path_to_model"] == files(BU) / "BU_agora.xml.gz"
        assert model_attrs["model_type"] == "agora"
        model = model_attrs["preprocess_model"]
        assert type(model) is Model
        assert len(model.metabolites) == 1500
        assert len(model.reactions) == 2418

        g.add_model("test_modelseed", files(BU) / "BU_modelSEED.sbml.gz", "modelseed")
        assert "test_modelseed" in g
        model_attrs = g.get_model_attrs("test_modelseed")
        assert model_attrs["path_to_model"] == files(BU) / "BU_modelSEED.sbml.gz"
        assert model_attrs["model_type"] == "modelseed"
        model = model_attrs["preprocess_model"]
        assert type(model) is Model
        assert len(model.metabolites) == 1271
        assert len(model.reactions) == 1161

        g.add_model("test_gapseq", files(BU) / "BU_gapseq.xml.gz", "gapseq")
        assert "test_gapseq" in g
        model_attrs = g.get_model_attrs("test_gapseq")
        assert model_attrs["path_to_model"] == files(BU) / "BU_gapseq.xml.gz"
        assert model_attrs["model_type"] == "gapseq"
        model = model_attrs["preprocess_model"]
        assert type(model) is Model
        assert len(model.metabolites) == 1520
        assert len(model.reactions) == 1891

        assert g.same_db_models["bigg"] == {"test_carveme": "carveme"}
        assert g.same_db_models["weird_bigg"] == {"test_agora": "agora"}
        assert g.same_db_models["modelseed"] == {
            "test_modelseed": "modelseed",
            "test_gapseq": "gapseq",
        }
        # TODO: finish the tests for gathered models
        g.run()
