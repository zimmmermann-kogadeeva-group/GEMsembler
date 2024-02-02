from importlib.resources import files

from gemsembler import load_sbml_model
from gemsembler.curation import get_duplicated_reactions, remove_b_type_exchange
from gemsembler.data import BU


class TestCuration:
    def test_remove_b_type_exchange(self):
        # Load the example modelseed model. Only modelseed models require
        # removal of `b` compartments from metabolites and reactions
        model = load_sbml_model(files(BU) / "BU_modelSEED.sbml.gz")

        # Get and check the numer of metabolites and reactions before removing
        # metabolites and reactions in the `b` compartment
        num_b_mbs = len([m for m in model.metabolites if m.id.endswith("_b")])
        assert num_b_mbs == 61

        num_b_reacs = len([r for r in model.reactions if r.id.endswith("_b")])
        assert num_b_reacs == 61

        # Remove the metabolites and reactions in the `b` compartment
        model_wo_bs = remove_b_type_exchange(model)

        # Get and check the numer of metabolites and reactions after removing
        # metabolites and reactions in the `b` compartment
        num_b_mbs = len([m for m in model_wo_bs.metabolites if m.id.endswith("_b")])
        assert num_b_mbs == 0

        num_b_reacs = len([r for r in model_wo_bs.reactions if r.id.endswith("_b")])
        assert num_b_reacs == 0

    def test_get_duplicated_reactions(self):
        # Check sizes of duplicated tables given agora model
        model = load_sbml_model(files(BU) / "BU_agora.xml.gz")
        dup = get_duplicated_reactions(model)
        assert len(dup) == 20

        # Check sizes of duplicated tables given carveme model
        model = load_sbml_model(files(BU) / "BU_carveme_hom.xml.gz")
        dup = get_duplicated_reactions(model)
        assert len(dup) == 204

        # Check sizes of duplicated tables given gapseq model
        model = load_sbml_model(files(BU) / "BU_gapseq.xml.gz")
        dup = get_duplicated_reactions(model)
        assert len(dup) == 110

        # Check sizes of duplicated tables given modelseed model
        model = load_sbml_model(files(BU) / "BU_modelSEED.sbml.gz")
        dup = get_duplicated_reactions(model)
        assert len(dup) == 0
