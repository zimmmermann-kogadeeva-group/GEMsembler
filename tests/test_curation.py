from supermodel import load_sbml_model
from supermodel.curation import remove_b_type_exchange, get_duplicated_reactions


class TestCuration:
    def test_remove_b_type_exchange(self):
        # Load the example modelseed model. Only modelseed models require
        # removal of `b` compartments from metabolites and reactions
        model = load_sbml_model("example/BU_modelSEED.sbml.gz")

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
        model = load_sbml_model("example/BU_agora.xml.gz")
        struct_dup, gpr_dup = get_duplicated_reactions(model)
        assert len(struct_dup) == 6
        assert len(gpr_dup) == 0

        # Check sizes of duplicated tables given carveme model
        model = load_sbml_model("example/BU_carveme_hom.xml.gz")
        struct_dup, gpr_dup = get_duplicated_reactions(model)
        assert len(struct_dup) == 96
        assert len(gpr_dup) == 32

        # Check sizes of duplicated tables given gapseq model
        model = load_sbml_model("example/BU_gapseq.xml.gz")
        struct_dup, gpr_dup = get_duplicated_reactions(model)
        assert len(struct_dup) == 78
        assert len(gpr_dup) == 16

        # Check sizes of duplicated tables given modelseed model
        model = load_sbml_model("example/BU_modelSEED.sbml.gz")
        struct_dup, gpr_dup = get_duplicated_reactions(model)
        assert len(struct_dup) == 0
        assert len(gpr_dup) == 0
