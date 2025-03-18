#!/usr/bin/env python3


def _get_diffs(model1, model2, attribute):
    if hasattr(model1, attribute) and hasattr(model2, attribute):
        model1_elems = {x.id for x in getattr(model1, attribute)}
        model2_elems = {x.id for x in getattr(model2, attribute)}
        return {
            "model1 - model2": model1_elems - model2_elems,
            "model2 - model1": model2_elems - model1_elems,
        }


def quickdiff(model1, model2):
    attrs = ["metabolites", "reactions", "genes"]

    diffs = {attr: _get_diffs(model1, model2, attr) for attr in attrs}

    if all(
        [
            len(x["model1 - model2"]) == 0 and len(x["model2 - model1"]) == 0
            for attr, x in diffs.items()
        ]
    ):
        return True

    return diffs
