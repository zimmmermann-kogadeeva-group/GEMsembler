#!/usr/bin/env python3

import pandas as pd
import re

from .dbs import (
    seed_orig_m,
    seed_orig_r,
    seed_addit_m,
    seed_addit_r,
    kegg_bigg_m,
    kegg_bigg_r,
    old_new_bigg_m,
    old_new_bigg_r,
    bigg_all_m,
    bigg_all_r,
    bigg_db_network_m,
    bigg_db_network_r,
)


class Converted(object):
    def __init__(
        self,
        check_db,
        compartment,
        annot=None,
        main=None,
        addit=None,
        pattern=None,
        no_conv=None,
    ):
        # Check whether in appropriate db
        self.compartment = compartment
        annot = set() if annot is None else {x for x in annot if x in check_db}
        main = set() if main is None else {x for x in main if x in check_db}
        addit = set() if addit is None else {x for x in addit if x in check_db}
        pattern = set() if pattern is None else {x for x in pattern if x in check_db}
        no_conv = set() if no_conv is None else {x for x in no_conv if x in check_db}

        # Do prioritisation
        self.annot_and_main = list(annot.intersection(main))
        self.annot = list(annot - main)
        self.main = list(main - annot)
        self.addit = list(addit - annot - main)
        self.pattern = list(pattern - annot - main - addit)
        self.no_conv = list(no_conv - annot - main - addit - pattern)

        # TODO: add back compartments

    def __repr__(self):
        return (
            "Converted class object\n"
            f"compartment: {self.compartment}\n"
            f"annot_and_main: {self.annot_and_main}\n"
            f"annot: {self.annot}\n"
            f"main: {self.main}\n"
            f"addit: {self.addit}\n"
            f"pattern: {self.pattern}\n"
            f"no_conv: {self.no_conv}\n"
        )

    def __str__(self):
        return (
            f"compartment: {self.compartment}\n"
            f"annot_and_main: {self.annot_and_main}\n"
            f"annot: {self.annot}\n"
            f"main: {self.main}\n"
            f"addit: {self.addit}\n"
            f"pattern: {self.pattern}\n"
            f"no_conv: {self.no_conv}\n"
        )


class ConvBase(object):
    def __init__(self, bigg_m=None, bigg_r=None):
        self.__bigg_m__ = bigg_m or set(bigg_all_m["universal_bigg_id"])
        self.__bigg_r__ = bigg_r or set(bigg_all_r["universal_bigg_id"])

    def convert_model(self, model):
        return {
            "metabolites": [self.convert_metabolite(x) for x in model.metabolites],
            "reactions": [self.convert_reaction(x) for x in model.reactions],
        }


class ConvGapseq(ConvBase):
    def __init__(
        self, main_map_m=None, main_map_r=None, additional_table_m=None,
        additional_table_r=None, bigg_m=None, bigg_r=None
    ):
        super().__init__(bigg_m, bigg_r)

        # TODO: checks that tables, if given, are of the appropriate format
        self.__main_map_m__ = main_map_m or seed_orig_m
        self.__main_map_r__ = main_map_r or seed_orig_r
        self.__addit_map_m__ = additional_table_m or seed_addit_m
        self.__addit_map_r__ = additional_table_r or seed_addit_r

        self.__annot_m__ = "bigg.metabolite"
        self.__annot_r__ = "bigg.reaction"
        self.__comp_regex__ = re.compile("_(c0|e0|p0)$")

    def convert_metabolite(self, metabolite):
        conv_annot = metabolite.annotation.get(self.__annot_m__, [])
        conv_annot = [conv_annot] if type(conv_annot) is str else conv_annot
        id_wo_comp = self.__comp_regex__.sub("", metabolite.id)
        conv_main = self.__main_map_m__.get(id_wo_comp, [])
        conv_addit = self.__addit_map_m__.get(id_wo_comp, [])

        return Converted(
            check_db=self.__bigg_m__,
            compartment=[metabolite.compartment.removesuffix("0")],
            annot=conv_annot,
            main=conv_main,
            addit=conv_addit,
        )

    def convert_reaction(self, reaction):
        conv_annot = reaction.annotation.get(self.__annot_r__, [])
        conv_annot = [conv_annot] if type(conv_annot) is str else conv_annot
        id_wo_comp = self.__comp_regex__.sub("", reaction.id)
        conv_main = self.__main_map_r__.get(id_wo_comp, [])
        conv_addit = self.__addit_map_r__.get(id_wo_comp, [])

        return Converted(
            check_db=self.__bigg_r__,
            compartment=[x.removesuffix("0") for x in reaction.compartments],
            annot=conv_annot,
            main=conv_main,
            addit=conv_addit,
        )


class ConvModelseed(ConvBase):
    def __init__(
        self, main_map_m=None, additional_table_m=None, main_map_r=None,
        additional_table_r=None, bigg_m=None, bigg_r=None
    ):
        super().__init__(bigg_m, bigg_r)

        # TODO: checks that tables, if given, are of the appropriate format
        self.__main_map_m__ = main_map_m or seed_orig_m
        self.__main_map_r__ = main_map_r or seed_orig_r
        self.__addit_map_m__ = additional_table_m or seed_addit_m
        self.__addit_map_m__ = additional_table_r or seed_addit_r
        self.__comp_regex__ = re.compile("_(c0|e0|b)$")

    def convert_metabolite(self, metabolite):
        id_wo_comp = self.__comp_regex__.sub("", metabolite.id)
        conv_main = self.__main_map_m__.get(id_wo_comp, [])
        conv_addit = self.__addit_map_m__.get(id_wo_comp, [])

        return Converted(
            check_db=self.__bigg_m__,
            compartment=[metabolite.compartment.removesuffix("0")],
            main=conv_main,
            addit=conv_addit,
        )

    def convert_reaction(self, reaction):
        id_wo_comp = self.__comp_regex__.sub("", reaction.id)
        conv_main = self.__main_map_r__.get(id_wo_comp, [])
        conv_addit = self.__addit_map_r__.get(id_wo_comp, [])

        return Converted(
            check_db=self.__bigg_r__,
            compartment=[x.removesuffix("0") for x in reaction.compartments],
            annot=conv_annot,
            main=conv_main,
            addit=conv_addit,
        )


class ConvAgora(ConvBase):
    def __init__(
        self, main_map_m=None, additional_table_m=None, main_map_r=None,
        additional_table_r=None, bigg_m=None, bigg_r=None
    ):
        super().__init__(bigg_m, bigg_r)

        # TODO: checks that tables, if given, are of the appropriate format
        self.__main_map_m__ = main_map_m or old_new_bigg_m
        self.__main_map_r__ = main_map_r or old_new_bigg_r
        self.__addit_map_m__ = additional_table_m or kegg_bigg_m
        self.__addit_map_r__ = additional_table_r or kegg_bigg_r

        self.__annot_m__ = "kegg.compound"
        self.__annot_r__ = "kegg.reaction"
        self.__comp_regex__ = re.compile("\[(c|e)\]$")

    def convert_metabolite(self, metabolite):
        id_wo_comp = self.__comp_regex__.sub("", metabolite.id)
        conv_main = self.__main_map_m__.get(id_wo_comp, [])
        conv_addit = [
            y
            for x in metabolite.annotation.get(self.__annot_m__, [])
            for y in self.__addit_map_m__.get(x, [])
        ]
        conv_pattern = ["_".join(id_wo_comp.rsplit("__", 1))]
        conv_noconv = [id_wo_comp]

        return Converted(
            check_db=self.__bigg_m__,
            compartment=[metabolite.compartment],
            main=conv_main,
            addit=conv_addit,
            pattern=conv_pattern,
            no_conv=conv_noconv,
        )

    def convert_reaction(self, reaction):
        id_wo_comp = self.__comp_regex__.sub("", reaction.id)
        conv_main = self.__main_map_r__.get(id_wo_comp, [])
        conv_addit = [
            y
            for x in reaction.annotation.get(self.__annot_r__, [])
            for y in self.__addit_map_r__.get(x, [])
        ]
        conv_pattern = ["_".join(id_wo_comp.rsplit("__", 1))]
        conv_noconv = [id_wo_comp]

        return Converted(
            check_db=self.__bigg_m__,
            compartment=reaction.compartments,
            main=conv_main,
            addit=conv_addit,
            pattern=conv_pattern,
            no_conv=conv_noconv,
        )


class ConvCarveme(ConvBase):
    def __init__(self, main_map_m=None, main_map_r=None, bigg_m=None, bigg_r=None):
        super().__init__(bigg_m, bigg_r)

        # TODO: checks that tables, if given, are of the appropriate format
        self.__main_map_m__ = main_map_m or old_new_bigg_m
        self.__main_map_r__ = main_map_r or old_new_bigg_r
        self.__comp_regex__ = re.compile("_(c|e|p)$")

    def convert_metabolite(self, metabolite):
        id_wo_comp = self.__comp_regex__.sub("", metabolite.id)
        if id_wo_comp not in self.__bigg_m__:
            conv_main = [
                x + metabolite.compartment.split("_")[1]
                for x in self.__main_map_m__.get(id_wo_comp, [])
            ]
        else:
            conv_main = [id_wo_comp]

        return Converted(
            check_db=self.__bigg_m__,
            compartment=[metabolite.compartment.split("_")[1]],
            main=conv_main,
        )

    def convert_reaction(self, reaction):
        id_wo_comp = reaction.id
        if id_wo_comp not in self.__bigg_r__:
            conv_main = [
                x + reaction.compartment.split("_")[1]
                for x in self.__main_map_r__.get(id_wo_comp, [])
            ]
        else:
            conv_main = [id_wo_comp]

        return Converted(
            check_db=self.__bigg_r__,
            compartment=[x.split("_")[1] for x in reaction.compartments],
            main=conv_main,
        )
