#!/usr/bin/env python3

import re
from abc import ABC, abstractmethod

from .dbs import (
    get_BiGG_lists,
    get_kegg_m,
    get_kegg_r,
    get_mnx_m,
    get_mnx_r,
    get_old_bigg_m,
    get_old_bigg_r,
    get_seed_addit_m,
    get_seed_addit_r,
    get_seed_orig_m,
    get_seed_orig_r,
)


# TODO: this could be a simple function
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
        metabolite=True,
    ):
        # Check whether in appropriate db
        self.compartments = compartment
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

        # If metabolite add compartment to metabolite id
        if metabolite:
            self.annot_and_main = [
                i + "_" + self.compartments[0] for i in self.annot_and_main
            ]
            self.annot = [i + "_" + self.compartments[0] for i in self.annot]
            self.main = [i + "_" + self.compartments[0] for i in self.main]
            self.addit = [i + "_" + self.compartments[0] for i in self.addit]
            self.pattern = [i + "_" + self.compartments[0] for i in self.pattern]
            self.no_conv = [i + "_" + self.compartments[0] for i in self.no_conv]

        # Find 1st not empty conversion and set it as highest available
        for attr in ["annot_and_main", "annot", "main", "addit", "pattern", "no_conv"]:
            if getattr(self, attr) or attr == "no_conv":
                self.highest = getattr(self, attr)
                self.level = attr
                break

    def __repr__(self):
        return (
            "Converted class object\n"
            f"compartments: {self.compartments}\n"
            f"highest: {self.highest}\n"
            f"level: {self.level}\n"
            f"annot_and_main: {self.annot_and_main}\n"
            f"annot: {self.annot}\n"
            f"main: {self.main}\n"
            f"addit: {self.addit}\n"
            f"pattern: {self.pattern}\n"
            f"no_conv: {self.no_conv}\n"
        )

    def __str__(self):
        return (
            f"compartments: {self.compartments}\n"
            f"highest: {self.highest}\n"
            f"level: {self.level}\n"
            f"annot_and_main: {self.annot_and_main}\n"
            f"annot: {self.annot}\n"
            f"main: {self.main}\n"
            f"addit: {self.addit}\n"
            f"pattern: {self.pattern}\n"
            f"no_conv: {self.no_conv}\n"
        )


class ConvBase(ABC):
    def __init__(self, bigg_m=None, bigg_r=None):
        self.__bigg_m__ = bigg_m or get_BiGG_lists(metabolites=True)
        self.__bigg_r__ = bigg_r or get_BiGG_lists(metabolites=False)

    @abstractmethod
    def convert_metabolite(self, metabolite):
        return NotImplemented

    @abstractmethod
    def convert_reaction(self, metabolite):
        return NotImplemented

    def convert_all_metabolites(self, model):
        return {x.id: self.convert_metabolite(x) for x in model.metabolites}

    def convert_all_reactions(self, model):
        return {x.id: self.convert_reaction(x) for x in model.reactions}

    def convert_model(self, model):
        return {
            "metabolites": self.convert_all_metabolites(model),
            "reactions": self.convert_all_reactions(model),
        }


class ConvGapseq(ConvBase):
    def __init__(
        self,
        main_map_m=None,
        main_map_r=None,
        additional_table_m=None,
        additional_table_r=None,
        bigg_m=None,
        bigg_r=None,
    ):
        super().__init__(bigg_m, bigg_r)

        # TODO: checks that tables, if given, are of the appropriate format
        self.__main_map_m__ = main_map_m or get_seed_orig_m()
        self.__main_map_r__ = main_map_r or get_seed_orig_r()
        self.__addit_map_m__ = additional_table_m or get_seed_addit_m()
        self.__addit_map_r__ = additional_table_r or get_seed_addit_r()

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
            metabolite=True,
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
            metabolite=False,
        )


class ConvModelseed(ConvBase):
    def __init__(
        self,
        main_map_m=None,
        additional_table_m=None,
        main_map_r=None,
        additional_table_r=None,
        bigg_m=None,
        bigg_r=None,
    ):
        super().__init__(bigg_m, bigg_r)

        # TODO: checks that tables, if given, are of the appropriate format
        self.__main_map_m__ = main_map_m or get_seed_orig_m()
        self.__main_map_r__ = main_map_r or get_seed_orig_r()
        self.__addit_map_m__ = additional_table_m or get_seed_addit_m()
        self.__addit_map_r__ = additional_table_r or get_seed_addit_r()
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
            metabolite=True,
        )

    def convert_reaction(self, reaction):
        id_wo_comp = self.__comp_regex__.sub("", reaction.id)
        conv_main = self.__main_map_r__.get(id_wo_comp, [])
        conv_addit = self.__addit_map_r__.get(id_wo_comp, [])

        return Converted(
            check_db=self.__bigg_r__,
            compartment=[x.removesuffix("0") for x in reaction.compartments],
            main=conv_main,
            addit=conv_addit,
            metabolite=False,
        )


class ConvMetanetx(ConvBase):
    def __init__(
        self, main_map_m=None, main_map_r=None, bigg_m=None, bigg_r=None,
    ):
        super().__init__(bigg_m, bigg_r)

        # TODO: checks that tables, if given, are of the appropriate format
        self.__main_map_m__ = main_map_m or get_mnx_m()
        self.__main_map_r__ = main_map_r or get_mnx_r()
        self.__comp_regex__ = re.compile("[c@].*$")

    def convert_metabolite(self, metabolite):
        id_wo_comp = self.__comp_regex__.sub("", metabolite.id)
        conv_main = self.__main_map_m__.get(id_wo_comp, [])

        return Converted(
            check_db=self.__bigg_m__,
            compartment=[metabolite.compartment],
            main=conv_main,
            metabolite=True,
        )

    def convert_reaction(self, reaction):
        id_wo_comp = self.__comp_regex__.sub("", reaction.id)
        conv_main = self.__main_map_r__.get(id_wo_comp, [])

        return Converted(
            check_db=self.__bigg_r__,
            compartment=list(reaction.compartments),
            main=conv_main,
            metabolite=False,
        )


class ConvAgora(ConvBase):
    def __init__(
        self,
        main_map_m=None,
        additional_table_m=None,
        main_map_r=None,
        additional_table_r=None,
        bigg_m=None,
        bigg_r=None,
    ):
        super().__init__(bigg_m, bigg_r)

        # TODO: checks that tables, if given, are of the appropriate format
        self.__main_map_m__ = main_map_m or get_old_bigg_m()
        self.__main_map_r__ = main_map_r or get_old_bigg_r()
        self.__addit_map_m__ = additional_table_m or get_kegg_m()
        self.__addit_map_r__ = additional_table_r or get_kegg_r()

        self.__annot_m__ = "kegg.compound"
        self.__annot_r__ = "kegg.reaction"
        self.__comp_regex__ = re.compile("\[([ce])\]$")

    def convert_metabolite(self, metabolite):
        id_wo_comp = self.__comp_regex__.sub("", metabolite.id)
        conv_main = self.__main_map_m__.get(id_wo_comp, [])
        annot = metabolite.annotation.get(self.__annot_m__, "")
        annot_bigg = metabolite.annotation.get("bigg.metabolite", [])

        # Make sure required annotation is a list
        if isinstance(annot, str):
            annot = [annot]

        if isinstance(annot_bigg, str):
            annot_bigg = [annot_bigg]

        if not isinstance(annot, list):
            raise ValueError(
                f"{self.__annot_m__} annotation of {metabolite.id} "
                f"metabolite is written in wrong type: {type(annot)}. "
                "It has to be either str or list. Please change it in "
                "the model or remove completely."
            )

        conv_addit = [y for x in annot for y in self.__addit_map_m__.get(x, [])]
        conv_pattern = ["__".join(id_wo_comp.rsplit("_", 1))]
        if annot_bigg:
            conv_for_annot = annot_bigg
        else:
            conv_for_annot = [id_wo_comp]

        return Converted(
            check_db=self.__bigg_m__,
            compartment=[metabolite.compartment],
            main=conv_main,
            addit=conv_addit,
            pattern=conv_pattern,
            annot=conv_for_annot,
            metabolite=True,
        )

    def convert_reaction(self, reaction):
        id_wo_comp = self.__comp_regex__.sub("", reaction.id)
        conv_main = self.__main_map_r__.get(id_wo_comp, [])
        conv_addit = [
            y
            for x in [reaction.annotation.get(self.__annot_r__, "")]
            for y in self.__addit_map_r__.get(x, [])
        ]
        conv_pattern = ["__".join(id_wo_comp.rsplit("_", 1))]
        conv_noconv = [id_wo_comp]

        return Converted(
            check_db=self.__bigg_r__,
            compartment=list(reaction.compartments),
            main=conv_main,
            addit=conv_addit,
            pattern=conv_pattern,
            annot=conv_noconv,
            metabolite=False,
        )


class ConvCarveme(ConvBase):
    def __init__(self, main_map_m=None, main_map_r=None, bigg_m=None, bigg_r=None):
        super().__init__(bigg_m, bigg_r)

        # TODO: checks that tables, if given, are of the appropriate format
        self.__main_map_m__ = main_map_m or get_old_bigg_m()
        self.__main_map_r__ = main_map_r or get_old_bigg_r()
        self.__comp_regex__ = re.compile("_([cep])$")

    def convert_metabolite(self, metabolite):
        id_wo_comp = self.__comp_regex__.sub("", metabolite.id)
        if id_wo_comp not in self.__bigg_m__:
            conv_main = [
                x + metabolite.compartment.split("_")[1]
                for x in self.__main_map_m__.get(id_wo_comp, [])
            ]
        else:
            conv_main = [id_wo_comp]

        comp = metabolite.compartment.split("_")
        if len(comp) == 1:
            raise ValueError(
                f"metabolite '{metabolite.id}' has incompatible "
                f"compartment: {metabolite.compartment}"
            )

        return Converted(
            check_db=self.__bigg_m__,
            compartment=[comp[1]],
            main=conv_main,
            metabolite=True,
        )

    def convert_reaction(self, reaction):
        id_wo_comp = reaction.id
        if id_wo_comp not in self.__bigg_r__:
            conv_main = [x for x in self.__main_map_r__.get(id_wo_comp, [])]
        else:
            conv_main = [id_wo_comp]

        return Converted(
            check_db=self.__bigg_r__,
            compartment=[x.split("_")[1] for x in reaction.compartments],
            main=conv_main,
            metabolite=False,
        )


class ConvBigg(ConvBase):
    def __init__(self, main_map_m=None, main_map_r=None, bigg_m=None, bigg_r=None):
        super().__init__(bigg_m, bigg_r)

        # TODO: checks that tables, if given, are of the appropriate format
        self.__main_map_m__ = main_map_m or get_old_bigg_m()
        self.__main_map_r__ = main_map_r or get_old_bigg_r()
        self.__comp_regex__ = re.compile("_([cep])$")

    def convert_metabolite(self, metabolite):
        id_wo_comp = self.__comp_regex__.sub("", metabolite.id)
        if id_wo_comp not in self.__bigg_m__:
            conv_main = [
                x + metabolite.compartment
                for x in self.__main_map_m__.get(id_wo_comp, [])
            ]
        else:
            conv_main = [id_wo_comp]

        return Converted(
            check_db=self.__bigg_m__,
            compartment=[metabolite.compartment],
            main=conv_main,
            metabolite=True,
        )

    def convert_reaction(self, reaction):
        id_wo_comp = reaction.id
        if id_wo_comp not in self.__bigg_r__:
            conv_main = [x for x in self.__main_map_r__.get(id_wo_comp, [])]
        else:
            conv_main = [id_wo_comp]

        return Converted(
            check_db=self.__bigg_r__,
            compartment=list(reaction.compartments),
            main=conv_main,
            metabolite=False,
        )


def remove_zero_for_notconv(notconv_id: str):
    return re.sub("0$", "", notconv_id)


def replace_square_brackets(notconv_id: str):
    if (notconv_id[-3] == "[") and (notconv_id[-1] == "]"):
        return notconv_id[:-3] + "_" + notconv_id[-2]
    else:
        return notconv_id


def no_changes_for_notconv(notconv_id: str):
    return notconv_id
