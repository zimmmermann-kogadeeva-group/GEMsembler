from abc import ABC, abstractmethod
import cobra
import pandas as pd


def intersection(a, b):
    out = []
    for el in a:
        if (el in b) & (el not in out):
            out.append(el)
    return out

def substraction(a, *args):
    out = []
    for el in a:
        for b in args:
            if (el not in b) & (el not in out):
                out.append(el)
    return out

# region Working with Compartments for different models
class Compartments(ABC):

    def getIDwoCompartment(self, original_id: str, pattern_to_remove: [str]) -> str:
        id_wo_compartment = original_id
        for pattern in pattern_to_remove:
            id_wo_compartment = id_wo_compartment.removesuffix(pattern)
        return id_wo_compartment

    @abstractmethod
    def getCompartment(self, original_compartment: str) -> str:
        """get compartments in standardised form"""

    def runGetCompartment(self, original: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction) -> [str]:
        if type(original) == cobra.core.metabolite.Metabolite:
            standard_compartments = [self.getCompartment(original.compartment)]
        if type(original) == cobra.core.reaction.Reaction:
            standard_compartments = []
            for comp in original.compartments:
                standard_compartments.append(self.getCompartment(comp))
        return standard_compartments


class CompartmentsForCarveme(Compartments):

    def getCompartment(self, original_compartment: str) -> str:
        standard_compartment = original_compartment.split("_")[1]
        return standard_compartment


class CompartmentsForGapseq(Compartments):

    def getCompartment(self, original_compartment: str) -> str:
        standard_compartment = original_compartment[:-1]
        return standard_compartment


class CompartmentsForModelseed(Compartments):

    def getCompartment(self, original_compartment: str) -> str:
        standard_compartment = original_compartment[:-1]
        return standard_compartment


class CompartmentsForAgora(Compartments):

    def getCompartment(self, original_compartment: str) -> str:
        standard_compartment = original_compartment
        return standard_compartment


# endregion

# region Working on converting ids for different models from different sources
class ConversionToBiGG(ABC):
    """Conversion IDs for metabolites and reactions to BiGG"""

    def __init__(self, met_bigg_to_check: list, react_bigg_to_check: list):
        self.met_bigg_to_check = met_bigg_to_check
        self.react_bigg_to_check = react_bigg_to_check

    def convertViaTable(self, id_wo_compartment: str, conversion_table: pd.core.frame.DataFrame) -> [str]:
        converted = conversion_table[conversion_table["old"] == id_wo_compartment]["new"]
        if converted.empty or converted.isnull().values.any():
            return []
        else:
            return converted.values[0].split(",")

    def convertViaAnnotation(self, annotation: dict, type_annotation: str) -> [str]:
        converted = annotation.get(type_annotation)
        if converted is not None:
            if type(converted) == list:
                return converted
            else:
                return [converted]
        else:
            return []

    def checkTMPinBigg(self, type_object: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction,
                       tmp_converted: dict) -> [str]:
        converted_bigg_ids = {}
        if type_object == cobra.core.metabolite.Metabolite:
            for key, value in tmp_converted.items():
                converted_bigg_ids[key] = []
                for v in value:
                    if v in self.met_bigg_to_check:
                        converted_bigg_ids[key].append(v)
        if type_object == cobra.core.reaction.Reaction:
            for key, value in tmp_converted.items():
                converted_bigg_ids[key] = []
                for v in value:
                    if v in self.react_bigg_to_check:
                        converted_bigg_ids[key].append(v)
        return converted_bigg_ids

    def prioritiseConvertedSET(self, converted_bigg_ids):
        prior_conv_ids = {
            "1-anno&orig": list(set(converted_bigg_ids["annotation"]) and set(converted_bigg_ids["originalDB"])),
            "2-anno": list(set(converted_bigg_ids["annotation"]) - set(converted_bigg_ids["originalDB"])),
            "3-orig": list(set(converted_bigg_ids["originalDB"]) - set(converted_bigg_ids["annotation"])),
            "4-addit": list(set(converted_bigg_ids["additionalDB"]) - set(converted_bigg_ids["annotation"]) - set(
                converted_bigg_ids["originalDB"])),
            "5-patt": list(set(converted_bigg_ids["pattern"]) - set(converted_bigg_ids["annotation"]) - set(
                converted_bigg_ids["originalDB"]) - set(converted_bigg_ids["additionalDB"])),
            "6-NOconv": list(set(converted_bigg_ids["NOconversion"]) - set(converted_bigg_ids["annotation"]) - set(
                converted_bigg_ids["originalDB"]) - set(converted_bigg_ids["additionalDB"]) - set(
                converted_bigg_ids["pattern"]))}
        return prior_conv_ids
    def prioritiseConverted(self, converted_bigg_ids):
        prior_conv_ids = {
            "1-anno&orig": intersection(converted_bigg_ids["annotation"], converted_bigg_ids["originalDB"]),
            "2-anno": substraction(converted_bigg_ids["annotation"], converted_bigg_ids["originalDB"]),
            "3-orig": substraction(converted_bigg_ids["originalDB"], converted_bigg_ids["annotation"]),
            "4-addit": substraction(converted_bigg_ids["additionalDB"], converted_bigg_ids["annotation"],
                converted_bigg_ids["originalDB"]),
            "5-patt": substraction(converted_bigg_ids["pattern"], converted_bigg_ids["annotation"],
                converted_bigg_ids["originalDB"], converted_bigg_ids["additionalDB"]),
            "6-NOconv": substraction(converted_bigg_ids["NOconversion"], converted_bigg_ids["annotation"],
                converted_bigg_ids["originalDB"], converted_bigg_ids["additionalDB"], converted_bigg_ids["pattern"])}
        return prior_conv_ids

    @abstractmethod
    def runConversion(self, do_compartments: Compartments,
                      original: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction, ) -> [str]:
        """Run conversion process for different ids sources"""


class ConversionForGapseq(ConversionToBiGG):
    """Conversion IDs from original modelSEED ids in gapseq model"""

    def __init__(self, m_original_convert_table: pd.core.frame.DataFrame,
                 r_original_convert_table: pd.core.frame.DataFrame,
                 m_additional_convert_table: pd.core.frame.DataFrame,
                 r_additional_convert_table: pd.core.frame.DataFrame,
                 met_bigg_to_check: list, react_bigg_to_check: list):
        super(ConversionForGapseq, self).__init__(met_bigg_to_check, react_bigg_to_check)
        self.compartment_pattern = ["_c0", "_e0", "_p0"]
        self.m_original_convert_table = m_original_convert_table
        self.r_original_convert_table = r_original_convert_table
        self.m_additional_convert_table = m_additional_convert_table
        self.r_additional_convert_table = r_additional_convert_table
        self.metabolite_annotation = "bigg.metabolite"
        self.reaction_annotation = "bigg.reaction"

    def runConversion(self, do_compartments: Compartments,
                      original: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction, ) -> [str]:
        tmp_converted = {}
        id_wo_compartment = do_compartments.getIDwoCompartment(original.id, self.compartment_pattern)
        if type(original) == cobra.core.metabolite.Metabolite:
            tmp_converted["annotation"] = super().convertViaAnnotation(original.annotation,
                                                                       self.metabolite_annotation)
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.m_original_convert_table)
            tmp_converted["additionalDB"] = super().convertViaTable(id_wo_compartment, self.m_additional_convert_table)
        if type(original) == cobra.core.reaction.Reaction:
            tmp_converted["annotation"] = super().convertViaAnnotation(original.annotation, self.reaction_annotation)
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.r_original_convert_table)
            tmp_converted["additionalDB"] = super().convertViaTable(id_wo_compartment, self.r_additional_convert_table)
        tmp_converted["pattern"] = []
        tmp_converted["NOconversion"] = []
        converted_bigg_ids = super().checkTMPinBigg(type(original), tmp_converted)
        prioritise_bigg_ids = super().prioritiseConverted(converted_bigg_ids)
        return [prioritise_bigg_ids, converted_bigg_ids]


class ConversionForModelseed(ConversionToBiGG):
    """Conversion IDs from original modelSEED ids in modelseed model"""

    def __init__(self, m_original_convert_table: pd.core.frame.DataFrame,
                 r_original_convert_table: pd.core.frame.DataFrame,
                 m_additional_convert_table: pd.core.frame.DataFrame,
                 r_additional_convert_table: pd.core.frame.DataFrame,
                 met_bigg_to_check: list, react_bigg_to_check: list):
        super(ConversionForModelseed, self).__init__(met_bigg_to_check, react_bigg_to_check)
        self.compartment_pattern = ["_c0", "_e0", "_b"]
        self.m_original_convert_table = m_original_convert_table
        self.r_original_convert_table = r_original_convert_table
        self.m_additional_convert_table = m_additional_convert_table
        self.r_additional_convert_table = r_additional_convert_table

    def runConversion(self, do_compartments: Compartments,
                      original: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction) -> [str]:
        tmp_converted = {"annotation": []}
        id_wo_compartment = do_compartments.getIDwoCompartment(original.id, self.compartment_pattern)
        if type(original) == cobra.core.metabolite.Metabolite:
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.m_original_convert_table)
            tmp_converted["additionalDB"] = super().convertViaTable(id_wo_compartment, self.m_additional_convert_table)
        if type(original) == cobra.core.reaction.Reaction:
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.r_original_convert_table)
            tmp_converted["additionalDB"] = super().convertViaTable(id_wo_compartment, self.r_additional_convert_table)
        tmp_converted["pattern"] = []
        tmp_converted["NOconversion"] = []
        converted_bigg_ids = super().checkTMPinBigg(type(original), tmp_converted)
        prioritise_bigg_ids = super().prioritiseConverted(converted_bigg_ids)
        return [prioritise_bigg_ids, converted_bigg_ids]


class ConversionForAgora(ConversionToBiGG):
    """Conversion IDs from original wierd BiGG (RECONE) ids in AGORA model"""

    def __init__(self, m_original_convert_table: pd.core.frame.DataFrame,
                 r_original_convert_table: pd.core.frame.DataFrame,
                 m_additional_convert_table: pd.core.frame.DataFrame,
                 r_additional_convert_table: pd.core.frame.DataFrame,
                 met_bigg_to_check: list, react_bigg_to_check: list):
        super(ConversionForAgora, self).__init__(met_bigg_to_check, react_bigg_to_check)
        self.compartment_pattern = ["(c)", "(e)"]
        self.m_original_convert_table = m_original_convert_table
        self.r_original_convert_table = r_original_convert_table
        self.m_additional_convert_table = m_additional_convert_table
        self.r_additional_convert_table = r_additional_convert_table
        self.metabolite_annotation = "kegg.compound"
        self.reaction_annotation = "kegg.reaction"

    def convertViaPattern(self, id_wo_compartment: str) -> [str]:
        converted = [id_wo_compartment[::-1].replace("_", "__", 1)[::-1]]
        return converted

    def runConversion(self, do_compartments: Compartments,
                      original: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction) -> [str]:
        tmp_converted = {"annotation": []}
        id_wo_compartment = do_compartments.getIDwoCompartment(original.id, self.compartment_pattern)
        if type(original) == cobra.core.metabolite.Metabolite:
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.m_original_convert_table)
            met_from_annotation = super().convertViaAnnotation(original.annotation, self.metabolite_annotation)
            tmp_converted["additionalDB"] = []
            for met in met_from_annotation:
                tmp_converted["additionalDB"] = tmp_converted["additionalDB"] + super().convertViaTable(met,
                                                                                                        self.m_additional_convert_table)
        if type(original) == cobra.core.reaction.Reaction:
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.r_original_convert_table)
            react_from_annotation = super().convertViaAnnotation(original.annotation, self.reaction_annotation)
            tmp_converted["additionalDB"] = []
            for react in react_from_annotation:
                tmp_converted["additionalDB"] = tmp_converted["additionalDB"] + super().convertViaTable(react,
                                                                                                        self.r_additional_convert_table)
        tmp_converted["pattern"] = self.convertViaPattern(id_wo_compartment)
        tmp_converted["NOconversion"] = [id_wo_compartment]
        converted_bigg_ids = super().checkTMPinBigg(type(original), tmp_converted)
        prioritise_bigg_ids = super().prioritiseConverted(converted_bigg_ids)
        return [prioritise_bigg_ids, converted_bigg_ids]


# endregion

def runConversionForALLmodels(model_types: [str], all_models: dict, CompartStrategies: dict, ConvertStrategies: dict,
                                obj_type: "metabolites" or "reactions", write_output=True, useroutname=None) -> dict:
    if write_output:
        if useroutname is not  None:
            output_ids = open("../Output/" + useroutname + "_" + obj_type + "_ids_conversion_output.tsv", "w")
            output_number = open("../Output/" + useroutname + "_" + obj_type + "_numbers_conversion_output.tsv", "w")
            output_NOprior = open("../Output/" + useroutname + "_" + obj_type + "_Noprior_conversion_output.tsv", "w")
        else:
            output_ids = open("../Output/" + obj_type + "_ids_conversion_output.tsv", "w")
            output_number = open("../Output/" + obj_type + "_numbers_conversion_output.tsv", "w")
            output_NOprior = open("../Output/" + obj_type + "_Noprior_conversion_output.tsv", "w")
        output_ids.write("Model_type\tID\t1-anno&orig\t2-anno\t3-orig\t4-addit\t5-patt\t6-NOconv\n")
        output_number.write("Model_type\tID\t1-anno&orig\t2-anno\t3-orig\t4-addit\t5-patt\t6-NOconv\n")
        output_NOprior.write("Model_type\tID\tanno\torig\taddit\tpatt\tNOconv\n")
    all_converted = {}
    for typ in model_types:
        all_converted.update({typ: {}})
        if obj_type == "metabolites":
            objects = all_models.get(typ).metabolites
        if obj_type == "reactions":
            objects = all_models.get(typ).reactions
        for obj in objects:
            ids = ConvertStrategies.get(typ).runConversion(CompartStrategies.get(typ), obj)
            bigg_ids = ids[0]
            noprior_ids = ids[1]
            compartments = CompartStrategies.get(typ).runGetCompartment(obj)
            if write_output:
                ids_string = f"{typ}\t{obj.id}\t{bigg_ids['1-anno&orig']}\t{bigg_ids['2-anno']}\t{bigg_ids['3-orig']}" \
                             f"\t{bigg_ids['4-addit']}\t{bigg_ids['5-patt']}\t{bigg_ids['6-NOconv']}\n"
                numbers_string = f"{typ}\t{obj.id}\t{len(bigg_ids['1-anno&orig'])}\t{len(bigg_ids['2-anno'])}\t" \
                                 f"{len(bigg_ids['3-orig'])}\t{len(bigg_ids['4-addit'])}\t{len(bigg_ids['5-patt'])}" \
                                 f"\t{len(bigg_ids['6-NOconv'])}\n"
                noprior_ids_string = f"{typ}\t{obj.id}\t{noprior_ids['annotation']}\t{noprior_ids['originalDB']}" \
                             f"\t{noprior_ids['additionalDB']}\t{noprior_ids['pattern']}\t{noprior_ids['NOconversion']}\n"
                output_ids.write(ids_string)
                output_number.write(numbers_string)
                output_NOprior.write(noprior_ids_string)
            all_converted.get(typ).update({obj.id: [compartments, bigg_ids]})
    if write_output:
        output_ids.close()
        output_number.close()
        output_NOprior.close()
    return all_converted

if __name__ == '__main__':
    pass
