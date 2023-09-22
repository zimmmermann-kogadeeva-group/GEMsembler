import re
from abc import ABC, abstractmethod
import cobra
import pandas as pd


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

    def prioritiseConverted(self, converted_bigg_ids):
        prior_conv_ids = {
            "1-anno&orig": list(
                set(converted_bigg_ids["annotation"]).intersection(set(converted_bigg_ids["originalDB"]))),
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

    def addCompartmentToConvertedM(self, compartment: str, converted_m: dict):
        converted_m_c = {}
        for level, ids in converted_m.items():
            converted_m_c[level] = [i + "_" + compartment for i in ids]
        return converted_m_c

    @abstractmethod
    def runConversion(self, original: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction) -> [str]:
        """Run conversion process for different ids sources"""

    def getIDwoCompartment(self, original_id: str, pattern_to_remove: str) -> str:
        id_wo_compartment = re.sub(pattern_to_remove, "", original_id)
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

class ConversionForGapseq(ConversionToBiGG):
    """Conversion IDs from original modelSEED ids in gapseq model"""

    def __init__(self, m_original_convert_table: pd.core.frame.DataFrame,
                 r_original_convert_table: pd.core.frame.DataFrame,
                 m_additional_convert_table: pd.core.frame.DataFrame,
                 r_additional_convert_table: pd.core.frame.DataFrame,
                 met_bigg_to_check: list, react_bigg_to_check: list):
        super(ConversionForGapseq, self).__init__(met_bigg_to_check, react_bigg_to_check)
        self.compartment_pattern = "_(c0|e0|p0)$"
        self.m_original_convert_table = m_original_convert_table
        self.r_original_convert_table = r_original_convert_table
        self.m_additional_convert_table = m_additional_convert_table
        self.r_additional_convert_table = r_additional_convert_table
        self.metabolite_annotation = "bigg.metabolite"
        self.reaction_annotation = "bigg.reaction"

    def getCompartment(self, original_compartment: str) -> str:
        standard_compartment = original_compartment[:-1]
        return standard_compartment

    def runConversion(self, original: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction) -> [str]:
        tmp_converted = {}
        id_wo_compartment = super().getIDwoCompartment(original.id, self.compartment_pattern)
        if type(original) == cobra.core.metabolite.Metabolite:
            tmp_converted["annotation"] = super().convertViaAnnotation(original.annotation,
                                                                       self.metabolite_annotation)
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.m_original_convert_table)
            tmp_converted["additionalDB"] = super().convertViaTable(id_wo_compartment, self.m_additional_convert_table)
        elif type(original) == cobra.core.reaction.Reaction:
            tmp_converted["annotation"] = super().convertViaAnnotation(original.annotation, self.reaction_annotation)
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.r_original_convert_table)
            tmp_converted["additionalDB"] = super().convertViaTable(id_wo_compartment, self.r_additional_convert_table)
        tmp_converted["pattern"] = []
        tmp_converted["NOconversion"] = []
        converted_bigg_ids = super().checkTMPinBigg(type(original), tmp_converted)
        prioritise_bigg_ids = super().prioritiseConverted(converted_bigg_ids)
        compartments = super().runGetCompartment(original)
        if type(original) == cobra.core.metabolite.Metabolite:
            prioritise_bigg_ids = super().addCompartmentToConvertedM(compartments[0], prioritise_bigg_ids)
        return [prioritise_bigg_ids, converted_bigg_ids, id_wo_compartment, compartments]


class ConversionForModelseed(ConversionToBiGG):
    """Conversion IDs from original modelSEED ids in modelseed model"""

    def __init__(self, m_original_convert_table: pd.core.frame.DataFrame,
                 r_original_convert_table: pd.core.frame.DataFrame,
                 m_additional_convert_table: pd.core.frame.DataFrame,
                 r_additional_convert_table: pd.core.frame.DataFrame,
                 met_bigg_to_check: list, react_bigg_to_check: list):
        super(ConversionForModelseed, self).__init__(met_bigg_to_check, react_bigg_to_check)
        self.compartment_pattern = "_(c0|e0|b)$"
        self.m_original_convert_table = m_original_convert_table
        self.r_original_convert_table = r_original_convert_table
        self.m_additional_convert_table = m_additional_convert_table
        self.r_additional_convert_table = r_additional_convert_table

    def getCompartment(self, original_compartment: str) -> str:
        standard_compartment = original_compartment[:-1]
        return standard_compartment

    def runConversion(self, original: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction) -> [str]:
        tmp_converted = {"annotation": []}
        id_wo_compartment = super().getIDwoCompartment(original.id, self.compartment_pattern)
        if type(original) == cobra.core.metabolite.Metabolite:
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.m_original_convert_table)
            tmp_converted["additionalDB"] = super().convertViaTable(id_wo_compartment, self.m_additional_convert_table)
        elif type(original) == cobra.core.reaction.Reaction:
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.r_original_convert_table)
            tmp_converted["additionalDB"] = super().convertViaTable(id_wo_compartment, self.r_additional_convert_table)
        tmp_converted["pattern"] = []
        tmp_converted["NOconversion"] = []
        converted_bigg_ids = super().checkTMPinBigg(type(original), tmp_converted)
        prioritise_bigg_ids = super().prioritiseConverted(converted_bigg_ids)
        compartments = super().runGetCompartment(original)
        if type(original) == cobra.core.metabolite.Metabolite:
            prioritise_bigg_ids = super().addCompartmentToConvertedM(compartments[0], prioritise_bigg_ids)
        return [prioritise_bigg_ids, converted_bigg_ids, id_wo_compartment, compartments]


class ConversionForAgora(ConversionToBiGG):
    """Conversion IDs from original wierd BiGG (RECONE) ids in AGORA model"""

    def __init__(self, m_original_convert_table: pd.core.frame.DataFrame,
                 r_original_convert_table: pd.core.frame.DataFrame,
                 m_additional_convert_table: pd.core.frame.DataFrame,
                 r_additional_convert_table: pd.core.frame.DataFrame,
                 met_bigg_to_check: list, react_bigg_to_check: list):
        super(ConversionForAgora, self).__init__(met_bigg_to_check, react_bigg_to_check)
        self.compartment_pattern = "\[(c|e)\]$"
        self.m_original_convert_table = m_original_convert_table
        self.r_original_convert_table = r_original_convert_table
        self.m_additional_convert_table = m_additional_convert_table
        self.r_additional_convert_table = r_additional_convert_table
        self.metabolite_annotation = "kegg.compound"
        self.reaction_annotation = "kegg.reaction"

    def getCompartment(self, original_compartment: str) -> str:
        standard_compartment = original_compartment
        return standard_compartment

    def convertViaPattern(self, id_wo_compartment: str) -> [str]:
        converted = [id_wo_compartment[::-1].replace("_", "__", 1)[::-1]]
        return converted

    def runConversion(self, original: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction) -> [str]:
        tmp_converted = {"annotation": []}
        id_wo_compartment = super().getIDwoCompartment(original.id, self.compartment_pattern)
        if type(original) == cobra.core.metabolite.Metabolite:
            tmp_converted["originalDB"] = super().convertViaTable(id_wo_compartment, self.m_original_convert_table)
            met_from_annotation = super().convertViaAnnotation(original.annotation, self.metabolite_annotation)
            tmp_converted["additionalDB"] = []
            for met in met_from_annotation:
                tmp_converted["additionalDB"] = tmp_converted["additionalDB"] + super().convertViaTable(met,
                                                                                                        self.m_additional_convert_table)
        elif type(original) == cobra.core.reaction.Reaction:
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
        compartments = super().runGetCompartment(original)
        if type(original) == cobra.core.metabolite.Metabolite:
            prioritise_bigg_ids = super().addCompartmentToConvertedM(compartments[0], prioritise_bigg_ids)
        return [prioritise_bigg_ids, converted_bigg_ids, id_wo_compartment, compartments]

# endregion

class ConversionForCarveMe(ConversionToBiGG):
    """Conversion IDs from original modelSEED ids in gapseq model"""

    def __init__(self, m_original_convert_table: pd.core.frame.DataFrame,
                 r_original_convert_table: pd.core.frame.DataFrame,
                 met_bigg_to_check: list, react_bigg_to_check: list):
        super(ConversionForCarveMe, self).__init__(met_bigg_to_check, react_bigg_to_check)
        self.m_original_convert_table = m_original_convert_table
        self.r_original_convert_table = r_original_convert_table
        self.compartment_pattern = "_(c|e|p)$"

    def getCompartment(self, original_compartment: str) -> str:
        standard_compartment = original_compartment.split("_")[1]
        return standard_compartment

    def runConversion(self, original: cobra.core.metabolite.Metabolite or cobra.core.reaction.Reaction) -> [str]:
        compartments = super().runGetCompartment(original)
        if type(original) == cobra.core.metabolite.Metabolite:
            id_wo_compartment = super().getIDwoCompartment(original.id, self.compartment_pattern)
            if id_wo_compartment not in self.met_bigg_to_check:
                bigg_ids = super().convertViaTable(id_wo_compartment, self.m_original_convert_table)
                if bigg_ids:
                    bigg_ids = [b+"_"+compartments[0] for b in bigg_ids]
                else:
                    bigg_ids = ["not_found_in_new_and_old_bigg"]
            else:
                bigg_ids = [original.id]
        elif type(original) == cobra.core.reaction.Reaction:
            id_wo_compartment = original.id
            if id_wo_compartment not in self.react_bigg_to_check:
                bigg_ids = super().convertViaTable(id_wo_compartment, self.r_original_convert_table)
                if not bigg_ids:
                    bigg_ids = ["not_found_in_new_and_old_bigg"]
            else:
                bigg_ids = [original.id]
        return [compartments, bigg_ids]



def summarizeConversion(model_types: [str], obj_type: "metabolites" or "reactions", useroutname=None):
    """ Summary for conversion. Number of converted for different levels. """
    levels = ['1-anno&orig', '2-anno', '3-orig', '4-addit', '5-patt', '6-NOconv']
    if useroutname is not None:
        outdata = pd.read_csv("../Output/" + useroutname + "_" + obj_type + "_numbers_conversion_output.tsv", sep="\t")
        summarydata = open("../Output/" + useroutname + "_" + obj_type + "_conversion_summary.txt", "w")
    else:
        outdata = pd.read_csv("../Output/" + obj_type + "_numbers_conversion_output.tsv", sep="\t")
        summarydata = open("../Output/" + obj_type + "_conversion_summary.txt", "w")
    for typ in model_types:
        number_uniq = []
        number_not_uniq = []
        number_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 1)].index))
        number_not_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] >= 1)].index))
        number_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                                       (outdata["2-anno"] == 1)].index))
        number_not_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                                      (outdata["2-anno"] >= 1)].index))
        number_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                                       (outdata["2-anno"] == 0) & (outdata["3-orig"] == 1)].index))
        number_not_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                                       (outdata["2-anno"] == 0) & (outdata["3-orig"] >= 1)].index))
        number_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                                       (outdata["2-anno"] == 0) & (outdata["3-orig"] == 0) &
                                       (outdata["4-addit"] == 1)].index))
        number_not_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                                       (outdata["2-anno"] == 0) & (outdata["3-orig"] == 0) &
                                       (outdata["4-addit"] >= 1)].index))
        number_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                                       (outdata["2-anno"] == 0) & (outdata["3-orig"] == 0) &
                                       (outdata["4-addit"] == 0) & (outdata["5-patt"] == 1)].index))
        number_not_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                                       (outdata["2-anno"] == 0) & (outdata["3-orig"] == 0) &
                                       (outdata["4-addit"] == 0) & (outdata["5-patt"] >= 1)].index))
        number_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                                       (outdata["2-anno"] == 0) & (outdata["3-orig"] == 0) &
                                       (outdata["4-addit"] == 0) & (outdata["5-patt"] == 0) & (
                                               outdata["6-NOconv"] == 1)].index))
        number_not_uniq.append(len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                                       (outdata["2-anno"] == 0) & (outdata["3-orig"] == 0) &
                                       (outdata["4-addit"] == 0) & (outdata["5-patt"] == 0) & (
                                               outdata["6-NOconv"] >= 1)].index))
        high = len(outdata[(outdata["Model_type"] == typ) & (outdata[["1-anno&orig", "2-anno", "3-orig"]] >= 1).any(
            axis=1)].index)
        middle = len(outdata[(outdata["Model_type"] == typ) & (outdata["1-anno&orig"] == 0) &
                             (outdata["2-anno"] == 0) & (outdata["3-orig"] == 0) &
                             (outdata[["4-addit", "5-patt", "6-NOconv"]] >= 1).any(axis=1)].index)
        low_no = len(outdata[(outdata["Model_type"] == typ) & (outdata[["1-anno&orig", "2-anno", "3-orig", "4-addit",
                                                                        "5-patt", "6-NOconv"]] == 0).all(axis=1)].index)
        uniq_id = len(set(outdata[outdata['Model_type'] == typ]['ID'].tolist()))
        print(f"{uniq_id} {obj_type} in {typ} model")
        summarydata.write(f"{uniq_id} {obj_type} in {typ} model\n")
        for i in range(len(levels)):
            print(f"For {typ} models {number_uniq[i]}  {obj_type} were converted uniquely with {levels[i]} level")
            summarydata.write(
                f"For {typ} models {number_uniq[i]}  {obj_type} were converted uniquely with {levels[i]} level\n")
            print(f"For {typ} models {number_not_uniq[i]}  {obj_type} were converted with {levels[i]} level")
            summarydata.write(
                f"For {typ} models {number_not_uniq[i]}  {obj_type} were converted with {levels[i]} level\n")
        print(f"For {typ} models {high}  {obj_type} were converted with 1st, 2d or 3d level")
        print(f"For {typ} models {middle}  {obj_type} were converted with 4th, 5th or 6th level")
        print(f"For {typ} models {low_no}  {obj_type} were not converted at all")
        summarydata.write(f"For {typ} models {high}  {obj_type} were converted with 1st, 2d or 3d level\n")
        summarydata.write(f"For {typ} models {middle}  {obj_type} were converted with 4th, 5th or 6th level\n")
        summarydata.write(f"For {typ} models {low_no}  {obj_type} were not converted at all\n")
    summarydata.close()


# def findManyToOne(model_types: [str], obj_type: "metabolites" or "reactions", useroutname: str):
#     levels = ['1-anno&orig', '2-anno', '3-orig', '4-addit', '5-patt', '6-NOconv']
#     if useroutname is not None:
#         outdata = pd.read_csv("../Output/" + useroutname + "_" + obj_type + "_ids_conversion_output.tsv", sep="\t")
#         file_name = "../Output/" + useroutname + "_" + obj_type + "_many_to_one_conversion_.tsv"
#     else:
#         outdata = pd.read_csv("../Output/" + obj_type + "_ids_conversion_output.tsv", sep="\t")
#         file_name = "../Output/" + obj_type + "_many_to_one_conversion_.tsv"
#     outdata = outdata.drop('ID_c', axis=1)
#     outdata = outdata.drop_duplicates(keep="first")
#     for l in levels:
#         outdata[l] = outdata[l].apply(literal_eval)
#     repeated_data = pd.DataFrame(
#         columns=['Model_type', 'ID', '1-anno&orig', '2-anno', '3-orig', '4-addit', '5-patt', '6-NOconv'])
#     for typ in model_types:
#         for lev in levels:
#             ids_list_in_list = outdata[(outdata["Model_type"] == typ)][lev].tolist()
#             ids_list = [item for sublist in ids_list_in_list for item in sublist]
#             ids_occurence = Counter(ids_list)
#             ids_repeated = general.findKeysByValue(ids_occurence, 1, operator.gt)
#             if ids_repeated != []:
#                 for rep in ids_repeated:
#                     repeated_record = outdata[
#                         (outdata["Model_type"] == typ) & (outdata[lev].str.contains(rep, regex=False))]
#                     orig_id = repeated_record["ID"].values[0]
#                     if orig_id not in repeated_data[repeated_data["Model_type"] == typ]["ID"].tolist():
#                         repeated_data = pd.concat([repeated_data, repeated_record], ignore_index=True)
#     repeated_data.to_csv(file_name, sep='\t')
#
# def addCompartmentToConvertedM(compartment: str, converted_m: dict):
#     converted_m_c = {}
#     for level, ids in converted_m.items():
#         converted_m_c[level] = [i + "_" + compartment for i in ids]
#     return converted_m_c

def runConversionForALLmodels(model_types: [str], all_models: dict, ConvertStrategies: dict,
                              obj_type: "metabolites" or "reactions",
                              write_output=True, do_summary=True, do_many_to_one=True,
                              useroutname=None) -> dict:
    """ Running conversion for all models with conversion needed. """
    all_converted = {}
    ids_strings = []
    numbers_strings = []
    noprior_ids_strings = []
    for typ in model_types:
        all_converted.update({typ: {}})

        if obj_type == "metabolites":
            objects = all_models.get(typ).metabolites
        elif obj_type == "reactions":
            objects = all_models.get(typ).reactions

        for obj in objects:
            ids = ConvertStrategies.get(typ).runConversion(obj)
            bigg_ids = ids[0]
            noprior_ids = ids[1]
            id_wo_compartment = ids[2]
            compartments = ids[3]
            # if obj_type == "metabolites":
            #     bigg_ids = addCompartmentToConvertedM(compartments[0], bigg_ids)
            all_converted.get(typ).update({obj.id: [compartments, bigg_ids]})
            ids_strings.append(
                f"{typ}\t{obj.id}\t{bigg_ids['1-anno&orig']}\t{bigg_ids['2-anno']}\t{bigg_ids['3-orig']}"
                f"\t{bigg_ids['4-addit']}\t{bigg_ids['5-patt']}\t{bigg_ids['6-NOconv']}\n")
            numbers_strings.append(
                f"{typ}\t{obj.id}\t{len(bigg_ids['1-anno&orig'])}\t{len(bigg_ids['2-anno'])}\t"
                f"{len(bigg_ids['3-orig'])}\t{len(bigg_ids['4-addit'])}\t{len(bigg_ids['5-patt'])}"
                f"\t{len(bigg_ids['6-NOconv'])}\n")
            noprior_ids_strings.append(
                f"{typ}\t{obj.id}\t{noprior_ids['annotation']}\t{noprior_ids['originalDB']}"
                f"\t{noprior_ids['additionalDB']}\t{noprior_ids['pattern']}\t{noprior_ids['NOconversion']}\n")
    if write_output:
        if useroutname is not None:
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
        output_ids.writelines(ids_strings)
        output_number.writelines(numbers_strings)
        output_NOprior.writelines(noprior_ids_strings)
        output_ids.close()
        output_number.close()
        output_NOprior.close()
        if do_summary:
            summarizeConversion(model_types, obj_type, useroutname)
        # if do_many_to_one:
        #     findManyToOne(model_types, obj_type, useroutname)
    return all_converted

def runNoneConversionChecking(model_types: [str], all_models: dict, ConvertStrategies: dict,
                              obj_type: "metabolites" or "reactions") -> dict:
    """ Checking IDs for all models with no need for conversion that supposed to have BiGG IDs originally,
    but maybe old or incorrect. """
    all_checked = {}
    all_not_in_bigg = {}
    for typ in model_types:
        all_checked.update({typ: {}})
        all_not_in_bigg.update({typ: {}})
        if obj_type == "metabolites":
            objects = all_models.get(typ).metabolites
        elif obj_type == "reactions":
            objects = all_models.get(typ).reactions

        for obj in objects:
            ids = ConvertStrategies.get(typ).runConversion(obj)
            if ids[1][0] != "not_found_in_new_and_old_bigg":
                all_checked.get(typ).update({obj.id: ids})
            else:
                all_not_in_bigg.get(typ).update({obj.id: ids})
    return all_checked, all_not_in_bigg
