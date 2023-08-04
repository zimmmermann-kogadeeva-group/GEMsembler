import re

from cobra import Model
import pandas as pd
import os
from sympy import symbols, sympify, expand
from general import is_float

def getGenome(ncbi_genome_name: str):
    genomedata = open(ncbi_genome_name, "r")
    genome_lines = genomedata.readlines()
    genomedata.close()
    genomes = {}
    genomeid = ""
    for line in genome_lines:
        if line.startswith(">"):
            if genomeid != "":
                genomes[genomeid] = genome
            genomeid = line.split()[0][1:]
            genome = ""
        else:
            genome = genome + line.strip()
    genomes[genomeid] = genome
    return genomes


def getGenesGapseq(gapseq_genes_name: str, gapseq_model: Model, genomes: dict):
    gene_gapseq_fasta = open(gapseq_genes_name, "w")
    for gene in gapseq_model.genes:
        gene_gapseq_fasta.write(">" + gene.id + "\n")
        start = gene.id.split("_")[-2]
        end = gene.id.split("_")[-1]
        genomeid = gene.id.removeprefix("gp_").removesuffix("_"+start+"_"+end)
        genomeid = ".".join(genomeid.rsplit("_", 1))
        start = int(start)
        end = int(end)
        if start < end:
            gene_gapseq_fasta.write(genomes.get(genomeid)[start:end] + "\n")
        else:
            gene_gapseq_fasta.write(genomes.get(genomeid)[end:start] + "\n")
    gene_gapseq_fasta.close()


def getLocusTagGenes(ncbi_cds_name:str, ncbi_protein_name: str, feature_table_name:str,
                     out_nt_fasta: str, out_aa_fasta: str, do_old_locus_tag=False):
    feature_table = pd.read_csv(feature_table_name, sep="\t")
    with open(ncbi_cds_name, "r") as cds:
        cds_lines = cds.readlines()
    out_nt = open(out_nt_fasta, "w")
    for cdsl in cds_lines:
        if cdsl.startswith(">"):
            new_locus_tag = cdsl.split("[locus_tag=")[1].split("]")[0]
            if do_old_locus_tag:
                attr = feature_table[(feature_table["locus_tag"]==new_locus_tag) & (feature_table["# feature"]=="gene")]["attributes"]
                if attr.empty:
                    old_locus_tag = new_locus_tag
                elif type(attr.values[0]) != str:
                    old_locus_tag = new_locus_tag
                elif "old_locus_tag" not in attr.values[0]:
                    old_locus_tag = new_locus_tag
                else:
                    old_locus_tag = attr.values[0].split("old_locus_tag=")[1]
                out_nt.write(">" + old_locus_tag + "\n")
            else:
                out_nt.write(">"+new_locus_tag+"\n")
        else:
            out_nt.write(cdsl)
    out_nt.close()
    with open(ncbi_protein_name, "r") as proteins:
        proteins_lines = proteins.readlines()
    out_aa = open(out_aa_fasta, "w")
    for protl in proteins_lines:
        if protl.startswith(">"):
            pp_id = protl.split(" ")[0][1:]
            pnew_locus_tag = \
                feature_table[(feature_table["product_accession"] == pp_id) & (feature_table["# feature"] == "CDS")][
                    "locus_tag"].values[0]
            if do_old_locus_tag:
                attr = \
                feature_table[(feature_table["locus_tag"] == pnew_locus_tag) & (feature_table["# feature"] == "gene")][
                    "attributes"]
                if attr.empty:
                    pold_locus_tag = pnew_locus_tag
                elif type(attr.values[0]) != str:
                    pold_locus_tag = pnew_locus_tag
                elif "old_locus_tag" not in attr.values[0]:
                    pold_locus_tag = pnew_locus_tag
                else:
                    pold_locus_tag = attr.values[0].split("old_locus_tag=")[1]
                out_aa.write(">" + pold_locus_tag + "\n")
            else:
                out_aa.write(">" + pnew_locus_tag + "\n")
        else:
            out_aa.write(protl)
    out_aa.close()



def runGenesConversion(model_type: [str], all_models: dict, input_genomes_names: dict, gapseq_genes_names:dict, ncbi_cds_name:str, ncbi_protein_name: str, feature_table_name:str,
                     out_nt_fasta: str, out_aa_fasta: str, do_old_locus_tag=False):
    os.chdir("../Data/")
    getLocusTagGenes(ncbi_cds_name, ncbi_protein_name, feature_table_name, out_nt_fasta, out_aa_fasta, do_old_locus_tag)
    os.system(f"makeblastdb -in {out_nt_fasta} -out ncbi_cds -dbtype nucl -title 'ncbi_cds' -parse_seqids")
    os.system(f"makeblastdb -in {out_aa_fasta} -out ncbi_proteins -dbtype prot -title 'ncbi_proteins' -parse_seqids")
    for typ in model_type:
        if typ == "gapseq":
            genomes = getGenome(input_genomes_names.get(typ))
            getGenesGapseq(gapseq_genes_names.get(typ), all_models.get(typ), genomes)
            os.system(
                f"blastn -query {gapseq_genes_names.get(typ)} -db ncbi_cds -max_target_seqs 1 -outfmt '6' -out {typ}_blast.tsv")
        elif typ == "agora":
            os.system(
                f"blastn -query {input_genomes_names.get(typ)} -db ncbi_cds -max_target_seqs 1 -outfmt '6' -out {typ}_blast.tsv")
        else:
            os.system(
                f"blastp -query {input_genomes_names.get(typ)} -db ncbi_proteins -max_target_seqs 1 -outfmt '6' -out {typ}_blast.tsv")


def makeNewGPR(gpr: str, g_id_convert: dict):
    new_gpr = gpr
    for old_id, new_id in g_id_convert.items():
        new_gpr = re.sub(rf'\b{old_id}\b', new_id, new_gpr)
        if new_id != "not_found":
            globals()[new_id] = symbols(new_id)
    new_gpr= re.sub(r'\bnot_found\b', "1", new_gpr)
    new_gpr = re.sub(r'\bor\b', "+", new_gpr)
    new_gpr = re.sub(r'\band\b', "*", new_gpr)
    equation = sympify(new_gpr)
    equation_expanded = expand(equation)
    if is_float(str(equation_expanded)):
        new_gpr_expanded = ""
    else:
        new_gpr_expanded = str(equation_expanded).split(" + ")
        for i in range(len(new_gpr_expanded)):
            genes_and =[]
            for g in new_gpr_expanded[i].split("*"):
                if g and (not is_float(g)):
                    genes_and.append(g)
            if not genes_and:
                new_gpr_expanded.pop(i)
            elif (len(genes_and) > 1) and (len(new_gpr_expanded) > 1):
                new_gpr_expanded[i] = "(" + " and ".join(sorted(genes_and)) + ")"
            else:
                new_gpr_expanded[i] = " and ".join(sorted(genes_and))
        new_gpr_expanded = " or ".join(sorted(new_gpr_expanded))
    return new_gpr_expanded


def uniteGPR(list_gpr: [str]):
    united_and = []
    for gpr in list_gpr:
        for genes_ands in gpr.split(" or "):
            united_and.append(genes_ands)
    unique_and = list(set(united_and))
    if len(unique_and) > 1:
        for i in range(len(unique_and)):
            if ("and" in unique_and[i]) and ("(" not in unique_and[i]):
                unique_and[i] = "(" + unique_and[i] + ")"
    united_gpr = " or ".join(unique_and)
    return united_gpr
