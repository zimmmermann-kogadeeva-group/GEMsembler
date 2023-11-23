import re
from cobra import Model
import pandas as pd
import os
from sympy import symbols, sympify, expand
from general import is_float


def checkNTorAA(path_fasta: str):
    """Check whether fasta file is nt or aa.
    Codes are taken from  https://web.cas.org/help/BLAST/topics/codes.htm"""
    nt_letters = [
        "A",
        "C",
        "G",
        "T",
        "U",
        "R",
        "Y",
        "K",
        "M",
        "S",
        "W",
        "B",
        "D",
        "H",
        "V",
        "N",
        "-",
    ]
    aa_letters = [
        "A",
        "B",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "J",
        "K",
        "L",
        "M",
        "N",
        "O",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "U",
        "V",
        "W",
        "X",
        "Y",
        "Z",
        "*",
        "-",
    ]
    aa_specific = list(set(aa_letters) - set(nt_letters))
    with open(path_fasta) as fasta_file:
        file_lines = fasta_file.readlines()
    sequence = ""
    for line in file_lines:
        if not line.startswith(">"):
            sequence = sequence + line.strip()
    aa_status = False
    for aa in aa_specific:
        if aa in sequence:
            aa_status = True
    return aa_status


def checkGeneIDs(path_fasta: str, model: Model):
    gene_ids = [g.id for g in model.genes]
    with open(path_fasta) as fasta_file:
        file_lines = fasta_file.readlines()
    ids_intersect = []
    for line in file_lines:
        if line.startswith(">"):
            id_seq = line.split()[0][1:]
            if id_seq in gene_ids:
                ids_intersect.append(id_seq)
    if ids_intersect:
        ids_same = True
    else:
        ids_same = False
    return ids_same


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


def getGenesGapseq(output_gapseq_genes_name: str, gapseq_model: Model, genomes: dict):
    gene_gapseq_fasta = open(output_gapseq_genes_name, "w")
    for gene in gapseq_model.genes:
        gene_gapseq_fasta.write(">" + gene.id + "\n")
        start = gene.id.split("_")[-2]
        end = gene.id.split("_")[-1]
        genomeid = gene.id.removeprefix("gp_").removesuffix("_" + start + "_" + end)
        genomeid = ".".join(genomeid.rsplit("_", 1))
        start = int(start)
        end = int(end)
        if start < end:
            gene_gapseq_fasta.write(genomes.get(genomeid)[start:end] + "\n")
        else:
            gene_gapseq_fasta.write(genomes.get(genomeid)[end:start] + "\n")
    gene_gapseq_fasta.close()


def getLocusTagGenes(
    ncbi_cds_name: str,
    ncbi_protein_name: str,
    feature_table_name: str,
    out_nt_fasta: str,
    out_aa_fasta: str,
    do_old_locus_tag=False,
):
    feature_table = pd.read_csv(feature_table_name, sep="\t")
    with open(ncbi_cds_name, "r") as cds:
        cds_lines = cds.readlines()
    out_nt = open(out_nt_fasta, "w")
    for cdsl in cds_lines:
        if cdsl.startswith(">"):
            new_locus_tag = cdsl.split("[locus_tag=")[1].split("]")[0]
            if do_old_locus_tag:
                attr = feature_table[
                    (feature_table["locus_tag"] == new_locus_tag)
                    & (feature_table["# feature"] == "gene")
                ]["attributes"]
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
                out_nt.write(">" + new_locus_tag + "\n")
        else:
            out_nt.write(cdsl)
    out_nt.close()
    with open(ncbi_protein_name, "r") as proteins:
        proteins_lines = proteins.readlines()
    out_aa = open(out_aa_fasta, "w")
    for protl in proteins_lines:
        if protl.startswith(">"):
            pp_id = protl.split(" ")[0][1:]
            pnew_locus_tag = feature_table[
                (feature_table["product_accession"] == pp_id)
                & (feature_table["# feature"] == "CDS")
            ]["locus_tag"].values[0]
            if do_old_locus_tag:
                attr = feature_table[
                    (feature_table["locus_tag"] == pnew_locus_tag)
                    & (feature_table["# feature"] == "gene")
                ]["attributes"]
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


def runGenesConversion(
    model_type: [str],
    all_models: dict,
    input_genomes_names: dict,
    gapseq_genes_names: dict,
    ncbi_cds_name: str,
    ncbi_protein_name: str,
    feature_table_name: str,
    out_nt_fasta: str,
    out_aa_fasta: str,
    do_old_locus_tag=True,
):
    os.chdir("../Data/")
    getLocusTagGenes(
        ncbi_cds_name,
        ncbi_protein_name,
        feature_table_name,
        out_nt_fasta,
        out_aa_fasta,
        do_old_locus_tag,
    )
    os.system(
        f"makeblastdb -in {out_nt_fasta} -out ncbi_cds -dbtype nucl -title 'ncbi_cds' -parse_seqids"
    )
    os.system(
        f"makeblastdb -in {out_aa_fasta} -out ncbi_proteins -dbtype prot -title 'ncbi_proteins' -parse_seqids"
    )
    for typ in model_type:
        if typ == "gapseq":
            genomes = getGenome(input_genomes_names.get(typ))
            getGenesGapseq(gapseq_genes_names.get(typ), all_models.get(typ), genomes)
            os.system(
                f"blastn -query {gapseq_genes_names.get(typ)} -db ncbi_cds -max_target_seqs 1 -outfmt '6' -out {typ}_blast.tsv"
            )
        elif typ == "agora":
            os.system(
                f"blastn -query {input_genomes_names.get(typ)} -db ncbi_cds -max_target_seqs 1 -outfmt '6' -out {typ}_blast.tsv"
            )
        else:
            os.system(
                f"blastp -query {input_genomes_names.get(typ)} -db ncbi_proteins -max_target_seqs 1 -outfmt '6' -out {typ}_blast.tsv"
            )


def makeNewGPR(gpr: str, g_id_convert: dict):
    new_gpr = gpr
    mix_gpr = gpr
    for old_id, new_id in g_id_convert.items():
        new_gpr = re.sub(rf"\b{old_id}\b", new_id, new_gpr)
        if new_id != "not_found":
            mix_gpr = re.sub(rf"\b{old_id}\b", new_id, mix_gpr)
            globals()[new_id] = symbols(new_id)
        else:
            old_if_formated = "g_" + old_id.replace(".", "_")
            mix_gpr = re.sub(rf"\b{old_id}\b", old_if_formated, mix_gpr)
            globals()[old_if_formated] = symbols(old_if_formated)
    new_gpr = re.sub(r"\bnot_found\b", "1", new_gpr)
    new_gpr = re.sub(r"\bor\b", "+", new_gpr)
    new_gpr = re.sub(r"\band\b", "*", new_gpr)
    mix_gpr = re.sub(r"\bor\b", "+", mix_gpr)
    mix_gpr = re.sub(r"\band\b", "*", mix_gpr)
    equation = sympify(new_gpr)
    equation_expanded = expand(equation)
    equation_mix = sympify(mix_gpr)
    equation_expanded_mix = expand(equation_mix)
    if is_float(str(equation_expanded)):
        new_gpr_expanded = ""
    else:
        new_gpr_expanded = str(equation_expanded).split(" + ")
        for i in range(len(new_gpr_expanded)):
            genes_and = []
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
    new_gpr_expanded_mix = str(equation_expanded_mix).split(" + ")
    for j in range(len(new_gpr_expanded_mix)):
        genes_and_mix = []
        for gm in new_gpr_expanded_mix[j].split("*"):
            if gm and (not is_float(gm)):
                genes_and_mix.append(gm)
        if not genes_and_mix:
            new_gpr_expanded_mix.pop(j)
        elif (len(genes_and_mix) > 1) and (len(new_gpr_expanded_mix) > 1):
            new_gpr_expanded_mix[j] = "(" + " and ".join(sorted(genes_and_mix)) + ")"
        else:
            new_gpr_expanded_mix[j] = " and ".join(sorted(genes_and_mix))
    new_gpr_expanded_mix = " or ".join(sorted(new_gpr_expanded_mix))
    return new_gpr_expanded, new_gpr_expanded_mix


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
