import gzip
import os
import re
import time
import warnings
from functools import partial
from mimetypes import guess_type
from pathlib import PosixPath

import ncbi_genome_download as ngd
import pandas as pd
from cobra import Model
from sympy import expand, symbols, sympify

from .general import is_float


def check_nt_or_aa(path_fasta: PosixPath):
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


def get_genome(ncbi_genome_name: str):
    encoding = guess_type(ncbi_genome_name)[1]  # uses file extension
    _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open

    with _open(ncbi_genome_name) as input_fasta:
        genome_lines = input_fasta.readlines()
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


def get_genes_gapseq(
    output_gene_folder: PosixPath,
    input_gapseq_genome_name: str,
    gapseq_model: Model,
    model_type: str,
    model_id: str,
):
    genomes = get_genome(input_gapseq_genome_name)
    head, tail = os.path.split(input_gapseq_genome_name)
    output_genes_name_gapseq = (
        output_gene_folder
        / f"{os.path.splitext(tail)[0]}_{model_type}_{model_id}_genes.faa"
    )
    gene_gapseq_fasta = open(output_genes_name_gapseq, "w")
    for gene in gapseq_model.genes:
        gene_gapseq_fasta.write(">" + gene.id + "\n")
        start = gene.id.split("_")[-2]
        end = gene.id.split("_")[-1]
        genomeid = gene.id.removeprefix("gp_").removesuffix("_" + start + "_" + end)
        if genomeid not in genomes.keys():
            genomeid = ".".join(genomeid.rsplit("_", 1))
        if genomeid not in genomes.keys():
            return False, False
        start = int(start)
        end = int(end)
        if start < end:
            gene_gapseq_fasta.write(genomes.get(genomeid)[start:end] + "\n")
        else:
            gene_gapseq_fasta.write(genomes.get(genomeid)[end:start] + "\n")
    gene_gapseq_fasta.close()
    aa_status = check_nt_or_aa(output_genes_name_gapseq)
    return output_genes_name_gapseq, aa_status


def get_genes_not_gapseq(
    output_gene_folder: PosixPath,
    input_genes_name: PosixPath,
    model: Model,
    model_type: str,
    model_id: str,
):
    encoding = guess_type(input_genes_name)[1]  # uses file extension
    _open = partial(gzip.open, mode="rt") if encoding == "gzip" else open

    with _open(input_genes_name) as input_fasta:
        lines = input_fasta.readlines()
    head, tail = os.path.split(input_genes_name)
    output_genes_name = (
        output_gene_folder
        / f"{os.path.splitext(tail)[0]}_{model_type}_{model_id}_genes.faa"
    )

    genes_fasta = open(output_genes_name, "w")
    new_id = ""
    for line in lines:
        if line.startswith(">"):
            if new_id != "":
                if new_id in model.genes:
                    genes_fasta.write(">" + new_id + "\n")
                    genes_fasta.write(gene)
            old_id = line.strip().split(" ")[0][1:]
            if model_type == "carveme":
                new_id = "_".join(old_id.rsplit(".", 1))
                new_id = new_id.replace(":", "_")
            elif model_type == "agora":
                new_id = old_id.removeprefix("fig|")
            elif model_type == "modelseed":
                new_id = old_id
            else:
                new_id = old_id
            gene = ""
        else:
            gene = gene + line
    if new_id in model.genes:
        genes_fasta.write(">" + new_id + "\n")
        genes_fasta.write(gene)
    genes_fasta.close()
    aa_status = check_nt_or_aa(output_genes_name)
    return output_genes_name, aa_status


def get_locus_tag_genes(
    ncbi_cds_name: str,
    ncbi_protein_name: str,
    feature_table_name: str,
    out_nt_fasta: str,
    out_aa_fasta: str,
    do_old_locus_tag: bool,
):
    feature_table = pd.read_csv(feature_table_name, sep="\t", compression="gzip")
    out_nt = open(out_nt_fasta, "w")
    with gzip.open(ncbi_cds_name, "rt") as cds:
        locus_tags = []
        for cdsl in cds:
            if cdsl.startswith(">"):
                if locus_tags != []:
                    for lt in locus_tags:
                        out_nt.write(">" + lt + "\n")
                        out_nt.write(one_seq)
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
                    locus_tags = old_locus_tag.split(",")
                else:
                    locus_tags = [new_locus_tag]
                one_seq = ""
            else:
                one_seq = one_seq + cdsl
    for lt in locus_tags:
        out_nt.write(">" + lt + "\n")
        out_nt.write(one_seq)
    out_nt.close()
    out_aa = open(out_aa_fasta, "w")
    with gzip.open(ncbi_protein_name, "rt") as proteins:
        locus_tags = []
        for protl in proteins:
            if protl.startswith(">"):
                if locus_tags != []:
                    for lt in locus_tags:
                        out_aa.write(">" + lt + "\n")
                        out_aa.write(one_seq)
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
                    locus_tags = pold_locus_tag.split(",")
                else:
                    locus_tags = pnew_locus_tag.split(",")
                one_seq = ""
            else:
                one_seq = one_seq + protl
        for lt in locus_tags:
            out_aa.write(">" + lt + "\n")
            out_aa.write(one_seq)
    out_aa.close()


def get_final_fasta_with_ncbi_assemble(
    output_folder: PosixPath, assembly_id: str, do_old_locus_tag=True
):
    gene_path = output_folder / "tmp_gene_conversion"
    path = gene_path / "ncbi_assembly"
    path.mkdir(parents=True, exist_ok=True)
    ngd.download(
        assembly_accessions=assembly_id,
        output=path,
        file_formats="cds-fasta,protein-fasta,features",
        flat_output=True,
    )
    if not os.listdir(path):
        warnings.warn("\nWarning! NCBI download failed. Trying second time")
        time.sleep(15)
        ngd.download(
            assembly_accessions=assembly_id,
            output=path,
            file_formats="cds-fasta,protein-fasta,features",
            flat_output=True,
        )
    if not os.listdir(path):
        warnings.warn("\nWarning! NCBI download failed second time. Trying third time")
        time.sleep(15)
        ngd.download(
            assembly_accessions=assembly_id,
            output=path,
            file_formats="cds-fasta,protein-fasta,features",
            flat_output=True,
        )
    cds_faa = ""
    prot_faa = ""
    feat = ""
    for file in os.listdir(path):
        if file.endswith("_cds_from_genomic.fna.gz"):
            cds_faa = os.path.join(path, file)
        if file.endswith("_protein.faa.gz"):
            prot_faa = os.path.join(path, file)
        if file.endswith("_feature_table.txt.gz"):
            feat = os.path.join(path, file)
    if not cds_faa:
        raise ValueError("CDS fasta file is not found")
    if not prot_faa:
        raise ValueError("Protein fasta file is not found")
    if not feat:
        raise ValueError("Feature table is not found")
    final_nt_faa = os.path.join(output_folder, f"final_nt_{assembly_id}.faa")
    final_aa_faa = os.path.join(output_folder, f"final_aa_{assembly_id}.faa")
    get_locus_tag_genes(
        cds_faa,
        prot_faa,
        feat,
        final_nt_faa,
        final_aa_faa,
        do_old_locus_tag=do_old_locus_tag,
    )
    return final_nt_faa, final_aa_faa


def makeNewGPR(gpr: str, g_id_convert: dict):
    new_gpr = gpr
    mix_gpr = gpr
    for old_id, new_id in g_id_convert.items():
        new_id = new_id.replace(".", "_").replace('"', "").replace(":", "_")
        old_id = old_id.replace(".", "_").replace('"', "").replace(":", "_")
        new_gpr = new_gpr.replace(".", "_").replace('"', "").replace(":", "_")
        mix_gpr = mix_gpr.replace(".", "_").replace('"', "").replace(":", "_")
        if new_id[0].isdigit():
            new_id = "g_" + new_id

        new_gpr = re.sub(rf"\b{old_id}\b", new_id, new_gpr)

        if new_id != "not_found":
            mix_gpr = re.sub(rf"\b{old_id}\b", new_id, mix_gpr)
            globals()[new_id] = symbols(new_id)
        else:
            if old_id[0].isdigit():
                old_if_formated = "g_" + old_id
                mix_gpr = re.sub(rf"\b{old_id}\b", old_if_formated, mix_gpr)
                globals()[old_if_formated] = symbols(old_if_formated)
            else:
                globals()[old_id] = symbols(old_id)
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
    united_gpr = " or ".join(sorted(unique_and))
    return united_gpr
