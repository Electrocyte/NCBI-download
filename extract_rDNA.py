#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 08:36:22 2023

@author: mangi
"""

from Bio import SeqIO
import glob
import re
import os
import argparse

def reverse_complement(seq):
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'I': 'I',  # Inosine remains as 'I' in the reverse complement
        # Include other IUPAC codes if needed
        'R': 'Y',
        'Y': 'R',
        'S': 'S',
        'W': 'W',
        'K': 'M',
        'M': 'K',
        'B': 'V',
        'D': 'H',
        'H': 'D',
        'V': 'B',
        'N': 'N'
    }
    try:
        return "".join(complement[base] for base in reversed(seq))
    except KeyError:
        print(f"Error: The base '{seq}' is not recognized.")
        return ""

iupac_dict = {
    'R': '[AG]',
    'Y': '[CT]',
    'S': '[GC]',
    'W': '[AT]',
    'K': '[GT]',
    'M': '[AC]',
    'B': '[CGT]',
    'D': '[AGT]',
    'H': '[ACT]',
    'V': '[ACG]',
    'N': '[ACGT]',
    'I': '[ACGT]'  # Inosine can pair with any of the canonical bases
}


gene_config = {
    "_16S": {
        "primerFnorm": "AGRGTTTGATCMTGGCTCAG",
        "primerRnorm": "GACATCGAGGTGCCAAAC",
        "save_folder": "BAC",
        "_type_": "16S-ITS-23S"
        # ... add other configurations if required
    },
    "rpoB_KA": {
        "primerFnorm": "CGATCCGAAGGACAACCTGTT",
        "primerRnorm": "TCRTCRTAIGGCATRTCYTC",
        "save_folder": "BAC-rpoB",
        "_type_": "rpoB-bac"
        # ... add other configurations if required
    },
    "rpoB_OA": {
        "primerFnorm": "GGYTWYGAAGTNCGHGACGTDCA",
        "primerRnorm": "TCRTCRTAIGGCATRTCYTC",
        "save_folder": "BAC-rpoB-OA",
        "_type_": "rpoB-OA-bac"
        # ... add other configurations if required
    },
    "rpoB_A": {
        "primerFnorm": "GCITTYATGCCITGGAAYGG",
        "primerRnorm": "TCRTCRTAIGGCATRTCYTC",
        "save_folder": "BAC-rpoB-A",
        "_type_": "rpoB-A-bac"
        # ... add other configurations if required
    },
    "cpn60": {
        "primerFnorm": "GAIIIIGCIGGIGAYGGIACIACIAC",
        "primerRnorm": "YKIYKITCICCRAAICCIGGIGCYTT",
        "save_folder": "BAC-cpn60",
        "_type_": "cpn60-bac"
        # ... add other configurations if required
    }
}


def setup_gene_config(gene, iupac_dict, out_head, fasta_directory):
    config = gene_config[gene]

    primerFnorm = config["primerFnorm"]
    primerRnorm = config["primerRnorm"]

    # Store original values
    primerFnorm_ori = primerFnorm
    primerRnorm_ori = primerRnorm

    # Convert IUPAC codes to regular expressions
    for key, value in iupac_dict.items():
        primerFnorm = primerFnorm.replace(key, value)
        primerRnorm = primerRnorm.replace(key, value)

    # Calculate reverse complement
    primerFrev = reverse_complement(primerFnorm_ori)
    primerRrev = reverse_complement(primerRnorm_ori)

    # Convert IUPAC codes in the reverse complements
    for key, value in iupac_dict.items():
        primerFrev = primerFrev.replace(key, value)
        primerRrev = primerRrev.replace(key, value)

    # Set up locations
    save_loc = f"{out_head}/operons/library/{config['save_folder']}/"
    fasta_file_paths = f"{fasta_directory}library/bacteria/*.fna"
    seqid_loc = f"{out_head}/operons/operon-{config['_type_']}-NCBI-seqids.csv"

    return primerFnorm, primerRnorm, primerFrev, primerRrev, save_loc, fasta_file_paths, seqid_loc


def write_sequence_to_file(sequence_name, _type_, extracted_sequence, filename):
    """Writes the provided sequence to a file.

    Args:
        sequence_name (str): The name of the sequence.
        _type_ (str): The type of the sequence.
        extracted_sequence (str): The sequence to be written.
        filename (str): Path to the file where the sequence will be written.

    Returns:
        int: Returns 1 if the sequence is successfully written, 0 otherwise.
    """
    print(filename)
    with open(filename, "w") as f:
        f.write(">" + sequence_name + _type_ + "\n" + extracted_sequence)
    return 1


def main(args):
    location = args.location
    local = args.local
    out_head = args.out_head

    bacterial = args.bacterial
    fungal = args.fungal

    _16S = args.r16S
    rpoB_KA = args.rpoB_KA
    rpoB_OA = args.rpoB_OA
    rpoB_A = args.rpoB_A
    cpn60 = args.cpn60

    fungal_centrifuge = args.fungal_centrifuge
    fungal_JGI = args.fungal_JGI
    protozoa_genbank = args.protozoa_genbank
    fungal_genbank = args.fungal_genbank

    if bacterial:
        if _16S or rpoB_KA or rpoB_OA or cpn60 or rpoB_A:
            fasta_directory = f"{location}/bacteria/"

    if fungal:
        if fungal_centrifuge:
            fasta_directory = f"{location}/any_fungus/centrifuge/"
        if fungal_JGI:
            _folder_ = f"{local}/JGI/fungal/"
        if protozoa_genbank:
            _folder_ = f"{local}/genbank/protozoa/"
        if fungal_genbank:
            _folder_ = f"{local}/genbank/fungal/"

    if fungal:
        # fungal
        primerFnorm = "CAGTAGTCATATGCTTGTC"  # 18S NS1 short F
        primerRnorm = "CTATGTTTTAATTAGACAGTCAG"  # 28S RCA95m R
        primerFrev = "GACAAGCATATGACTACTG" # rev comp
        primerRrev = "CTGACTGTCTAATTAAAACATAG" # rev comp
        _type_ = "18S-ITS-28S"

        if fungal_centrifuge:
            save_loc = f"{out_head}/operons/library/EUK/"
            fasta_file_paths = f"{fasta_directory}library/fungi/*.fna"
            seqid_loc = f"{out_head}/operons/operon-{_type_}-NCBI-seqids.csv"
            multiline = False

        if fungal_JGI:
            save_loc = f"{_folder_}/EUK/"
            fasta_file_paths = f"{_folder_}/renamed/*.fasta"
            seqid_loc = f"{_folder_}/operon-{_type_}-JGI-seqids.csv"
            multiline = True

        if protozoa_genbank:
            save_loc = f"{_folder_}/EUK-PROTO/"
            fasta_file_paths = f"{_folder_}/protozoa_fa/*.fasta"
            seqid_loc = f"{_folder_}/protozoa-operon-{_type_}-NCBI-gbff-seqids.csv"
            multiline = False

        if fungal_genbank:
            save_loc = f"{_folder_}/EUK-FUNGI/"
            fasta_file_paths = f"{_folder_}/fungi_fa/*.fasta"
            seqid_loc = f"{_folder_}/fungi-operon-{_type_}-NCBI-gbff-seqids.csv"
            multiline = False

    if bacterial:
        if _16S:
            primerFnorm, primerRnorm, primerFrev, primerRrev, save_loc, fasta_file_paths, seqid_loc = setup_gene_config("_16S", iupac_dict, out_head, fasta_directory)

        if rpoB_KA:
            primerFnorm, primerRnorm, primerFrev, primerRrev, save_loc, fasta_file_paths, seqid_loc = setup_gene_config("rpoB_KA", iupac_dict, out_head, fasta_directory)

        if rpoB_OA:
            primerFnorm, primerRnorm, primerFrev, primerRrev, save_loc, fasta_file_paths, seqid_loc = setup_gene_config("rpoB_OA", iupac_dict, out_head, fasta_directory)

        if rpoB_A:
            primerFnorm, primerRnorm, primerFrev, primerRrev, save_loc, fasta_file_paths, seqid_loc = setup_gene_config("rpoB_A", iupac_dict, out_head, fasta_directory)

        if cpn60:
            primerFnorm, primerRnorm, primerFrev, primerRrev, save_loc, fasta_file_paths, seqid_loc = setup_gene_config("cpn60", iupac_dict, out_head, fasta_directory)

    complement_pair = (primerFnorm, primerRrev)
    rev_complement_pair = (primerFrev, primerRnorm)

    pairs = {"complement": complement_pair, "revcomp": rev_complement_pair}

    count = 0

    for name, p in pairs.items():

        primerF = p[0]
        primerR = p[1]

        # Convert primer sequences to regex pattern
        pattern = primerF + "(.*)" + primerR

        # Define input fasta file path

        print(fasta_file_paths)

        fasta_file_globs = glob.glob(fasta_file_paths)
        print(len(fasta_file_globs))

        os.makedirs(save_loc, exist_ok = True)

        for fasta_file_path in fasta_file_globs:

            # Parse fasta file
            with open(fasta_file_path, "r") as handle:
                fasta_sequences = list(SeqIO.parse(handle, 'fasta'))


            # Loop through fasta sequences
            for i, fasta in enumerate(fasta_sequences):
                sequence = str(fasta.seq)
                sequence_name = "_".join(fasta.description.split()[0:4])  # Take the first word of the fasta header as the sequence name
                # Use regex to remove punctuation from string
                sequence_name = re.sub(r'(?<!\d)[.,;](?!\d)', '', sequence_name)
                sequence_name = sequence_name.replace("/", "")

                if multiline:
                    base_name = os.path.basename(fasta_file_path)
                    sequence_name = base_name.split("_Assembly")[0]
                    match = re.search(pattern, sequence, re.DOTALL)
                else:
                    match = re.search(pattern, sequence)

                if match:
                    # Extract sequence
                    extracted_sequence = match.group()
                    print(fasta_file_path, sequence_name)
                    print(match, len(extracted_sequence))

                    filename = f"{save_loc}{sequence_name}_{_type_}_{name}.fasta"

                    # Decision to write to file based on sequence length and type
                    should_write = False
                    if _type_ == "cpn60-bac" and 450 <= len(extracted_sequence) <= 1000:
                        should_write = True
                    elif _type_ == "rpoB-OA-bac" and 800 <= len(extracted_sequence) <= 3000:
                        should_write = True
                    elif _type_ == "rpoB-A-bac" and 600 <= len(extracted_sequence) <= 1200:
                        should_write = True
                    elif 1000 <= len(extracted_sequence) <= 10000:
                        should_write = True

                    print(should_write, _type_)
                    if should_write:
                        count += write_sequence_to_file(sequence_name, _type_, extracted_sequence, filename)
                    else:
                        pattern = primerF + "(.*?)" + primerR

                        matches = re.findall(pattern, extracted_sequence)
                        minimized_matches = []
                        print()
                        print(fasta_file_path)
                        print(sequence_name)
                        for match in matches:
                            no_n = match.count('N')
                            print(len(match), no_n, len(match)%3)
                            # if len(match) < 10000:
                            minimized_matches.append((match, no_n, len(match)%3))

                        # Filter the list to keep only tuples where the first item is less than 10000
                        filtered_matches = [match for match in minimized_matches if len(match[0]) < 10000]

                        # If filtered list is not empty, find the minimum
                        if filtered_matches:
                            best_match = min(filtered_matches, key=lambda x: (x[1], x[2]))
                            extracted_sequence = best_match[0]
                            print(f"Found rDNA length: {len(extracted_sequence)}")
                            if 1000 <= len(extracted_sequence) <= 10000:
                                count += 1
                                # Write to file
                                print(filename)
                                with open(filename, "w") as f:
                                    f.write(">" + sequence_name + "-" + _type_ + "\n" + extracted_sequence)
                            print()
                        else:
                            print("No matches found with first element less than 10000.")

    original_file_count = count / len(fasta_file_globs) * 100
    print(f"Number of rDNA operon regions found: {count}, {original_file_count:.2f}; {count/original_file_count:.2f}")

    gen_files = glob.glob(f"{save_loc}*a")

    sseqids = []
    for gen_file in gen_files:

        split = gen_file.split("/")[-1]
        print(split)
        if multiline:
            sseqid = "_".join(split.split("_")[2:3])
            name = " ".join(split.split(".")[0].split("_")[0:2]) + " " + " ".join(split.split(".")[-2].split("_")[-2:])
        else:
            sseqid = "_".join(split.split("_")[0:2])
            name = " ".join(split.split(".")[1].split("_")[1:])
        sseqids.append([sseqid, name])

    with open(seqid_loc, "w") as f_out:
        print(seqid_loc)
        f_out.write("code,name\n")
        f_out.write("\n".join([f"{i},{j}" for i, j in sseqids]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract rDNA Operon Regions")

    parser.add_argument("-l", "--location", required=True, help="Location of sequencing data Centrifuge libraries.")
    parser.add_argument("-z", "--local", required=True, help="Local path of sequencing data.")
    parser.add_argument("-o", "--out_head", required=True, help="Specify the save location.")

    parser.add_argument("-b", "--bacterial", action="store_true", help="Specify if bacterial.")
    parser.add_argument("-s", "--r16S", action="store_true", help="Specify if 16S-23S.")
    parser.add_argument("-r", "--rpoB_KA", action="store_true", help="Specify if rpoB KA.")
    parser.add_argument("-a", "--rpoB_OA", action="store_true", help="Specify if rpoB OA.")
    parser.add_argument("-t", "--rpoB_A", action="store_true", help="Specify if rpoB A.")
    parser.add_argument("-n", "--cpn60", action="store_true", help="Specify if cpn60.")

    parser.add_argument("-f", "--fungal", action="store_true", help="Specify if fungal.")
    parser.add_argument("-c", "--fungal_centrifuge", action="store_true", help="Specify if using fungal centrifuge.")
    parser.add_argument("-j", "--fungal_JGI", action="store_true", help="Specify if using fungal JGI.")
    parser.add_argument("-p", "--protozoa_genbank", action="store_true", help="Specify if using protozoa genbank.")
    parser.add_argument("-g", "--fungal_genbank", action="store_true", help="Specify if using fungal genbank.")

    args = parser.parse_args()
    main(args)
