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

#### EDIT HERE ####
bacterial = False
if bacterial:
    fasta_directory = "SequencingData/Centrifuge_libraries/bacteria/"
    out_head = "" # e.g. /mnt/e/

fungal = True

fungal_centrifuge = False
if fungal_centrifuge:
    fasta_directory = "SequencingData/Centrifuge_libraries/any_fungus/"
    out_head = "" # e.g. /home/yourname

fungal_JGI = False
if fungal_JGI:
    _folder_ = "/mnt/e/SequencingData/JGI/"

protozoa_genbank = False
if protozoa_genbank:
    _folder_ = "/mnt/e/SequencingData/genbank/"

fungal_genbank = True
if fungal_genbank:
    _folder_ = "/mnt/e/SequencingData/genbank/"
#### EDIT HERE ####

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
    # bacterial
    primerFnorm = "AGRGTTTGATCMTGGCTCAG"  # 16S 27F expanded
    # Convert primer sequences to regex pattern
    primerFnorm = primerFnorm.replace("R", "[AG]").replace("M", "[AC]")
    primerRnorm = "GACATCGAGGTGCCAAAC"  # 23S 2490R
    primerFrev = "CTGAGCCAKGATCAAACYCT" # rev comp
    # Convert primer sequences to regex pattern
    primerFrev = primerFrev.replace("Y", "[TC]").replace("K", "[GT]")
    primerRrev = "GTTTGGCACCTCGATGTC" # rev comp

    _type_ = "16S-ITS-23S"
    save_loc = "{out_head}/operons/library/BAC/"
    fasta_file_paths = f"{fasta_directory}library/bacteria/*.fna"
    seqid_loc = f"{out_head}/operons/operon-{_type_}-NCBI-seqids.csv"
    multiline = False


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
        fasta_sequences = SeqIO.parse(open(fasta_file_path), 'fasta')

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

                # Check if the length of the sequence is between 1000 and 10000
                if 1000 <= len(extracted_sequence) <= 10000:
                    count += 1
                    # Write to file
                    print(filename)
                    with open(filename, "w") as f:
                        f.write(">" + sequence_name + _type_ + "\n" + extracted_sequence)
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

# time python extract_rDNA.py
