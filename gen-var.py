#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 10:36:00 2023

@author: mangi
"""

import os
import glob
import subprocess
import pandas as pd
from Bio import SeqIO
from typing import Dict, List, Union


def generate_tetranucleotides() -> List[str]:
    """Generate all possible tetranucleotides."""
    bases = ['A', 'C', 'G', 'T']
    return [a+b+c+d for a in bases for b in bases for c in bases for d in bases]


def get_complimentary_tetra(tetra: str) -> str:
    """Return the complimentary tetranucleotide."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([comp[base] for base in reversed(tetra)])


def tetranucleotide_frequencies(seq: str) -> Dict[str, int]:
    """Calculate tetranucleotide frequencies in a sequence."""
    freqs = {tetra: 0 for tetra in generate_tetranucleotides()}
    for i in range(0, len(seq) - 3):
        tetra = seq[i:i+4]
        freqs[tetra] += 1
    return freqs


def calculate_variance(freqs: Dict[str, int]) -> float:
    """Calculate variance of complimentary tetranucleotide frequencies."""
    values = []
    for tetra, freq in freqs.items():
        comp = get_complimentary_tetra(tetra)
        combined_freq = freq + freqs[comp]
        values.append(combined_freq)

    mean = sum(values) / len(values) if len(values) > 0 else 0

    # If mean is zero, return zero to avoid ZeroDivisionError
    if mean == 0:
        return 0

    s_squared = sum([(x - mean)**2 for x in values]) / len(values)
    return s_squared / mean


def calculate_gc_content(seq: str) -> float:
    """Calculate GC content of a sequence."""
    g_count = seq.count('G')
    c_count = seq.count('C')
    return (g_count + c_count) / len(seq) * 100  # return as a percentage


def get_total_reads(file_path: str) -> int:
    """Get the total number of reads in the file (FASTA or FASTQ)."""
    result = subprocess.run(["wc", "-l", file_path], capture_output=True, text=True)
    total_lines = int(result.stdout.split()[0])
    if file_path.endswith(".fastq") or file_path.endswith(".fq"):
        return total_lines // 4
    elif file_path.endswith(".fasta") or file_path.endswith(".fa"):
        return total_lines // 2
    else:
        raise ValueError("Unknown file format. Expected .fasta or .fastq.")


def main(file_path: str) -> List[List[Union[str, float, int]]]:
    """Main function."""
    fname = " ".join(file_path.split('/')[-1].split('.')[0].split('-')[0:3]).replace(' BAC', '')
    # print(fname)

    total_reads = get_total_reads(file_path)
    # print(f"Total reads: {total_reads}")

    # Detect file format for SeqIO parsing
    if file_path.endswith(".fastq") or file_path.endswith(".fq"):
        file_format = "fastq"
    elif file_path.endswith(".fasta") or file_path.endswith(".fa"):
        file_format = "fasta"
    else:
        raise ValueError("Unknown file format. Expected .fasta or .fastq.")

    n=0
    list_var = []
    for record in SeqIO.parse(file_path, file_format):
        dna_seq = str(record.seq)
        variance = calculate_variance(tetranucleotide_frequencies(dna_seq))
        gc_content = calculate_gc_content(dna_seq)
        n+=1
        list_var.append([fname, record.id, variance, gc_content, total_reads])
        # print(f"Read number: {n}; Read ID: {record.id}, Variance: {variance:.4f}, GC Content: {gc_content:.2f}%")

    print(f"Taxon: {fname}, Total reads: {total_reads}")
    return list_var


directory_path = '/mnt/d/SequencingData/Harmonisation/DNA/analysis/sample_data/20230823_aDNA_Pa-9027-16S-23S-2X-non-bisulfite-primers-5_1000CFU_924_12/trimmed/no*q'
fq = glob.glob(directory_path)
fq = []

# New path for FASTA files
fasta_directory_path = '/mnt/d/SequencingData/Harmonisation/DNA/analysis/sample_data/20230823_aDNA_Pa-9027-16S-23S-2X-non-bisulfite-primers-5_1000CFU_924_12/trimmed/minimap2-no-host/*Centrifuge.fasta'
fa = glob.glob(fasta_directory_path)

# Using the function
# filtered_fa = [file for file in fa if 'BAC-operons' in file]
# print(filtered_fa)

list_var = []

if len(fq) > 0:
    for i in fq:
        dict_var = main(i)
        list_var.append(dict_var)

elif len(fa) > 0:
    for j in fa:
        # print(j)
        dict_var = main(j)
        list_var.append(dict_var)

# Convert list of dictionaries to DataFrame
df = pd.DataFrame([item for sublist in list_var for item in sublist],
                  columns=['Taxon', 'Read ID', 'Variance', 'GC Content', 'Total Reads'])

# Check if var folder exists and create it if it doesn't
var_folder = "/mnt/d/SequencingData/Harmonisation/DNA/analysis/sample_data/20230823_aDNA_Pa-9027-16S-23S-2X-non-bisulfite-primers-5_1000CFU_924_12/var"
if not os.path.exists(var_folder):
    os.makedirs(var_folder)

# Save DataFrame to CSV file in var folder
csv_file = os.path.join(var_folder, "variances.csv")
df.to_csv(csv_file, index=False)

# Print success message
print(f"CSV file saved to {csv_file}")
