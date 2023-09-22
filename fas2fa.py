#!/usr/bin/env python3

import os
import re
import argparse
from typing import List

def fasta_splitter(input_file: str, output_dir: str) -> None:
    with open(input_file, 'r') as f:
        content = f.read().strip()
        sequences = content.split('>')[1:]  # Split by '>' but skip the first empty item

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for sequence in sequences:
        lines = sequence.split('\n')
        header = lines[0]
        # Clean up punctuation, but leave underscores
        clean_header = re.sub(r'[^\w\s_]', '', header).replace(' ', '_')
        seq = '\n'.join(lines[1:])

        output_file = os.path.join(output_dir, clean_header + '.fasta')
        with open(output_file, 'w') as f:
            f.write('>' + header + '\n' + seq)

def main():
    parser = argparse.ArgumentParser(description="Splits a FASTA file with multiple headers into individual FASTA files.")
    parser.add_argument("-f", "--input_filepath", help="Path to the input FASTA file.")
    parser.add_argument("-o", "--output_directory", help="Path to the output directory where individual FASTA files will be saved.")

    args = parser.parse_args()

    fasta_splitter(args.input_filepath, args.output_directory)

if __name__ == '__main__':
    main()
