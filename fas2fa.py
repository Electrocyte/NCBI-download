#!/usr/bin/env python3

import os
import re
import argparse
from typing import List

def fasta_splitter(input_file: str, output_dir: str, limit_filename_length: bool = False) -> None:
    with open(input_file, 'r') as f:
        content = f.read().strip()
        sequences = content.split('>')[1:]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for sequence in sequences:
        lines = sequence.split('\n')
        header = lines[0]
        clean_header = re.sub(r'[^\w\s_]', '', header).replace(' ', '_')

        # If filename length should be limited
        if limit_filename_length:
            # Take 10 first and 20 last characters
            clean_header = clean_header[:10] + clean_header[-30:]

        seq = '\n'.join(lines[1:])
        output_file = os.path.join(output_dir, clean_header + '.fasta')
        with open(output_file, 'w') as f:
            f.write('>' + header + '\n' + seq)

def main():
    parser = argparse.ArgumentParser(description="Splits a FASTA file into individual FASTA files.")
    parser.add_argument("-f", "--input_filepath", help="Path to the input FASTA file.")
    parser.add_argument("-o", "--output_directory", help="Path to the output directory.")
    parser.add_argument("-l", "--limit", action='store_true', help="Limit the length of the output file names.")

    args = parser.parse_args()

    fasta_splitter(args.input_filepath, args.output_directory, args.limit)

if __name__ == '__main__':
    main()
