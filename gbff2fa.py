import os
import glob
import json
import argparse
import subprocess
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Convert gbff files to fasta from NCBI')

parser.add_argument('-t', '--target', help='Input kingdom target from NCBI')
parser.add_argument('-o', '--directory', help='directory top level')


args = parser.parse_args()

target = args.target
directory = args.directory


# # the root directory where you want to start your search
# directory = "/mnt/e/SequencingData/genbank/"
# target = "protozoa"

root_dir = f"{directory}{target}"

failures_file = f"{directory}{target}-failures.json"

# create output directory if it doesn't exist
output_dir = f"{directory}{target}_fa"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# check if failures.json exists and use it as input
if os.path.exists(failures_file):
    with open(failures_file, 'r') as f:
        failures_dict = json.load(f)
    files = list(failures_dict.keys())
else:
    # walk through the directory structure, finding .gbff.gz and .gbff files
    files = list(glob.iglob(root_dir + '/**/*_genomic.gbff*', recursive=True))

failures = {}

for i, filepath in enumerate(files):
    try:
        # if it's a gzipped file, decompress it
        if filepath.endswith('.gz'):
            try:
                subprocess.check_call(['gunzip', '-k', '-f', filepath])
                filepath = filepath[:-3]  # get the name of the uncompressed file
            except subprocess.CalledProcessError as e:
                print(f"Error decompressing file {filepath}: {e}")
                failures[filepath] = str(e)
                if os.path.exists(filepath):
                    os.remove(filepath)  # remove the corrupt gzipped file
                continue  # skip this file and move to the next one

        # derive output file path
        fa_file = os.path.basename(filepath).replace(".gbff", ".fasta")
        output_file = os.path.join(output_dir, fa_file)

        # check if fasta file already exists
        if os.path.exists(output_file):
            print(f"Fasta file {output_file} already exists. Skipping conversion.")
            continue

        # convert the file
        count = SeqIO.convert(filepath, "genbank", output_file, "fasta")
        print(f"Converted {count} records from {filepath} to {output_file} --- {i/len(files)*100:.2f}% complete")

    except Exception as e:
        print(f"Error with file {filepath}: {e}")
        failures[filepath] = str(e)
        if os.path.exists(filepath):
            if filepath.endswith('.gz'):
               os.remove(filepath)  # remove the corrupt gzipped file

# save failures to json
if os.path.exists(failures_file):
    with open(failures_file.replace("failures", "failures-1"), 'w') as f:
        json.dump(failures, f)
else:
    with open(failures_file, 'w') as f:
        json.dump(failures, f)

# time python3 gbff2fa.py
