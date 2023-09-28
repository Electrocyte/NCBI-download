import argparse
import os
import pandas as pd
import glob
from itertools import zip_longest
import shutil


# saw this in the python docs. looks like exactly what you need
def grouper(iterable, n, fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def adjust_q_scores(line: str) -> str:
    return ''.join([chr(ord(char) - 1) for char in line.strip()]) + '\n'


def process_fastq(input_path: str, output_path: str) -> None:
    output_dir = os.path.dirname(output_path)  # Get the directory part of output_path
    os.makedirs(output_dir, exist_ok=True)  # Create the directory if it does not exist

    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for readID, sequence, plus, quality in grouper(infile, 4, None):
            outfile.write(readID)
            outfile.write(sequence)
            outfile.write(plus)
            outfile.write(adjust_q_scores(quality))


def process_samples(input_directory: str, csv_path: str) -> None:
    df = pd.read_csv(csv_path, sep=",")
    for index, row in df.iterrows():
        sample_fields = [row['date'], row['NA'], row['strain'], row['concentration_CFU'], row['batch'], row['duration_h']]
        sample = "_".join(str(field) for field in sample_fields)
        barcode = row['Barcode']
        input_file_paths = glob.glob(f"{input_directory}/{sample}/{barcode}/fastq_pass/*.fastq")
        print(f"Processing {len(input_file_paths)} FASTQ files for sample {sample}...")

        # Store the temporary output file paths
        temp_output_file_paths = []
        for n, input_file_path in enumerate(input_file_paths):
            os.makedirs(f"{input_directory}/{sample}/adjusted/", exist_ok=True)
            temp_output_file_path = f"{input_directory}/{sample}/adjusted/{n}-{barcode}_adjusted-1Q.fastq"
            print(temp_output_file_path)
            process_fastq(input_file_path, temp_output_file_path)
            temp_output_file_paths.append(temp_output_file_path)

        os.makedirs(f"{input_directory}/analysis/sample_data/{sample}/trimmed/", exist_ok=True)
        output_file_path = f"{input_directory}/analysis/sample_data/{sample}/trimmed/{barcode}_adjusted-1Q.fastq"

        # Concatenate all temporary output files to the final output file path
        print(f"Concatenating {len(temp_output_file_paths)} FASTQ files to {output_file_path}...")
        with open(output_file_path, 'wb') as outfile:
            for file_path in temp_output_file_paths:
                with open(file_path, 'rb') as infile:
                    shutil.copyfileobj(infile, outfile)

        # # Optionally, you can remove the temporary files after concatenation.
        # for file_path in temp_output_file_paths:
        #     os.remove(file_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Adjust Q scores in FASTQ files.")
    parser.add_argument("-i", "--input_directory", type=str, required=True, help="Path to the directory containing FASTQ files.")
    parser.add_argument("-c", "--csv_path", type=str, required=True, help="Path to the CSV file containing sample metadata.")
    args = parser.parse_args()

    process_samples(args.input_directory, args.csv_path)


if __name__ == '__main__':
    main()
