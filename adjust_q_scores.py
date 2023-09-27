import argparse
import os
import pandas as pd
import glob


def adjust_q_scores(line: str) -> str:
    return ''.join([chr(ord(char) - 1) for char in line.strip()]) + '\n'


def process_fastq(input_path: str, output_path: str) -> None:
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for index, line in enumerate(infile):
            if (index + 1) % 4 == 0:
                outfile.write(adjust_q_scores(line))
            else:
                outfile.write(line)


def process_samples(input_directory: str, csv_path: str) -> None:
    df = pd.read_csv(csv_path, sep="\t")
    for index, row in df.iterrows():
        sample_fields = [row['date'], row['NA'], row['strain'], row['concentration_CFU'], row['batch'], row['duration_h']]
        sample = "_".join(str(field) for field in sample_fields)
        barcode = row['Barcode']
        input_file_paths = glob.glob(f"{input_directory}/{sample}/{barcode}/fastq_pass/*.fastq")
        print(f"Processing {len(input_file_paths)} FASTQ files for sample {sample}...")
        for input_file_path in input_file_paths:
            output_file_path = f"{input_directory}/analysis/sample_data/{sample}/trimmed/adjusted-1Q.fastq"
            process_fastq(input_file_path, output_file_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Adjust Q scores in FASTQ files.")
    parser.add_argument("-i", "--input_directory", type=str, required=True, help="Path to the directory containing FASTQ files.")
    parser.add_argument("-c", "--csv_path", type=str, required=True, help="Path to the CSV file containing sample metadata.")
    args = parser.parse_args()

    process_samples(args.input_directory, args.csv_path)


if __name__ == '__main__':
    main()
