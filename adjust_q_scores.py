import argparse
import os
import pandas as pd
import glob


def adjust_q_scores(line: str) -> str:
    return ''.join([chr(ord(char) - 1) for char in line.strip()]) + '\n'


def process_fastq(input_path: str, output_path: str) -> None:
    output_dir = os.path.dirname(output_path)  # Get the directory part of output_path
    os.makedirs(output_dir, exist_ok=True)  # Create the directory if it does not exist

    with open(input_path, 'r') as infile, open(output_path, 'a') as outfile:
        for index, line in enumerate(infile):
            if (index + 1) % 4 == 0:
                outfile.write(adjust_q_scores(line))
            else:
                outfile.write(line)


def process_samples(input_directory: str, csv_path: str) -> None:
    df = pd.read_csv(csv_path, sep=",")
    for index, row in df.iterrows():
        sample_fields = [row['date'], row['NA'], row['strain'], row['concentration_CFU'], row['batch'], row['duration_h']]
        sample = "_".join(str(field) for field in sample_fields)
        barcode = row['Barcode']
        input_file_paths = glob.glob(f"{input_directory}/{sample}/{barcode}/fastq_pass/*.fastq")
        print(f"Processing {len(input_file_paths)} FASTQ files for sample {sample}...")
        for input_file_path in input_file_paths:
            os.makedirs(f"{input_directory}/analysis/", exist_ok=True)
            os.makedirs(f"{input_directory}/analysis/sample_data/", exist_ok=True)
            os.makedirs(f"{input_directory}/analysis/sample_data/{sample}/", exist_ok=True)
            os.makedirs(f"{input_directory}/analysis/sample_data/{sample}/trimmed/", exist_ok=True)
            output_file_path = f"{input_directory}/analysis/sample_data/{sample}/trimmed/{barcode}_adjusted-1Q.fastq"
            # If it's the first file, ensure the output file does not exist
            if input_file_path == input_file_paths[0] and os.path.exists(output_file_path):
                os.remove(output_file_path)
            process_fastq(input_file_path, output_file_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Adjust Q scores in FASTQ files.")
    parser.add_argument("-i", "--input_directory", type=str, required=True, help="Path to the directory containing FASTQ files.")
    parser.add_argument("-c", "--csv_path", type=str, required=True, help="Path to the CSV file containing sample metadata.")
    args = parser.parse_args()

    process_samples(args.input_directory, args.csv_path)


if __name__ == '__main__':
    main()
