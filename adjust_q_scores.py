import argparse
import os


def adjust_q_scores(line: str) -> str:
    return ''.join([chr(ord(char) - 1) for char in line.strip()]) + '\n'


def process_fastq(input_path: str, output_path: str) -> None:
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for index, line in enumerate(infile):
            if (index + 1) % 4 == 0:
                outfile.write(adjust_q_scores(line))
            else:
                outfile.write(line)


def process_fastq_directory(input_directory: str, output_directory: str) -> None:
    for file_name in os.listdir(input_directory):
        if file_name.endswith('.fastq'):
            input_file_path = os.path.join(input_directory, file_name)
            output_file_path = os.path.join(output_directory, file_name)
            process_fastq(input_file_path, output_file_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Adjust Q scores in FASTQ files.")
    parser.add_argument("-i", "--input_directory", type=str, required=True, help="Path to the directory containing FASTQ files.")
    parser.add_argument("-o", "--output_directory", type=str, required=True, help="Path to the directory to save adjusted FASTQ files.")
    args = parser.parse_args()

    process_fastq_directory(args.input_directory, args.output_directory)


if __name__ == '__main__':
    main()
