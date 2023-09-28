import subprocess
import argparse
import os

def run_nanofilt(quality: int, input_file: str, output_directory: str) -> None:
    """
    Runs NanoFilt on the input file with the specified quality,
    and writes the output to the output directory.

    :param quality: Quality threshold for NanoFilt.
    :param input_file: Path to the input FASTQ file.
    :param output_directory: Path to the output directory.
    """
    output_file = os.path.join(output_directory, f"q{quality}.fastq")

    # Run the NanoFilt command
    print(output_file)
    with open(output_file, 'w') as outfile:
        subprocess.run(['NanoFilt', '-q', str(quality), input_file], stdout=outfile)


def process_files(directory: str) -> None:
    """
    Process all FASTQ files in the specified directory with NanoFilt.

    :param directory: Path to the directory containing FASTQ files.
    """
    # Loop through all subdirectories and files in the directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.startswith('trimmed') and file.endswith('.fastq'):
                file_path = os.path.join(root, file)

                # Set output_directory as the current root directory
                output_directory = root
                os.makedirs(output_directory, exist_ok=True)

                # Run NanoFilt with quality 14 and 17
                run_nanofilt(14, file_path, output_directory)
                run_nanofilt(17, file_path, output_directory)


def main():
    parser = argparse.ArgumentParser(description="Process FASTQ files with NanoFilt.")
    parser.add_argument("-d", "--directory", type=str, help="Path to the directory containing FASTQ files.")

    args = parser.parse_args()
    process_files(args.directory)


if __name__ == '__main__':
    main()
