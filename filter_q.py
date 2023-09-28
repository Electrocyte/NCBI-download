import subprocess
import os

def run_nanofilt(quality: int, input_file: str, output_file: str) -> None:
    """
    Runs NanoFilt on the input file with the specified quality,
    and writes the output to the output file.

    :param quality: Quality threshold for NanoFilt.
    :param input_file: Path to the input FASTQ file.
    :param output_file: Path to the output FASTQ file.
    """
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
                base, ext = os.path.splitext(file)

                # Run NanoFilt with quality 14 and 17
                run_nanofilt(14, file_path, f"{base}_q14{ext}")
                run_nanofilt(17, file_path, f"{base}_q17{ext}")


# Set the directory containing your FASTQ files
dir_path = '/mnt/usersData/DNA_uHQ/analysis/sample_data'

# Process all files in the directory
process_files(dir_path)
