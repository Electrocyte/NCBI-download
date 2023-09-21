import argparse
import os
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from typing import List
from Bio.Application import ApplicationError


def concatenate_and_split(files: List[str], target_size: int, prefix: str) -> List[str]:
    """
    Concatenate sequences from multiple files until the file reaches a target size.
    """
    batch_count = 1
    current_size = 0
    current_content = []
    output_files = []

    for filepath in files:
        with open(filepath, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                record_size = len(record.format("fasta").encode('utf-8'))

                if (current_size + record_size) > target_size * 1024:
                    output_path = f"{prefix}-{batch_count}.fasta"
                    print(f"Writing concatenated sequences to {output_path}...")
                    with open(output_path, "w") as output_file:
                        SeqIO.write(current_content, output_file, "fasta")
                    output_files.append(output_path)

                    batch_count += 1
                    current_size = 0  # Reset current size
                    current_content = []  # Reset current content

                current_content.append(record)
                current_size += record_size

    # Write the remaining sequences if there's any left
    if current_content:
        output_path = f"{prefix}-{batch_count}.fasta"
        print(f"Writing concatenated sequences to {output_path}...")
        with open(output_path, "w") as output_file:
            SeqIO.write(current_content, output_file, "fasta")
        output_files.append(output_path)

    return output_files



def run_msa(input_dir: str, savename: str, batch_size:int=50) -> None:
    filenames = [f for f in os.listdir(f"{input_dir}/{savename}") if f.endswith(".fasta")]

    combined_folder = os.path.join(input_dir, "cat_files")
    os.makedirs(combined_folder, exist_ok=True)

    num_batches = (len(filenames) + batch_size - 1) // batch_size

    for batch_num in range(num_batches):
        print(f"Processing batch {batch_num+1}/{num_batches}...")

        start_idx = batch_num * batch_size
        end_idx = start_idx + batch_size
        batch_files = filenames[start_idx:end_idx]

        combined_file_path = os.path.join(combined_folder, f"{savename}_batch{batch_num}.fasta")

        with open(combined_file_path, "w") as outfile:
            for filename in batch_files:
                with open(os.path.join(f"{input_dir}/{savename}", filename), "r") as infile:
                    outfile.write(infile.read())

        # Check the combined file size
        file_size_kb = os.path.getsize(combined_file_path) / 1024

        # If it exceeds 90KB, split the data using concatenate_and_split
        if file_size_kb > 90:
            batch_files_paths = [os.path.join(f"{input_dir}/{savename}", filename) for filename in batch_files]
            split_files = concatenate_and_split(batch_files_paths, 90, combined_file_path.replace('.fasta', ''))
            print(split_files)

            for n, split_file in enumerate(split_files):
                # Check the number of sequences in the split_file
                num_sequences = sum(1 for _ in SeqIO.parse(split_file, "fasta"))
                if num_sequences < 2:
                    print(f"Skipping alignment for {split_file} as it contains less than two sequences.")
                    continue

                output_file = split_file.replace(".fasta", f"-{n}-clustal.aln")
                output_file = output_file.replace("cat_files/", "")

                if os.path.exists(output_file):
                    print(f"Output file {output_file} already exists. Skipping this batch...")
                    continue

                clustalomega_cline = ClustalOmegaCommandline(infile=split_file, outfile=output_file, verbose=True, auto=True)
                try:
                    stdout, stderr = clustalomega_cline()
                    print(f"MSA completed for {split_file}. Check the output in {output_file}.")
                except ApplicationError as e:
                    print(f"Error running Clustal Omega for {split_file}: {e}")
        else:
            # If it doesn't exceed 90KB, simply run Clustal Omega for the combined file
            output_file = os.path.join(input_dir, f"{savename}-batch{batch_num}-clustal.aln")

            if os.path.exists(output_file):
                print(f"Output file {output_file} already exists. Skipping this batch...")
                continue

            clustalomega_cline = ClustalOmegaCommandline(infile=combined_file_path, outfile=output_file, verbose=True, auto=True)
            stdout, stderr = clustalomega_cline()
            print(f"MSA completed for batch {batch_num+1}/{num_batches}. Check the output in {output_file}.")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform MSA on FASTA sequences.")

    parser.add_argument("-d", "--directory", required=True, help="Directory to save the FASTA files.")
    parser.add_argument("-s", "--savename", required=True, help="Batch ID.")
    args = parser.parse_args()

    run_msa(args.directory, args.savename)
