import argparse
import os
from Bio.Align.Applications import ClustalOmegaCommandline


def run_msa(input_dir, savename):
    # Combining all fasta files into one for MSA
    combined_file_path = os.path.join(input_dir, "pcna_combined.fasta")
    output_file = os.path.join(input_dir, f"{savename}-clustal.aln")

    # Check if combined file already exists to avoid overwriting
    if not os.path.exists(combined_file_path):
        with open(combined_file_path, "w") as outfile:
            for filename in os.listdir(input_dir):
                if filename.endswith(".fasta"):
                    with open(os.path.join(input_dir, filename), "r") as infile:
                        outfile.write(infile.read())

    # Running MSA using Clustal Omega
    try:
        clustalomega_cline = ClustalOmegaCommandline(infile=combined_file_path, outfile=output_file, verbose=True, auto=True)
        stdout, stderr = clustalomega_cline()
    except Exception as e:
        print(f"Error running Clustal Omega: {e}")
        return

    print(f"MSA completed. Check the output in {output_file}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform MSA on FASTA sequences.")

    parser.add_argument("-d", "--directory", required=True, help="Directory to save the FASTA files.")
    parser.add_argument("-s", "--savename", required=True, help="Batch ID.")
    args = parser.parse_args()

    run_msa(args.directory, args.savename)
