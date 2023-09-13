import argparse
import os
from Bio import Entrez
import time


def fetch_sequences(email, uid_file, directory, savename):
    # Set the provided email
    Entrez.email = email

    # Ensure the directory exists
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Load UIDs from the file
    uid_file = f"{directory}/{uid_file}"
    print(f"Loading UIDs from {uid_file}...")
    with open(uid_file, "r") as file:
        ids = [line.strip() for line in file.readlines()]

    # Fetch and save sequence for each UID
    for uid in ids:
        if not uid:  # skip if UID is empty
            continue
        time.sleep(1)  # pause for 1 second to be nice to NCBI

        output_path = os.path.join(directory, f"{savename}_{uid}.fasta")

        # Check if file already exists
        if os.path.exists(output_path):
            print(f"File {output_path} already exists. Skipping UID {uid}.")
            continue

        try:
            seq_handle = Entrez.efetch(db="nucleotide", id=uid, rettype="fasta", retmode="text")
        except Exception as e:
            print(f"Error fetching UID {uid}: {e}")
            continue

        print(f"Saving {uid} to {output_path}...")
        with open(output_path, "w") as out_file:
            out_file.write(seq_handle.read())
        seq_handle.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch FASTA sequences from NCBI using a list of UIDs.")
    parser.add_argument("-e", "--email", required=True, help="Your email for NCBI E-utilities access.")
    parser.add_argument("-f", "--file", required=True, dest="uid_file", help="Path to the file containing the list of UIDs.")
    parser.add_argument("-d", "--directory", required=True, help="Directory to save the FASTA files.")
    parser.add_argument("-s", "--savename", required=True, help="Batch ID.")
    args = parser.parse_args()

    fetch_sequences(args.email, args.uid_file, args.directory, args.savename)
