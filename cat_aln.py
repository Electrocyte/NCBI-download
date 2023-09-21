import glob
import argparse
from Bio import AlignIO
# from Bio.Align import MultipleSeqAlignment
# from Bio.SeqRecord import SeqRecord
from typing import List
import re

def extract_consecutive_regions(filename: str) -> None:
    # Generate all permutations of "ATCG" repeated 4 times
    permutations = [
        a + b + c + d
        for a in "ATCG"
        for b in "ATCG"
        for c in "ATCG"
        for d in "ATCG"
    ]

    # Create a regex pattern to search for any of the permutations
    pattern = "|".join(permutations)

    with open(filename, 'r') as f:
        lines = f.readlines()

    new_content = []
    for i, line in enumerate(lines):
        line = line.strip()

        if line.startswith(">"):  # If it's a header, save it
            new_content.append(line + '\n')
            continue

        match = re.search(pattern, line)
        if match:  # If we found a match
            start = max(0, match.start() - 20)  # 20 positions before, but not negative
            end = min(len(line), match.end() + 20)  # Matched region and 20 positions after

            new_content.append(line[start:end] + '\n')

    # Write the modified content to a new .txt file
    output_filename = filename.rsplit(".", 1)[0] + "_modified.txt"  # Construct the output filename
    print(output_filename)
    with open(output_filename, 'w') as f:
        f.writelines(new_content)


# def add_clustal_header(filename: str) -> None:
#     """Add CLUSTAL header to an alignment file."""
#     with open(filename, 'r') as f:
#         content = f.read()
#     with open(filename, 'w') as f:
#         f.write("CLUSTAL W multiple sequence alignment\n\n" + content)


def find_conserved_regions(alignment, window_size=20, threshold=0.95) -> List[int]:
    """
    Find conserved regions from an alignment.

    Parameters:
    - alignment: A Biopython alignment object.
    - window_size: The size of the region to check for conservation.
    - threshold: Fraction of sequences that need to have the same base to consider it conserved.

    Returns:
    - List of start positions of conserved regions.
    """
    conserved_positions = []
    num_sequences = len(alignment)
    for i in range(len(alignment[0]) - window_size + 1):  # loop over positions in alignment
        window = alignment[:, i:i+window_size]  # extract window from alignment
        if sum(window.count(c) for c in set(window[0]))/num_sequences >= threshold:
            conserved_positions.append(i)
    return conserved_positions


# def fix_sequence_ids(filename: str) -> None:
#     """Replace spaces in sequence IDs with underscores."""
#     with open(filename, 'r') as f:
#         lines = f.readlines()

#     with open(filename, 'w') as f:
#         for line in lines:
#             if line.startswith('>'):
#                 parts = line.split()
#                 fixed_id = "_".join(parts[:2])  # Assuming the ID and the next part need to be merged
#                 rest_of_line = " ".join(parts[2:])
#                 f.write(f"{fixed_id} {rest_of_line}\n")
#             else:
#                 f.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform MSA on FASTA sequences.")

    parser.add_argument("-d", "--directory", required=True, help="Directory to save the FASTA files.")
    parser.add_argument("-s", "--savename", required=True, help="Batch ID.")
    args = parser.parse_args()

    search_aln = f"{args.directory}/{args.savename}*-clustal.aln"

    alignment_files = glob.glob(search_aln)

    # # Correct alignment files by adding headers
    # for alignment_file in alignment_files:
    #     fix_sequence_ids(alignment_file)
    #     add_clustal_header(alignment_file)

    all_conserved_regions = {}

    # Then, process the alignment files as you intended
    for alignment_file in alignment_files:
        extract_consecutive_regions(alignment_file)
    #     alignment = AlignIO.read(alignment_file, "clustal")
    #     conserved_positions = find_conserved_regions(alignment)

    #     # Store conserved regions for each alignment file
    #     all_conserved_regions[alignment_file] = conserved_positions

    #     # Print conserved regions for the current alignment
    #     print(f"Conserved regions for {alignment_file}:")
    #     print(conserved_positions)
    #     print('-' * 80)  # print a separator line
    # print(len(all_conserved_regions))


# def concatenate_alignments(alignment_files):
#     """
#     Concatenate multiple alignment files.

#     Parameters:
#     - alignment_files: List of paths to alignment files.

#     Returns:
#     - A concatenated Biopython alignment object.
#     """
#     # Load the first alignment
#     master_alignment = AlignIO.read(alignment_files[0], "clustal")
#     master_ids = [record.id for record in master_alignment]

#     # Loop through the remaining alignments and concatenate
#     for aln_file in alignment_files[1:]:
#         aln = AlignIO.read(aln_file, "clustal")
#         aln_ids = [record.id for record in aln]

#         # Ensure the sequences are in the same order across all alignments
#         if master_ids != aln_ids:
#             raise ValueError("The sequence IDs do not match across alignment files!")

#         # Concatenate the sequences
#         concatenated_sequences = []
#         for master_record, aln_record in zip(master_alignment, aln):
#             concatenated_seq = master_record.seq + aln_record.seq
#             concatenated_record = SeqRecord(concatenated_seq, id=master_record.id)
#             concatenated_sequences.append(concatenated_record)

#         master_alignment = MultipleSeqAlignment(concatenated_sequences)

#     return master_alignment






#     combined_alignment = concatenate_alignments(alignment_files)

#     # Save the combined alignment to a file
#     output_path = f"{args.directory}/concat_{args.savename}_alignment.aln"
#     print(output_path)
#     AlignIO.write(combined_alignment, output_path, "clustal")
#     print(f"Combined alignment saved to {output_path}")
