import re
import argparse
from typing import List, Dict

# Step 1: Parse fasta to get organism names
def parse_fasta(filename: str) -> List[str]:
    """
    Parses a FASTA file to extract organism names from headers.

    :param filename: Path to the FASTA file.
    :return: A list of organism names from the FASTA headers.
    """
    with open(filename, 'r') as f:
        organisms = [line.strip().split('|')[4] for line in f if line.startswith(">")]
    return organisms


# Step 2: Fetch taxonomy IDs from names.dmp
def get_taxonomy_id_from_names(seq_names: List[str], names_file: str = "names.dmp") -> Dict[str, str]:
    """
    Fetches taxonomy IDs corresponding to organism names from the names.dmp file.

    :param seq_names: List of organism names.
    :param names_file: Path to the names.dmp file.
    :return: Dictionary with organism names as keys and their taxonomy IDs as values.
    """
    tax_ids = {}
    with open(names_file, 'r') as f:
        for line in f:
            for name in seq_names:
                if re.search(f"\t\|\t{name}\t\|\t", line) and "scientific name" in line:
                    tax_id = line.split("\t|\t")[0]
                    tax_ids[name] = tax_id
    return tax_ids


# Step 3: Extract relevant lines from nodes.dmp based on taxonomy IDs
def extract_from_nodes(tax_ids: Dict[str, str], nodes_file: str = "nodes.dmp") -> List[str]:
    """
    Extracts lines from the nodes.dmp file corresponding to given taxonomy IDs.

    :param tax_ids: Dictionary of organism names and their taxonomy IDs.
    :param nodes_file: Path to the nodes.dmp file.
    :return: A list of relevant lines from the nodes.dmp file.
    """
    relevant_nodes = []
    with open(nodes_file, 'r') as f:
        for line in f:
            for _, tax_id in tax_ids.items():
                if line.startswith(f"{tax_id}\t|"):
                    relevant_nodes.append(line)
    return relevant_nodes


def main(args):
    organism_names = parse_fasta(args.fasta_file)
    tax_ids = get_taxonomy_id_from_names(organism_names, names_file=args.names_file)
    relevant_nodes = extract_from_nodes(tax_ids, nodes_file=args.nodes_file)

    print("Relevant Taxonomy IDs:", tax_ids)
    print("\nRelevant Nodes from nodes.dmp:")
    for line in relevant_nodes:
        print(line, end='')

    # Saving relevant data to new files at the specified save location
    with open(f'{args.save_location}/rvdb_names.dmp', 'w') as names_out:
        for name, tid in tax_ids.items():
            names_out.write(f"{tid}\t|\t{name}\t|\t\n")

    with open(f'{args.save_location}/rvdb_nodes.dmp', 'w') as nodes_out:
        nodes_out.writelines(relevant_nodes)


# Main function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract taxonomy information from given files.")

    parser.add_argument("-f", "--fasta_file", type=str, help="Path to the FASTA file.")
    parser.add_argument("-n", "--names_file", type=str, help="Path to the names.dmp file.")
    parser.add_argument("-o", "--nodes_file", type=str, help="Path to the nodes.dmp file.")
    parser.add_argument("-s", "--save_location", type=str, default=".", help="Directory to save the new .dmp files. Default is the current directory.")

    args = parser.parse_args()

    main(args)
