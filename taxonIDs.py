def create_taxid_to_index_dict(nodes_file):
    taxid_to_index = {}
    with open(nodes_file, 'r') as f:
        for line in f:
            fields = line.split('\t|\t')
            taxid = fields[0]
            taxid_to_index[taxid] = taxid  # set the taxid as the index
    return taxid_to_index


def create_index_to_name_dict(names_file):
    index_to_name = {}
    with open(names_file, 'r') as f:
        for line in f:
            fields = line.split('\t|\t')
            index = fields[0]
            name = fields[1]
            name_class = fields[3].replace("|", "").replace(" ", "").strip()
            if name_class.strip() == "scientificname":
                index_to_name[index] = name
    return index_to_name

#### EDIT HERE ####
taxonomy_dir = "/mnt/e/SequencingData/genbank/taxonomy/"
output_dir = "/mnt/e/SequencingData/genbank/"
#### EDIT HERE ####

nodes_file = f"{taxonomy_dir}nodes.dmp"  # replace with your actual path
names_file = f"{taxonomy_dir}names.dmp"  # replace with your actual path

taxid_to_index = create_taxid_to_index_dict(nodes_file)
index_to_name = create_index_to_name_dict(names_file)

print(taxid_to_index['65357']) # Albugo candida
print(index_to_name['65357']) # Albugo candida

print(len(taxid_to_index), len(index_to_name))

print("gunzip -r nucl_gb.accession2taxid.gz")

with open(f'{output_dir}nucl_gb.accession2taxid', 'r') as input_file, open(f'{output_dir}seqid2taxid-names.map', 'w') as output_file:
    print(f'{output_dir}seqid2taxid-names.map')
    for line in input_file:
        accession, accession_version, taxid, gi = line.strip().split('\t')
        index = taxid_to_index.get(taxid, None)
        if index:
            tax_name = index_to_name.get(index, 'Unknown')
            output_file.write(f"{accession_version}\t{taxid}\t{tax_name}\n")
        else:
            print(f"Could not retrieve index for taxid: {taxid}")
