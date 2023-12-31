# Building your custom NCBI database

# These do not work, sadly.
# wget --recursive --no-parent --timestamping --cut-dirs=3 --no-host-directories --accept "*.gbff.gz" -P fungi/ ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/
# wget --recursive --no-parent --timestamping --cut-dirs=3 --no-host-directories --accept "*.gbff.gz" -P protozoa/ ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

# Example run:

cd GitHub/NCBI-download/

# #########################
# download genbank .gbff.gz

time ./download-genbank.py -t "fungi" -o "/mnt/genbank/"
time ./download-genbank.py -t "fungi" -o "/mnt/genbank/" -s # skip existing folders (only useful for initial run that breaks)
time python download-genbank.py -t "fungi" -o "/mnt/genbank/"

time ./download-genbank.py -t "protozoa" -o "/mnt/genbank/"
time python download-genbank.py -t "protozoa" -o "/mnt/genbank/" -f "/mnt/genbank/protozoa-failures.json"
# remember to change target to organisms desired, found here: ftp://ftp.ncbi.nlm.nih.gov/genomes/

# #########################
# convert gbff to fa and extract rDNA

time python3 gbff2fa.py -t "protozoa" -o "/mnt/genbank/"
time python extract_rDNA.py

time python3 gbff2fa.py -t "fungi" -o "/mnt/genbank/"
time python extract_rDNA.py

# #########################
# Meanwhile, start building the seqid2taxid.map file.
# Requires: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz, fastest way is to download it from website.
# this will convert nucl_gb.accession2taxid.gz into a .map, however it is incomplete.
time python NCBI_acc2taxid_seqid2taxid.py

# to make it complete, we need:
# https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
time python taxonIDs.py
# this will use the taxonomy names.dmp and nodes.dmp to fill in the missing taxonIDs.
