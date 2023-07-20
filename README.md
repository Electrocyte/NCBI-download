# Building your custom NCBI database

## Prerequisites

- Python 3.10 
- Pipenv (Install it using `pip install pipenv`)

## Prepare the Pipenv environment to run the scripts

1. Clone the repository:

```bash
git clone https://github.com/Electrocyte/NCBI-download.git
```

2. Navigate to the project directory:

```bash
cd NCBI-download/
```

3. Install dependencies using Pipenv:

```bash
pipenv install
```

This will get you a starting point for taxonomy:
```bash
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
```

### Example run:


## Running the Project

### Activate the Pipenv shell:

```bash
cd NCBI-download/
pipenv shell
```

# #########################
### Download genbank .gbff.gz

```bash
time ./download-genbank.py -t "fungi" -o "/mnt/genbank/"
time ./download-genbank.py -t "fungi" -o "/mnt/genbank/" -s # skip existing folders (only useful for initial run that breaks)
time python download-genbank.py -t "fungi" -o "/mnt/genbank/"

time ./download-genbank.py -t "protozoa" -o "/mnt/genbank/"
time python download-genbank.py -t "protozoa" -o "/mnt/genbank/" -f "/mnt/genbank/protozoa-failures.json"
```

### Remember to change target to organisms desired, found here: ftp://ftp.ncbi.nlm.nih.gov/genomes/


# #########################
### Convert gbff to fa and extract rDNA

```bash
time python3 gbff2fa.py -t "protozoa" -o "/mnt/genbank/"
time python extract_rDNA.py
```
extract_rDNA.py will require you to manually set the parameters within the script.

```bash
time python3 gbff2fa.py -t "fungi" -o "/mnt/genbank/"
time python extract_rDNA.py
```


# #########################
### Meanwhile, start building the seqid2taxid.map file.
Requires: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz, fastest way is to download it from website.
This will convert of nucl_gb.accession2taxid.gz into a .map, however it is incomplete.
```bash
time python NCBI_acc2taxid_seqid2taxid.py
```


### To make it complete, we need:

https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
```bash
time python taxonIDs.py
```
This will use the taxonomy names.dmp and nodes.dmp to fill in the missing taxonIDs.


### These do not work, sadly.

```bash
wget --recursive --no-parent --timestamping --cut-dirs=3 --no-host-directories --accept "*.gbff.gz" -P fungi/ ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/
wget --recursive --no-parent --timestamping --cut-dirs=3 --no-host-directories --accept "*.gbff.gz" -P protozoa/ ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/
```
