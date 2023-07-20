# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 13:29:27 2023

@author: mangi
"""

#### EDIT HERE ####
output_dir = "/mnt/genbank/"
#### EDIT HERE ####

print("gunzip -r nucl_gb.accession2taxid.gz")

with open(f'{output_dir}nucl_gb.accession2taxid', 'r') as input_file, open(f'{output_dir}seqid2taxid.map', 'w') as output_file:
    print(f'{output_dir}seqid2taxid.map')
    for line in input_file:
        accession, accession_version, taxid, gi = line.strip().split('\t')
        output_file.write(f"{accession_version}\t{taxid}\n")
