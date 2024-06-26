#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from skbio import DistanceMatrix
from skbio.stats.distance import permanova

# Assume 'bc_array' is your precomputed Bray-Curtis dissimilarity matrix
# and 'sample_df' is your metadata DataFrame

# For each categorical variable you are interested in, perform PERMANOVA:

#### EDIT HERE ####
# directory = "/mnt/path/to/Confidence/"
directory = "/mnt/e/SequencingData/CRAAM/analysis/Confidence/"
#### EDIT HERE ####

coords_file = f"{directory}/coord-df.csv"
npfile = f"{directory}/bc-array.npy"
sample_df = pd.read_csv(f"{directory}/sample-df-labels.csv")

with open(npfile, 'rb') as f:
    bc_array = np.load(f)

plotCols = ['BSA', 'Kit', 'PCR', 'Spike', 'TC-Spike', 'Polymerase',
       'Cycles', 'Flow-Cell']

for col in plotCols:
    # Load coordinates and sample labels
    coords_df = pd.read_csv(coords_file)
    coords_df["Sample"] = sample_df[col]  # Add sample names

    # Check if there is more than one unique group in the 'Sample' column
    if len(coords_df["Sample"].unique()) > 1:
        # Calculate PERMANOVA
        dm = DistanceMatrix(bc_array)
        results = permanova(dm, coords_df["Sample"])
        # Print p-value
        print(f"{col} p-value: ", results['p-value'])
    else:
        print(f"No variation OR DATA in {col}, skipping PERMANOVA analysis.")

# for col in plotCols:
#     # Convert coordinates into a pandas DataFrame
#     coords_df = pd.read_csv(coords_file)
#     coords_df["Sample"] = sample_df[col]  # Add sample names

#     # Calculate PERMANOVA
#     dm = DistanceMatrix(bc_array)
#     results = permanova(dm, coords_df["Sample"])

#     # Print p-value
#     print(f"{col} p-value: ", results['p-value'])

# time ./permanova.py
