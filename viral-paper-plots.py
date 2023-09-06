#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 12:52:40 2023

@author: mangi
"""

from string import whitespace
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob



def preprocess_data(file: str) -> pd.DataFrame:
    """Load and preprocess the data from a CSV file."""
    df = pd.read_csv(file)
    
    # Remove duplicates
    non_dupe_cols = list(set(df.columns) - set(['spp']))
    df = df.drop_duplicates(subset=non_dupe_cols)
    
    # Extract CFU and species details
    df["CFU"] = df["sample"].str.split("_", expand=True)[3].str.replace("CFU", "", regex=True).astype(float)
    species_split = df["spp"].str.split(" ", expand=True)
    df["species"] = species_split[0] + " " + species_split[1]
    df["genus"] = species_split[0]
    
    return df

def fix_time_values(df: pd.DataFrame) -> pd.DataFrame:
    """Fix any broken time values in the dataset."""
    b = df[df["time"] < 2880]
    a = df[df["time"] > 2880]
    min_a = a.groupby(["sample"])["time"].min().to_dict()

    dfs = [b]
    for k, v in min_a.items():
        aa = a[a["sample"] == k]
        aa['time'] -= v
        dfs.append(aa)

    return pd.concat(dfs)

def filter_species_of_interest(df: pd.DataFrame, species: dict, take_TP: bool) -> pd.DataFrame:
    """Filter the dataset for species of interest."""
    spps_of_int = []

    for key, values in species.items():
        for value in values:
            if take_TP:
                spps = df[(df["spp"].str.contains(key)) & (df["sample"].str.contains(value))]
            else:
                spps = df[(df["spp"].str.contains(key)) & (~df["sample"].str.contains(value))]
            spps_of_int.append(spps)

    return pd.concat(spps_of_int)

def generate_data_for_plot(species: dict, file: str, plot: bool, take_TP: bool, zeroes: bool) -> pd.DataFrame:
    """Process and generate data for plotting."""
    df = preprocess_data(file)
    
    # Fix broken time values
    df = fix_time_values(df)
    
    # TODO: Add plot generation logic if needed

    # Filter for non-zero CFU values
    if zeroes:
        df = df[df["CFU"] > 0]
        df.reset_index(drop=True, inplace=True)
    
    # Filter for species of interest
    df = filter_species_of_interest(df, species, take_TP)
    
    return df

def filter_data(species: dict, files: list, plot: bool, take_TP: bool, zeroes: bool) -> pd.DataFrame:
    """Filter data across multiple files."""
    dfs = [generate_data_for_plot(species, file, plot, take_TP, zeroes) for file in files]
    return pd.concat(dfs).reset_index(drop=True)

# Rest of the code stays as is.



folder = "D:/Dropbox/AA SMART/fourier_data/2023-08-31/"

species = {"Feline leukemia":["FeLV"], "Porcine circovirus": ["PCV"], "Minute virus":["MVM"]}

zeroes = False # include negative controls.
plot = True
experiment_time = 2880

expt = "reps"


sample_rename3 = {
            "20230325_ssDNA_TC-FeLV-PCV1-MVM-10-1-rep1_10CFU_45_36":"10-1",
            "20230325_ssDNA_TC-FeLV-PCV1-MVM-10-1-rep2_10CFU_45_36":"10-1",
            "20230325_ssDNA_TC-FeLV-PCV1-MVM-10-1-rep3_10CFU_45_36":"10-1",
            "20230330_ssDNA_TC-FeLV-PCV1-MVM-1-1-rep1_10CFU_45_36":"1-1",
            "20230330_ssDNA_TC-FeLV-PCV1-MVM-1-1-rep2_10CFU_45_36":"1-1",
            "20230330_ssDNA_TC-FeLV-PCV1-MVM-1-1-rep3_10CFU_45_36":"1-1",
            "20230402_ssDNA_TC-FeLV-PCV1-MVM-0pt1-1-rep1_10CFU_45_36":"0.1-1",
            "20230402_ssDNA_TC-FeLV-PCV1-MVM-0pt1-1-rep2_10CFU_45_36":"0.1-1",
            "20230402_ssDNA_TC-FeLV-PCV1-MVM-0pt1-1-rep3_10CFU_45_36":"0.1-1",
            "20230410_ssDNA_TC-FeLV-PCV1-MVM-0pt01-1-rep1_10CFU_45_36":"0.01-1",
            "20230410_ssDNA_TC-FeLV-PCV1-MVM-0pt01-1-rep2_10CFU_45_36":"0.01-1",
            "20230410_ssDNA_TC-FeLV-PCV1-MVM-0pt01-1-rep3_10CFU_45_36":"0.01-1",
            "20230421_ssDNA_TC-FeLV-PCV1-MVM-0pt01-1-rep1-rep2_10CFU_46_36":"0.01-1",
            "20230421_ssDNA_TC-FeLV-PCV1-MVM-0pt01-1-rep2-rep2_10CFU_46_36":"0.01-1",
            "20230421_ssDNA_TC-FeLV-PCV1-MVM-0pt01-1-rep3-rep2_10CFU_46_36":"0.01-1",
            "20230426_ssDNA_TC-FeLV-PCV1-MVM-0pt001-1-rep1_10CFU_46_36":"0.001-1",
            "20230426_ssDNA_TC-FeLV-PCV1-MVM-0pt001-1-rep2_10CFU_46_36":"0.001-1",
            "20230426_ssDNA_TC-FeLV-PCV1-MVM-0pt001-1-rep3_10CFU_46_36":"0.001-1",
            "20230512_ssDNA_TC-FeLV-PCV1-MVM-0pt0001-1-rep1_10CFU_46_36":"0.0001-1",
            "20230512_ssDNA_TC-FeLV-PCV1-MVM-0pt0001-1-rep2_10CFU_46_36":"0.0001-1",
            "20230512_ssDNA_TC-FeLV-PCV1-MVM-0pt0001-1-rep3_10CFU_46_36":"0.0001-1",
            "20230512_ssDNA_TC-FeLV-PCV1-MVM-0pt0001-1-rep4_10CFU_46_36":"0.0001-1",
            "20230518_ssDNA_TCell-ctrl-1_10CFU_46_36":"T-Cells only",
            "20230518_ssDNA_TCell-ctrl-2_10CFU_46_36":"T-Cells only",
            "20230518_ssDNA_TCell-ctrl-3_10CFU_46_36":"T-Cells only",
            "20230524_ssDNA_Medium-ctrl-1_10CFU_46_36":"Medium control",
            "20230524_ssDNA_Medium-ctrl-2_10CFU_46_36":"Medium control",
            "20230524_ssDNA_Medium-ctrl-3_10CFU_46_36":"Medium control",
            "20230526_ssDNA_Viral-ctrl-1_10CFU_46_36":"Positive control",
            "20230526_ssDNA_Viral-ctrl-2_10CFU_46_36":"Positive control",
            "20230526_ssDNA_Viral-ctrl-3_10CFU_46_36":"Positive control",
            "20230607_ssDNA_Viral-ctrl2-1_10CFU_46_36":"Positive control",
            "20230607_ssDNA_Viral-ctrl2-2_10CFU_46_36":"Positive control",
            "20230607_ssDNA_Viral-ctrl2-3_10CFU_46_36":"Positive control",
            "20230803_ssDNA_TC-FeLV-PCV1-MVM-1-1-rep1_10CFU_47_36":"1-1 Baseline",
            "20230803_ssDNA_TC-FeLV-PCV1-MVM-1-1-rep2_10CFU_47_36":"1-1 Baseline",
            "20230803_ssDNA_TC-FeLV-PCV1-MVM-1-1-rep3_10CFU_47_36":"1-1 Baseline",
            "20230828_ssDNA_TC-FeLV-PCV1-MVM-0pt1-1-rep1_10CFU_47_36":"0.1-1 Baseline",
            "20230828_ssDNA_TC-FeLV-PCV1-MVM-0pt1-1-rep2_10CFU_47_36":"0.1-1 Baseline",
            "20230828_ssDNA_TC-FeLV-PCV1-MVM-0pt1-1-rep3_10CFU_47_36":"0.1-1 Baseline"
        }

files = glob.glob(f"{folder}/new_complete_concatenated_CDFs.csv")

TPs = filter_data(species, files, plot, True, zeroes)
TPs["samp_name"] = TPs["sample"]

TPs = TPs.replace({"samp_name": sample_rename3})

# Dictionary to hold results for different species
species_dfs = {}

# Outer loop over unique 'samp_name'
for samp_name in TPs['sample'].unique():
    
    # Filter the dataframe for the current 'samp_name'
    samp_df = TPs[TPs['sample'] == samp_name]
    
    # Inner loop over unique species within the current 'samp_name'
    for species in samp_df['species'].unique():
        
        # Filter the dataframe for the current species within the current 'samp_name'
        group = samp_df[samp_df['species'] == species]
        
        # Sort by time
        sorted_group = group.sort_values(by='time', ascending=True)
        
        # Compute cumulative frequency
        sorted_group['cumulative_frequency'] = sorted_group['frequency'].cumsum()
        
        # Add the result to the dictionary using the species as the key
        if species not in species_dfs:
            species_dfs[species] = []
        species_dfs[species].append(sorted_group)


# Dictionary to hold the mean results for each species
mean_species_dfs = {}

for species, dataframes in species_dfs.items():
    # Concatenate all dataframes associated with the current species
    concatenated_df = pd.concat(dataframes)
    
    # Group by 'samp_name' and 'time', and compute mean for cdf and cumulative_frequency
    grouped = concatenated_df.groupby(['samp_name', "time"]).agg({
        'cdf': 'mean',
        'cumulative_frequency': 'mean',
        'species': 'first'
    }).reset_index()

    mean_species_dfs[species] = grouped


hue_order = ['10-1', '1-1', '0.1-1', '0.01-1', '0.001-1', '0.0001-1', '1-1 Baseline', '0.1-1 Baseline']


# Loop through each species in the mean_species_dfs dictionary
for species, plot_df in mean_species_dfs.items():
    
    plot_df = plot_df[plot_df['samp_name'].isin(['1-1', '1-1 Baseline'])]
    hue_order = ['1-1','1-1 Baseline',]
    
    f, ax = plt.subplots(figsize=(15, 15))
    
    sns.scatterplot(x='time', y='cumulative_frequency', hue='samp_name', data=plot_df, hue_order=hue_order, s = 100)
    plt.rcParams['legend.title_fontsize'] = 30    
    plt.title(f"Cumulative reads over time for {species}", fontsize=40)
    plt.xlim(0, 360)  # Set x-axis limit to 6h (360 minutes)
    plt.xlabel("Time (in minutes)", fontsize=35)
    plt.ylabel("Cumulative Reads", fontsize=35)
    plt.legend(title="Spike \nconcentration", bbox_to_anchor=(1.01, 0.7), loc=2, borderaxespad=0., prop={"size": 25})

    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25) 

    ax.set_yscale('log')  # Convert y-axis to logarithmic scale
    
    plt.show()


import json

# Convert each DataFrame in the dictionary to a JSON string
json_data = {species: df.to_json(orient="records") for species, df in mean_species_dfs.items()}

# Write the JSON data to a file
with open(f'{folder}/mean_species_dfs.json', 'w') as outfile:
    json.dump(json_data, outfile)


import pandas as pd
import json

# Load the JSON data from the file
with open(f'{folder}mean_species_dfs.json', 'r') as infile:
    json_data = json.load(infile)

# Convert each JSON string back into a DataFrame
loaded_mean_species_dfs = {species: pd.read_json(json_str, orient="records") for species, json_str in json_data.items()}
