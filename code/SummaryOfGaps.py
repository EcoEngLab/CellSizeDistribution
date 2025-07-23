"""
SummaryOfGaps.py

This script processes CSV files from family-level output directories to create a comprehensive
summary of host species, taxonomy, and trait data.

- Reads all CSV files from a specified directory
- Extracts host species names from filenames
- Groups data by taxonomy columns and study_accession
- Calculates average reads across studies for each family
- Combines results into a long-format CSV
"""

import pandas as pd
import os
import glob
import re
from pathlib import Path

# Configuration
input_dir = "family_level_output_archaea"
output_file = "host_family_traits_archaea.csv"

# Check if input directory exists
if not os.path.exists(input_dir):
    exit()

# Get all CSV files in the directory
pattern = os.path.join(input_dir, "*.csv")
csv_files = glob.glob(pattern)

if not csv_files:
    exit()

# Process each file
all_data = []

for file_path in csv_files:
    try:
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Extract host species from filename
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        parts = base_name.split('_')
        
        # Find the pattern where we have Genus_species followed by L and a number
        host_species = base_name
        for i in range(len(parts) - 2):
            if (parts[i+1] and 
                parts[i+2].startswith('L') and 
                parts[i+2][1:].isdigit()):
                genus = parts[i]
                species = parts[i+1]
                host_species = f"{genus} {species}"
                break
        
        # Fallback: try to extract first two parts as genus and species
        if host_species == base_name and len(parts) >= 2:
            host_species = f"{parts[0]} {parts[1]}"
        
        # Identify taxonomy columns
        taxonomy_cols = ['domain', 'phylum', 'class', 'order', 'family']
        available_taxonomy_cols = [col for col in taxonomy_cols if col in df.columns]
        
        # Identify median trait columns
        median_trait_cols = [col for col in df.columns if '_median' in col]
        
        # Check if required columns exist
        required_cols = ['study_accession', 'read_count'] + available_taxonomy_cols
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            continue
        
        # Calculate average reads for each family across studies
        taxonomy_grouped = df.groupby(available_taxonomy_cols, as_index=False).agg({
            'read_count': 'mean',  # Average across studies
            **{col: 'first' for col in median_trait_cols}  # Keep first value for traits
        }).rename(columns={'read_count': 'average_reads'})
        
        # Add host species column
        taxonomy_grouped['host_species'] = host_species
        
        # Reorder columns to match desired output format
        output_cols = ['host_species'] + available_taxonomy_cols + ['average_reads'] + median_trait_cols
        taxonomy_grouped = taxonomy_grouped[output_cols]
        
        all_data.append(taxonomy_grouped)
        
    except Exception as e:
        continue

if not all_data:
    exit()

# Combine all processed data
combined_df = pd.concat(all_data, ignore_index=True)

# Sort by host_species and taxonomy
sort_cols = ['host_species', 'domain', 'phylum', 'class', 'order', 'family']
available_sort_cols = [col for col in sort_cols if col in combined_df.columns]
combined_df = combined_df.sort_values(available_sort_cols)

# Save the combined data
combined_df.to_csv(output_file, index=False)
