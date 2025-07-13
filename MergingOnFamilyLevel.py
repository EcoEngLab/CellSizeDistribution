"""
This script processes pooled amplicon sequencing data and merges it with bacteria cell size data at the family level.

- Loads all pooled-by-study TSV files from "pooled_output" directory
- Pools read counts at the family level for each study
- Calculates median cell size values for each family from cell size database
- Merges family-level read counts with family-level median cell size values
- Outputs merged tables as TSV files
"""

import pandas as pd
import numpy as np
import os
import glob
from pathlib import Path

# Configuration
pooled_input_dir = "pooled_output"
cell_size_file = "VolumeOutputBacteria.csv"
output_dir = "family_level_output"

# Create output directory if it doesn't exist
Path(output_dir).mkdir(exist_ok=True)

# Load cell size data and calculate median values by family
try:
    bacteria_df = pd.read_csv(cell_size_file)
    
    # Clean up 'NA' values
    taxonomy_cols = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']
    bacteria_df[taxonomy_cols] = bacteria_df[taxonomy_cols].replace('NA', np.nan)
    
    # Ensure cell size columns are numeric
    cell_size_cols = ['LengthMin', 'LengthMax', 'WidthMin', 'WidthMax', 
                     'VolumeMin', 'VolumeMax', 'AvgVolume', 'GeoMeanVolume', 'SurfaceArea']
    for col in cell_size_cols:
        bacteria_df[col] = pd.to_numeric(bacteria_df[col], errors='coerce')
    
    # Calculate median values by family
    family_median_data = {
        'LengthMin': 'median',
        'LengthMax': 'median',
        'WidthMin': 'median',
        'WidthMax': 'median',
        'VolumeMin': 'median',
        'VolumeMax': 'median',
        'AvgVolume': 'median',
        'GeoMeanVolume': 'median',
        'SurfaceArea': 'median',
        'Cell shape': lambda x: x.mode()[0] if not x.mode().empty else np.nan
    }
    
    # Calculate comprehensive family-level statistics
    family_stats = bacteria_df.groupby('Family').agg({
        'LengthMin': ['median', 'std', 'count'],
        'LengthMax': ['median', 'std', 'count'],
        'WidthMin': ['median', 'std', 'count'],
        'WidthMax': ['median', 'std', 'count'],
        'VolumeMin': ['median', 'std', 'count'],
        'VolumeMax': ['median', 'std', 'count'],
        'AvgVolume': ['median', 'std', 'count'],
        'GeoMeanVolume': ['median', 'std', 'count'],
        'SurfaceArea': ['median', 'std', 'count'],
        'Cell shape': ['nunique', 'count', lambda x: x.mode()[0] if not x.mode().empty else np.nan],
        'Species': 'count',
        'Genus': 'nunique'
    }).reset_index()
    
    # Flatten column names
    family_stats.columns = ['Family'] + [f'{col[0]}_{col[1]}' for col in family_stats.columns[1:]]
    
    # Calculate standard errors and confidence intervals
    volume_cols = ['AvgVolume', 'GeoMeanVolume', 'VolumeMin', 'VolumeMax']
    for col in volume_cols:
        if f'{col}_std' in family_stats.columns and f'{col}_count' in family_stats.columns:
            # Standard Error of the Mean
            family_stats[f'{col}_sem'] = family_stats[f'{col}_std'] / np.sqrt(family_stats[f'{col}_count'])
            
            # Coefficient of Variation
            family_stats[f'{col}_cv'] = family_stats[f'{col}_std'] / family_stats[f'{col}_median']
            
            # 95% Confidence Interval (assuming normal distribution)
            family_stats[f'{col}_ci_lower'] = family_stats[f'{col}_median'] - (1.96 * family_stats[f'{col}_sem'])
            family_stats[f'{col}_ci_upper'] = family_stats[f'{col}_median'] + (1.96 * family_stats[f'{col}_sem'])
    
    # Calculate shape diversity metrics
    family_stats['shape_diversity'] = family_stats['Cell shape_nunique'] / family_stats['Cell shape_count']
    family_stats['shape_consistency'] = 1 - family_stats['shape_diversity']
    
    # Coverage statistics
    family_stats['species_coverage'] = family_stats['Species_count']
    family_stats['genus_coverage'] = family_stats['Genus_nunique']
    
    # Rename columns for clarity
    family_stats = family_stats.rename(columns={
        'Cell shape_<lambda_0>': 'Cell shape_dominant',
        'Cell shape_nunique': 'Cell shape_count_unique',
        'Cell shape_count': 'Cell shape_total'
    })
    
    # Create the family summary for merging (keeping original structure for compatibility)
    family_summary_df = family_stats[['Family'] + [col for col in family_stats.columns if 'median' in col]].copy()
    
    # Rename columns to indicate they are family-level medians
    rename_dict = {col: f'{col.replace("_median", "")}_family_median' for col in family_summary_df.columns if col != 'Family'}
    family_summary_df = family_summary_df.rename(columns=rename_dict)
    
except Exception as e:
    print(f"Error loading cell size data: {e}")
    exit()

# Find all pooled TSV files
pattern = os.path.join(pooled_input_dir, "*_pooled_by_study.tsv")
pooled_files = glob.glob(pattern)

if not pooled_files:
    print("No pooled files found")
    exit()

# Process each pooled file
successful = 0
failed = 0

for file_path in pooled_files:
    try:
        # Load the pooled file
        df = pd.read_csv(file_path, sep='\t')
        
        # Filter for Bacteria domain
        df = df[df['domain'].astype(str) == 'Bacteria'].copy()
        
        # Check if required columns exist
        required_cols = ['family', 'study_accession', 'read_count']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            failed += 1
            continue
        
        # Pool data by family and study_accession (keep study_accession separate)
        # Include taxonomic columns in grouping to retain them
        grouping_cols = ['domain', 'phylum', 'class', 'order', 'family', 'study_accession']
        family_pooled_df = df.groupby(grouping_cols, as_index=False)['read_count'].sum()
        
        # Merge with cell size data
        merged_df = pd.merge(
            family_pooled_df,
            family_summary_df,
            left_on='family',
            right_on='Family',
            how='left'
        )
        
        # Add match type column
        merged_df['MatchType'] = 'Family'
        merged_df.loc[merged_df['GeoMeanVolume_family_median'].isna(), 'MatchType'] = 'No Match'
        
        # Calculate relative abundance for each family within each study
        # Calculate total reads per study
        study_totals = merged_df.groupby('study_accession')['read_count'].sum()
        
        # Calculate relative abundance
        merged_df['relative_abundance'] = merged_df.apply(
            lambda row: row['read_count'] / study_totals[row['study_accession']] 
            if study_totals[row['study_accession']] > 0 else 0, 
            axis=1
        )
        
        # Remove the 'Family' column from output
        if 'Family' in merged_df.columns:
            merged_df = merged_df.drop(columns=['Family'])
        
        # Calculate summary statistics
        num_total_families = merged_df['family'].nunique()
        num_matched_families = merged_df[merged_df['GeoMeanVolume_family_median'].notna()]['family'].nunique()
        num_unmatched_families = num_total_families - num_matched_families
        percent_matched = num_matched_families / num_total_families if num_total_families > 0 else 0
        
        # Calculate read count based statistics
        total_reads = merged_df['read_count'].sum()
        matched_reads = merged_df[merged_df['GeoMeanVolume_family_median'].notna()]['read_count'].sum()
        unmatched_reads = total_reads - matched_reads
        percent_matched_by_reads = matched_reads / total_reads if total_reads > 0 else 0
        
        # Create output filename
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_filename = f"{base_name}_family_with_cellsize.csv"
        output_path = os.path.join(output_dir, output_filename)
        
        # Save merged data
        merged_df.to_csv(output_path, index=False)
        
        # Add summary statistics to the file
        with open(output_path, 'a') as f:
            f.write('\n# Summary Statistics\n')
            f.write(f"# Total unique families: {num_total_families}\n")
            f.write(f"# Matched families: {num_matched_families}\n")
            f.write(f"# Unmatched families: {num_unmatched_families}\n")
            f.write(f"# Percentage matched families: {percent_matched:.2%}\n")
            f.write(f"# Total reads: {total_reads}\n")
            f.write(f"# Reads in matched families: {matched_reads}\n")
            f.write(f"# Reads in unmatched families: {unmatched_reads}\n")
            f.write(f"# Percentage matched by read count: {percent_matched_by_reads:.2%}\n")
        
        successful += 1
        
    except Exception as e:
        failed += 1

# Save comprehensive family-level statistics
family_stats_output_path = os.path.join(output_dir, "family_level_statistics.csv")
family_stats.to_csv(family_stats_output_path, index=False)