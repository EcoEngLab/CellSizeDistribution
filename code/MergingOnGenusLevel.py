"""
This script processes pooled amplicon sequencing data and merges it with bacteria and archaea cell size data at the genus level.

- Loads all pooled-by-study TSV files from "pooled_output" directory
- Pools read counts at the genus level for each study
- Calculates median cell size values for each genus from cell size database
- Merges genus-level read counts with genus-level median cell size values
- Outputs merged tables as TSV files for both Bacteria and Archaea separately
"""

import pandas as pd
import numpy as np
import os
import glob
from pathlib import Path

# Configuration
pooled_input_dir = "pooled_output"
bacteria_cell_size_file = "VolumeOutputBacteria.csv"
archaea_cell_size_file = "VolumeOutputArchaea.csv"
bacteria_output_dir = "genus_level_output_bacteria"
archaea_output_dir = "genus_level_output_archaea"

# Create output directories if they don't exist
Path(bacteria_output_dir).mkdir(exist_ok=True)
Path(archaea_output_dir).mkdir(exist_ok=True)

def process_domain(cell_size_file, output_dir, domain_name):
    """
    Process data for a specific domain (Bacteria or Archaea)
    """
    print(f"Processing {domain_name} data...")
    
    # Load cell size data and calculate median values by genus
    try:
        domain_df = pd.read_csv(cell_size_file)
        
        # Clean up 'NA' values
        taxonomy_cols = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']
        domain_df[taxonomy_cols] = domain_df[taxonomy_cols].replace('NA', np.nan)
        
        # Ensure cell size columns are numeric
        cell_size_cols = ['LengthMin', 'LengthMax', 'WidthMin', 'WidthMax', 
                         'VolumeMin', 'VolumeMax', 'AvgVolume', 'GeoMeanVolume', 'SurfaceArea']
        for col in cell_size_cols:
            domain_df[col] = pd.to_numeric(domain_df[col], errors='coerce')
        
        # Calculate median values by genus
        genus_stats = domain_df.groupby('Genus').agg({
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
            'Family': 'nunique'
        }).reset_index()
        
        # Flatten column names
        genus_stats.columns = ['Genus'] + [f'{col[0]}_{col[1]}' for col in genus_stats.columns[1:]]
        
        # Calculate standard errors and confidence intervals
        volume_cols = ['AvgVolume', 'GeoMeanVolume', 'VolumeMin', 'VolumeMax']
        for col in volume_cols:
            if f'{col}_std' in genus_stats.columns and f'{col}_count' in genus_stats.columns:
                # Standard Error of the Mean
                genus_stats[f'{col}_sem'] = genus_stats[f'{col}_std'] / np.sqrt(genus_stats[f'{col}_count'])
                
                # Coefficient of Variation
                genus_stats[f'{col}_cv'] = genus_stats[f'{col}_std'] / genus_stats[f'{col}_median']
                
                # 95% Confidence Interval (assuming normal distribution)
                genus_stats[f'{col}_ci_lower'] = genus_stats[f'{col}_median'] - (1.96 * genus_stats[f'{col}_sem'])
                genus_stats[f'{col}_ci_upper'] = genus_stats[f'{col}_median'] + (1.96 * genus_stats[f'{col}_sem'])
        
        # Calculate shape diversity metrics
        genus_stats['shape_diversity'] = genus_stats['Cell shape_nunique'] / genus_stats['Cell shape_count']
        genus_stats['shape_consistency'] = 1 - genus_stats['shape_diversity']
        
        # Coverage statistics
        genus_stats['species_coverage'] = genus_stats['Species_count']
        genus_stats['family_coverage'] = genus_stats['Family_nunique']
        
        # Rename columns for clarity
        genus_stats = genus_stats.rename(columns={
            'Cell shape_<lambda_0>': 'Cell shape_dominant',
            'Cell shape_nunique': 'Cell shape_count_unique',
            'Cell shape_count': 'Cell shape_total'
        })
        
        # Create the genus summary for merging (keeping original structure for compatibility)
        genus_summary_df = genus_stats[['Genus'] + [col for col in genus_stats.columns if 'median' in col]].copy()
        
        # Rename columns to indicate they are genus-level medians
        rename_dict = {col: f'{col.replace("_median", "")}_genus_median' for col in genus_summary_df.columns if col != 'Genus'}
        genus_summary_df = genus_summary_df.rename(columns=rename_dict)
        
        return genus_summary_df, genus_stats
        
    except Exception as e:
        print(f"Error loading {domain_name} cell size data: {e}")
        return None, None

# Process both domains
bacteria_genus_summary, bacteria_genus_stats = process_domain(bacteria_cell_size_file, bacteria_output_dir, "Bacteria")
archaea_genus_summary, archaea_genus_stats = process_domain(archaea_cell_size_file, archaea_output_dir, "Archaea")

# Find all pooled TSV files
pattern = os.path.join(pooled_input_dir, "*_pooled_by_study.tsv")
pooled_files = glob.glob(pattern)

if not pooled_files:
    print("No pooled files found")
    exit()

# Process each pooled file for both domains
bacteria_successful = 0
bacteria_failed = 0
archaea_successful = 0
archaea_failed = 0

for file_path in pooled_files:
    try:
        # Load the pooled file
        df = pd.read_csv(file_path, sep='\t')
        
        # Check if required columns exist
        required_cols = ['genus', 'study_accession', 'read_count', 'domain']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Missing required columns in {file_path}: {missing_cols}")
            bacteria_failed += 1
            archaea_failed += 1
            continue
        
        # Process Bacteria
        if bacteria_genus_summary is not None:
            try:
                bacteria_df = df[df['domain'].astype(str) == 'Bacteria'].copy()
                
                if not bacteria_df.empty:
                    # Pool data by genus and study_accession (keep study_accession separate)
                    grouping_cols = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'study_accession']
                    bacteria_genus_pooled_df = bacteria_df.groupby(grouping_cols, as_index=False)['read_count'].sum()
                    
                    # Merge with cell size data
                    bacteria_merged_df = pd.merge(
                        bacteria_genus_pooled_df,
                        bacteria_genus_summary,
                        left_on='genus',
                        right_on='Genus',
                        how='left'
                    )
                    
                    # Add match type column
                    bacteria_merged_df['MatchType'] = 'Genus'
                    bacteria_merged_df.loc[bacteria_merged_df['GeoMeanVolume_genus_median'].isna(), 'MatchType'] = 'No Match'
                    
                    # Calculate relative abundance for each genus within each study
                    study_totals = bacteria_merged_df.groupby('study_accession')['read_count'].sum()
                    bacteria_merged_df['relative_abundance'] = bacteria_merged_df.apply(
                        lambda row: row['read_count'] / study_totals[row['study_accession']] 
                        if study_totals[row['study_accession']] > 0 else 0, 
                        axis=1
                    )
                    
                    # Remove the 'Genus' column from output
                    if 'Genus' in bacteria_merged_df.columns:
                        bacteria_merged_df = bacteria_merged_df.drop(columns=['Genus'])
                    
                    # Create output filename
                    base_name = os.path.splitext(os.path.basename(file_path))[0]
                    bacteria_output_filename = f"{base_name}_bacteria_genus_with_cellsize.csv"
                    bacteria_output_path = os.path.join(bacteria_output_dir, bacteria_output_filename)
                    
                    # Save merged data
                    bacteria_merged_df.to_csv(bacteria_output_path, index=False)
                    
                    # Calculate and add summary statistics
                    num_total_genera = bacteria_merged_df['genus'].nunique()
                    num_matched_genera = bacteria_merged_df[bacteria_merged_df['GeoMeanVolume_genus_median'].notna()]['genus'].nunique()
                    num_unmatched_genera = num_total_genera - num_matched_genera
                    percent_matched = num_matched_genera / num_total_genera if num_total_genera > 0 else 0
                    
                    total_reads = bacteria_merged_df['read_count'].sum()
                    matched_reads = bacteria_merged_df[bacteria_merged_df['GeoMeanVolume_genus_median'].notna()]['read_count'].sum()
                    unmatched_reads = total_reads - matched_reads
                    percent_matched_by_reads = matched_reads / total_reads if total_reads > 0 else 0
                    
                    with open(bacteria_output_path, 'a') as f:
                        f.write('\n# Summary Statistics (Bacteria)\n')
                        f.write(f"# Total unique genera: {num_total_genera}\n")
                        f.write(f"# Matched genera: {num_matched_genera}\n")
                        f.write(f"# Unmatched genera: {num_unmatched_genera}\n")
                        f.write(f"# Percentage matched genera: {percent_matched:.2%}\n")
                        f.write(f"# Total reads: {total_reads}\n")
                        f.write(f"# Reads in matched genera: {matched_reads}\n")
                        f.write(f"# Reads in unmatched genera: {unmatched_reads}\n")
                        f.write(f"# Percentage matched by read count: {percent_matched_by_reads:.2%}\n")
                    
                    bacteria_successful += 1
                    
            except Exception as e:
                print(f"Error processing Bacteria for {file_path}: {e}")
                bacteria_failed += 1
        
        # Process Archaea
        if archaea_genus_summary is not None:
            try:
                archaea_df = df[df['domain'].astype(str) == 'Archaea'].copy()
                
                if not archaea_df.empty:
                    # Pool data by genus and study_accession (keep study_accession separate)
                    grouping_cols = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'study_accession']
                    archaea_genus_pooled_df = archaea_df.groupby(grouping_cols, as_index=False)['read_count'].sum()
                    
                    # Merge with cell size data
                    archaea_merged_df = pd.merge(
                        archaea_genus_pooled_df,
                        archaea_genus_summary,
                        left_on='genus',
                        right_on='Genus',
                        how='left'
                    )
                    
                    # Add match type column
                    archaea_merged_df['MatchType'] = 'Genus'
                    archaea_merged_df.loc[archaea_merged_df['GeoMeanVolume_genus_median'].isna(), 'MatchType'] = 'No Match'
                    
                    # Calculate relative abundance for each genus within each study
                    study_totals = archaea_merged_df.groupby('study_accession')['read_count'].sum()
                    archaea_merged_df['relative_abundance'] = archaea_merged_df.apply(
                        lambda row: row['read_count'] / study_totals[row['study_accession']] 
                        if study_totals[row['study_accession']] > 0 else 0, 
                        axis=1
                    )
                    
                    # Remove the 'Genus' column from output
                    if 'Genus' in archaea_merged_df.columns:
                        archaea_merged_df = archaea_merged_df.drop(columns=['Genus'])
                    
                    # Create output filename
                    base_name = os.path.splitext(os.path.basename(file_path))[0]
                    archaea_output_filename = f"{base_name}_archaea_genus_with_cellsize.csv"
                    archaea_output_path = os.path.join(archaea_output_dir, archaea_output_filename)
                    
                    # Save merged data
                    archaea_merged_df.to_csv(archaea_output_path, index=False)
                    
                    # Calculate and add summary statistics
                    num_total_genera = archaea_merged_df['genus'].nunique()
                    num_matched_genera = archaea_merged_df[archaea_merged_df['GeoMeanVolume_genus_median'].notna()]['genus'].nunique()
                    num_unmatched_genera = num_total_genera - num_matched_genera
                    percent_matched = num_matched_genera / num_total_genera if num_total_genera > 0 else 0
                    
                    total_reads = archaea_merged_df['read_count'].sum()
                    matched_reads = archaea_merged_df[archaea_merged_df['GeoMeanVolume_genus_median'].notna()]['read_count'].sum()
                    unmatched_reads = total_reads - matched_reads
                    percent_matched_by_reads = matched_reads / total_reads if total_reads > 0 else 0
                    
                    with open(archaea_output_path, 'a') as f:
                        f.write('\n# Summary Statistics (Archaea)\n')
                        f.write(f"# Total unique genera: {num_total_genera}\n")
                        f.write(f"# Matched genera: {num_matched_genera}\n")
                        f.write(f"# Unmatched genera: {num_unmatched_genera}\n")
                        f.write(f"# Percentage matched genera: {percent_matched:.2%}\n")
                        f.write(f"# Total reads: {total_reads}\n")
                        f.write(f"# Reads in matched genera: {matched_reads}\n")
                        f.write(f"# Reads in unmatched genera: {unmatched_reads}\n")
                        f.write(f"# Percentage matched by read count: {percent_matched_by_reads:.2%}\n")
                    
                    archaea_successful += 1
                    
            except Exception as e:
                print(f"Error processing Archaea for {file_path}: {e}")
                archaea_failed += 1
                
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        bacteria_failed += 1
        archaea_failed += 1

# Save comprehensive genus-level statistics for both domains
if bacteria_genus_stats is not None:
    bacteria_stats_output_path = os.path.join(bacteria_output_dir, "bacteria_genus_level_statistics.csv")
    bacteria_genus_stats.to_csv(bacteria_stats_output_path, index=False)

if archaea_genus_stats is not None:
    archaea_stats_output_path = os.path.join(archaea_output_dir, "archaea_genus_level_statistics.csv")
    archaea_genus_stats.to_csv(archaea_stats_output_path, index=False)
