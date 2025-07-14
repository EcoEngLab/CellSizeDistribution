"""
This script processes pooled amplicon sequencing data and merges it with bacteria and archaea cell size data at the family level.

- Loads all pooled-by-study TSV files from "pooled_output" directory
- Pools read counts at the family level for each study
- Calculates median cell size values for each family from cell size database
- Merges family-level read counts with family-level median cell size values
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
bacteria_output_dir = "family_level_output_bacteria"
archaea_output_dir = "family_level_output_archaea"

# Create output directories if they don't exist
Path(bacteria_output_dir).mkdir(exist_ok=True)
Path(archaea_output_dir).mkdir(exist_ok=True)

def process_domain(cell_size_file, output_dir, domain_name):
    """
    Process data for a specific domain (Bacteria or Archaea)
    """
    print(f"Processing {domain_name} data...")
    
    # Load cell size data and calculate median values by family
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
        family_stats = domain_df.groupby('Family').agg({
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
        
        return family_summary_df, family_stats
        
    except Exception as e:
        print(f"Error loading {domain_name} cell size data: {e}")
        return None, None

# Process both domains
bacteria_family_summary, bacteria_family_stats = process_domain(bacteria_cell_size_file, bacteria_output_dir, "Bacteria")
archaea_family_summary, archaea_family_stats = process_domain(archaea_cell_size_file, archaea_output_dir, "Archaea")

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
        required_cols = ['family', 'study_accession', 'read_count', 'domain']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Missing required columns in {file_path}: {missing_cols}")
            bacteria_failed += 1
            archaea_failed += 1
            continue
        
        # Process Bacteria
        if bacteria_family_summary is not None:
            try:
                bacteria_df = df[df['domain'].astype(str) == 'Bacteria'].copy()
                
                if not bacteria_df.empty:
                    # Pool data by family and study_accession (keep study_accession separate)
                    grouping_cols = ['domain', 'phylum', 'class', 'order', 'family', 'study_accession']
                    bacteria_family_pooled_df = bacteria_df.groupby(grouping_cols, as_index=False)['read_count'].sum()
                    
                    # Merge with cell size data
                    bacteria_merged_df = pd.merge(
                        bacteria_family_pooled_df,
                        bacteria_family_summary,
                        left_on='family',
                        right_on='Family',
                        how='left'
                    )
                    
                    # Add match type column
                    bacteria_merged_df['MatchType'] = 'Family'
                    bacteria_merged_df.loc[bacteria_merged_df['GeoMeanVolume_family_median'].isna(), 'MatchType'] = 'No Match'
                    
                    # Calculate relative abundance for each family within each study
                    study_totals = bacteria_merged_df.groupby('study_accession')['read_count'].sum()
                    bacteria_merged_df['relative_abundance'] = bacteria_merged_df.apply(
                        lambda row: row['read_count'] / study_totals[row['study_accession']] 
                        if study_totals[row['study_accession']] > 0 else 0, 
                        axis=1
                    )
                    
                    # Remove the 'Family' column from output
                    if 'Family' in bacteria_merged_df.columns:
                        bacteria_merged_df = bacteria_merged_df.drop(columns=['Family'])
                    
                    # Create output filename
                    base_name = os.path.splitext(os.path.basename(file_path))[0]
                    bacteria_output_filename = f"{base_name}_bacteria_family_with_cellsize.csv"
                    bacteria_output_path = os.path.join(bacteria_output_dir, bacteria_output_filename)
                    
                    # Save merged data
                    bacteria_merged_df.to_csv(bacteria_output_path, index=False)
                    
                    # Calculate and add summary statistics
                    num_total_families = bacteria_merged_df['family'].nunique()
                    num_matched_families = bacteria_merged_df[bacteria_merged_df['GeoMeanVolume_family_median'].notna()]['family'].nunique()
                    num_unmatched_families = num_total_families - num_matched_families
                    percent_matched = num_matched_families / num_total_families if num_total_families > 0 else 0
                    
                    total_reads = bacteria_merged_df['read_count'].sum()
                    matched_reads = bacteria_merged_df[bacteria_merged_df['GeoMeanVolume_family_median'].notna()]['read_count'].sum()
                    unmatched_reads = total_reads - matched_reads
                    percent_matched_by_reads = matched_reads / total_reads if total_reads > 0 else 0
                    
                    with open(bacteria_output_path, 'a') as f:
                        f.write('\n# Summary Statistics (Bacteria)\n')
                        f.write(f"# Total unique families: {num_total_families}\n")
                        f.write(f"# Matched families: {num_matched_families}\n")
                        f.write(f"# Unmatched families: {num_unmatched_families}\n")
                        f.write(f"# Percentage matched families: {percent_matched:.2%}\n")
                        f.write(f"# Total reads: {total_reads}\n")
                        f.write(f"# Reads in matched families: {matched_reads}\n")
                        f.write(f"# Reads in unmatched families: {unmatched_reads}\n")
                        f.write(f"# Percentage matched by read count: {percent_matched_by_reads:.2%}\n")
                    
                    bacteria_successful += 1
                    
            except Exception as e:
                print(f"Error processing Bacteria for {file_path}: {e}")
                bacteria_failed += 1
        
        # Process Archaea
        if archaea_family_summary is not None:
            try:
                archaea_df = df[df['domain'].astype(str) == 'Archaea'].copy()
                
                if not archaea_df.empty:
                    # Pool data by family and study_accession (keep study_accession separate)
                    grouping_cols = ['domain', 'phylum', 'class', 'order', 'family', 'study_accession']
                    archaea_family_pooled_df = archaea_df.groupby(grouping_cols, as_index=False)['read_count'].sum()
                    
                    # Merge with cell size data
                    archaea_merged_df = pd.merge(
                        archaea_family_pooled_df,
                        archaea_family_summary,
                        left_on='family',
                        right_on='Family',
                        how='left'
                    )
                    
                    # Add match type column
                    archaea_merged_df['MatchType'] = 'Family'
                    archaea_merged_df.loc[archaea_merged_df['GeoMeanVolume_family_median'].isna(), 'MatchType'] = 'No Match'
                    
                    # Calculate relative abundance for each family within each study
                    study_totals = archaea_merged_df.groupby('study_accession')['read_count'].sum()
                    archaea_merged_df['relative_abundance'] = archaea_merged_df.apply(
                        lambda row: row['read_count'] / study_totals[row['study_accession']] 
                        if study_totals[row['study_accession']] > 0 else 0, 
                        axis=1
                    )
                    
                    # Remove the 'Family' column from output
                    if 'Family' in archaea_merged_df.columns:
                        archaea_merged_df = archaea_merged_df.drop(columns=['Family'])
                    
                    # Create output filename
                    base_name = os.path.splitext(os.path.basename(file_path))[0]
                    archaea_output_filename = f"{base_name}_archaea_family_with_cellsize.csv"
                    archaea_output_path = os.path.join(archaea_output_dir, archaea_output_filename)
                    
                    # Save merged data
                    archaea_merged_df.to_csv(archaea_output_path, index=False)
                    
                    # Calculate and add summary statistics
                    num_total_families = archaea_merged_df['family'].nunique()
                    num_matched_families = archaea_merged_df[archaea_merged_df['GeoMeanVolume_family_median'].notna()]['family'].nunique()
                    num_unmatched_families = num_total_families - num_matched_families
                    percent_matched = num_matched_families / num_total_families if num_total_families > 0 else 0
                    
                    total_reads = archaea_merged_df['read_count'].sum()
                    matched_reads = archaea_merged_df[archaea_merged_df['GeoMeanVolume_family_median'].notna()]['read_count'].sum()
                    unmatched_reads = total_reads - matched_reads
                    percent_matched_by_reads = matched_reads / total_reads if total_reads > 0 else 0
                    
                    with open(archaea_output_path, 'a') as f:
                        f.write('\n# Summary Statistics (Archaea)\n')
                        f.write(f"# Total unique families: {num_total_families}\n")
                        f.write(f"# Matched families: {num_matched_families}\n")
                        f.write(f"# Unmatched families: {num_unmatched_families}\n")
                        f.write(f"# Percentage matched families: {percent_matched:.2%}\n")
                        f.write(f"# Total reads: {total_reads}\n")
                        f.write(f"# Reads in matched families: {matched_reads}\n")
                        f.write(f"# Reads in unmatched families: {unmatched_reads}\n")
                        f.write(f"# Percentage matched by read count: {percent_matched_by_reads:.2%}\n")
                    
                    archaea_successful += 1
                    
            except Exception as e:
                print(f"Error processing Archaea for {file_path}: {e}")
                archaea_failed += 1
                
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        bacteria_failed += 1
        archaea_failed += 1

# Save comprehensive family-level statistics for both domains
if bacteria_family_stats is not None:
    bacteria_stats_output_path = os.path.join(bacteria_output_dir, "bacteria_family_level_statistics.csv")
    bacteria_family_stats.to_csv(bacteria_stats_output_path, index=False)

if archaea_family_stats is not None:
    archaea_stats_output_path = os.path.join(archaea_output_dir, "archaea_family_level_statistics.csv")
    archaea_family_stats.to_csv(archaea_stats_output_path, index=False)
