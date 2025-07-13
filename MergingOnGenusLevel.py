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
output_dir = "genus_level_output"  # Change output dir for clarity

Path(output_dir).mkdir(exist_ok=True)

# Load cell size data and calculate median values by genus
try:
    bacteria_df = pd.read_csv(cell_size_file)
    taxonomy_cols = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']
    bacteria_df[taxonomy_cols] = bacteria_df[taxonomy_cols].replace('NA', np.nan)
    cell_size_cols = ['LengthMin', 'LengthMax', 'WidthMin', 'WidthMax', 
                     'VolumeMin', 'VolumeMax', 'AvgVolume', 'GeoMeanVolume', 'SurfaceArea']
    for col in cell_size_cols:
        bacteria_df[col] = pd.to_numeric(bacteria_df[col], errors='coerce')
    # Calculate median values by genus
    genus_stats = bacteria_df.groupby('Genus').agg({
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
    genus_stats.columns = ['Genus'] + [f'{col[0]}_{col[1]}' for col in genus_stats.columns[1:]]
    # Standard errors, CIs, etc. (same as before, just for genus)
    volume_cols = ['AvgVolume', 'GeoMeanVolume', 'VolumeMin', 'VolumeMax']
    for col in volume_cols:
        if f'{col}_std' in genus_stats.columns and f'{col}_count' in genus_stats.columns:
            genus_stats[f'{col}_sem'] = genus_stats[f'{col}_std'] / np.sqrt(genus_stats[f'{col}_count'])
            genus_stats[f'{col}_cv'] = genus_stats[f'{col}_std'] / genus_stats[f'{col}_median']
            genus_stats[f'{col}_ci_lower'] = genus_stats[f'{col}_median'] - (1.96 * genus_stats[f'{col}_sem'])
            genus_stats[f'{col}_ci_upper'] = genus_stats[f'{col}_median'] + (1.96 * genus_stats[f'{col}_sem'])
    genus_stats['shape_diversity'] = genus_stats['Cell shape_nunique'] / genus_stats['Cell shape_count']
    genus_stats['shape_consistency'] = 1 - genus_stats['shape_diversity']
    genus_stats['species_coverage'] = genus_stats['Species_count']
    genus_stats['family_coverage'] = genus_stats['Family_nunique']
    genus_stats = genus_stats.rename(columns={
        'Cell shape_<lambda_0>': 'Cell shape_dominant',
        'Cell shape_nunique': 'Cell shape_count_unique',
        'Cell shape_count': 'Cell shape_total'
    })
    genus_summary_df = genus_stats[['Genus'] + [col for col in genus_stats.columns if 'median' in col]].copy()
    rename_dict = {col: f'{col.replace("_median", "")}_genus_median' for col in genus_summary_df.columns if col != 'Genus'}
    genus_summary_df = genus_summary_df.rename(columns=rename_dict)
except Exception as e:
    print(f"Error loading cell size data: {e}")
    exit()

# Find all pooled TSV files
pattern = os.path.join(pooled_input_dir, "*_pooled_by_study.tsv")
pooled_files = glob.glob(pattern)
if not pooled_files:
    print("No pooled files found")
    exit()

successful = 0
failed = 0

for file_path in pooled_files:
    try:
        df = pd.read_csv(file_path, sep='\t')
        df = df[df['domain'].astype(str) == 'Bacteria'].copy()
        required_cols = ['genus', 'study_accession', 'read_count']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            failed += 1
            continue
        grouping_cols = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'study_accession']
        genus_pooled_df = df.groupby(grouping_cols, as_index=False)['read_count'].sum()
        merged_df = pd.merge(
            genus_pooled_df,
            genus_summary_df,
            left_on='genus',
            right_on='Genus',
            how='left'
        )
        merged_df['MatchType'] = 'Genus'
        merged_df.loc[merged_df['GeoMeanVolume_genus_median'].isna(), 'MatchType'] = 'No Match'
        study_totals = merged_df.groupby('study_accession')['read_count'].sum()
        merged_df['relative_abundance'] = merged_df.apply(
            lambda row: row['read_count'] / study_totals[row['study_accession']] 
            if study_totals[row['study_accession']] > 0 else 0, 
            axis=1
        )
        if 'Genus' in merged_df.columns:
            merged_df = merged_df.drop(columns=['Genus'])
        num_total_genera = merged_df['genus'].nunique()
        num_matched_genera = merged_df[merged_df['GeoMeanVolume_genus_median'].notna()]['genus'].nunique()
        num_unmatched_genera = num_total_genera - num_matched_genera
        percent_matched = num_matched_genera / num_total_genera if num_total_genera > 0 else 0
        total_reads = merged_df['read_count'].sum()
        matched_reads = merged_df[merged_df['GeoMeanVolume_genus_median'].notna()]['read_count'].sum()
        unmatched_reads = total_reads - matched_reads
        percent_matched_by_reads = matched_reads / total_reads if total_reads > 0 else 0
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_filename = f"{base_name}_genus_with_cellsize.csv"
        output_path = os.path.join(output_dir, output_filename)
        merged_df.to_csv(output_path, index=False)
        with open(output_path, 'a') as f:
            f.write('\n# Summary Statistics\n')
            f.write(f"# Total unique genera: {num_total_genera}\n")
            f.write(f"# Matched genera: {num_matched_genera}\n")
            f.write(f"# Unmatched genera: {num_unmatched_genera}\n")
            f.write(f"# Percentage matched genera: {percent_matched:.2%}\n")
            f.write(f"# Total reads: {total_reads}\n")
            f.write(f"# Reads in matched genera: {matched_reads}\n")
            f.write(f"# Reads in unmatched genera: {unmatched_reads}\n")
            f.write(f"# Percentage matched by read count: {percent_matched_by_reads:.2%}\n")
        successful += 1
    except Exception as e:
        failed += 1

# Save comprehensive genus-level statistics
genus_stats_output_path = os.path.join(output_dir, "genus_level_statistics.csv")
genus_stats.to_csv(genus_stats_output_path, index=False)