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
import re

# --- CONFIGURATION ---
input_dir = "qiime_taxon_tables_split (1)"  # Directory with host species files
cell_size_file = "VolumeOutputBacteria.csv"  # Family-level cell size data
output_dir = "nopool_family_level_output_bacteria"
os.makedirs(output_dir, exist_ok=True)

# Read cell size data and prepare family-level medians
cell_size_df = pd.read_csv(cell_size_file)
taxonomy_cols = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']
cell_size_df[taxonomy_cols] = cell_size_df[taxonomy_cols].replace('NA', pd.NA)
cell_size_cols = ['LengthMin', 'LengthMax', 'WidthMin', 'WidthMax', 
                  'VolumeMin', 'VolumeMax', 'AvgVolume', 'GeoMeanVolume', 'SurfaceArea']
for col in cell_size_cols:
    cell_size_df[col] = pd.to_numeric(cell_size_df[col], errors='coerce')
family_stats = cell_size_df.groupby('Family').agg({
    'LengthMin': 'median',
    'LengthMax': 'median',
    'WidthMin': 'median',
    'WidthMax': 'median',
    'VolumeMin': 'median',
    'VolumeMax': 'median',
    'AvgVolume': 'median',
    'GeoMeanVolume': 'median',
    'SurfaceArea': 'median',
    'Cell shape': lambda x: x.mode()[0] if not x.mode().empty else pd.NA
}).reset_index()
family_stats = family_stats.rename(columns={col: f"{col}_family_median" for col in family_stats.columns if col != "Family"})

# Process each host species file
for file_path in glob.glob(os.path.join(input_dir, "*.tsv")):
    df = pd.read_csv(file_path, sep='\t')
    # Ensure 'family' column exists
    if 'family' not in df.columns:
        continue

    # Remove contaminants
    contaminant_keywords = ['Mitochondria', 'Chloroplast']
    tax_cols = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for col in tax_cols:
        if col in df.columns:
            df = df[~df[col].astype(str).str.contains('|'.join(contaminant_keywords), case=False, na=False)]

    # Flexible family matching
    def match_family(abundance_family, trait_families):
        if not isinstance(abundance_family, str) or not abundance_family:
            return None
        # Try exact match first
        if abundance_family in trait_families:
            return abundance_family
        for trait_family in trait_families:
            if not isinstance(trait_family, str) or not trait_family:
                continue
            # Delimiters: end of string, underscore, parenthesis, space, digit, or capital letter
            pattern = r"^" + re.escape(trait_family) + r"($|_|[(\s\dA-Z])"
            try:
                if re.match(pattern, abundance_family):
                    return trait_family
            except re.error:
                continue  # skip problematic patterns
        return None

    trait_families = set(family_stats['Family'])
    df['family_matched'] = df['family'].apply(lambda fam: match_family(fam, trait_families))

    # Identify sample columns: only those starting with 'SRR', 'ERR', or 'DRR'
    sample_cols = [col for col in df.columns if col.startswith(('SRR', 'ERR', 'DRR'))]

    # Merge with family-level cell size data
    merged = pd.merge(df, family_stats, left_on='family_matched', right_on='Family', how='left')
    if 'Family' in merged.columns:
        merged = merged.drop(columns=['Family'])

    # Convert all sample columns to numeric (coerce errors to NaN)
    for sample in sample_cols:
        merged[sample] = pd.to_numeric(merged[sample], errors='coerce')

    # Calculate relative abundance for each sample column
    for sample in sample_cols:
        total_reads = merged[sample].sum()
        if total_reads > 0:
            merged[f"{sample}_rel_abund"] = merged[sample] / total_reads
        else:
            merged[f"{sample}_rel_abund"] = 0

    # Drop 'family_matched' and 'family_matched_rel_abund' columns if they exist
    cols_to_drop = [col for col in merged.columns if col == "family_matched" or col.startswith("family_matched_rel_abund")]
    merged = merged.drop(columns=cols_to_drop, errors='ignore')

    # Recompute sample_cols to only include columns that still exist in merged and start with SRR/ERR/DRR
    sample_cols = [col for col in merged.columns if col.startswith(('SRR', 'ERR', 'DRR'))]

    # Ensure all sample columns are numeric before summing
    for sample in sample_cols:
        merged[sample] = pd.to_numeric(merged[sample], errors='coerce')

    # Optionally, filter sample_cols to only numeric columns
    sample_cols = [col for col in sample_cols if pd.api.types.is_numeric_dtype(merged[col])]

    # Before grouping, save the list of sample columns (read count columns)
    raw_sample_cols = sample_cols.copy()

    # Identify columns to group by and columns to sum
    group_cols = ['family']  # or ['family_matched'] if you want to use the matched name
    cellsize_cols = [col for col in merged.columns if col.endswith('_family_median') or col == 'Cell shape_family_median']
    # Only keep one row per family, summing sample columns and keeping median cell size columns
    grouped = merged.groupby(group_cols, as_index=False)[sample_cols + cellsize_cols].agg(
        {**{col: 'sum' for col in sample_cols}, **{col: 'first' for col in cellsize_cols}}
    )

    # After grouping, recalculate relative abundance for each sample column (read counts only)
    for sample in raw_sample_cols:
        total_reads = grouped[sample].sum()
        if total_reads > 0:
            grouped[f"{sample}_rel_abund"] = grouped[sample] / total_reads
        else:
            grouped[f"{sample}_rel_abund"] = 0

    # Remove any _rel_abund columns not associated with the original sample columns
    allowed_rel_abund_cols = {f"{sample}_rel_abund" for sample in raw_sample_cols}
    rel_abund_cols_to_drop = [col for col in grouped.columns if col.endswith('_rel_abund') and col not in allowed_rel_abund_cols]
    grouped = grouped.drop(columns=rel_abund_cols_to_drop, errors='ignore')

    # Remove any columns that end with '_rel_abund_rel_abund' before output
    grouped = grouped.drop(columns=[col for col in grouped.columns if col.endswith('_rel_abund_rel_abund')], errors='ignore')

    # Calculate base-10 logarithm of the median GeoMeanVolume for each row
    grouped['LogGeoMeanVolume'] = np.log10(grouped['GeoMeanVolume_family_median'])

    # Now, grouped only has one row per family, with summed abundances and family-level medians
    # You can drop genus/species columns if you want

    # Output merged file
    base = os.path.splitext(os.path.basename(file_path))[0]
    output_path = os.path.join(output_dir, f"{base}_family_with_cellsize.csv")
    grouped.to_csv(output_path, index=False)

    # Calculate summary statistics
    num_total_families = grouped['family'].nunique()
    num_matched_families = grouped[grouped['GeoMeanVolume_family_median'].notna()]['family'].nunique()
    num_unmatched_families = num_total_families - num_matched_families
    percent_matched = num_matched_families / num_total_families if num_total_families > 0 else 0

    # For total reads, sum across all sample columns
    total_reads = grouped[sample_cols].sum().sum()
    matched_reads = grouped[grouped['GeoMeanVolume_family_median'].notna()][sample_cols].sum().sum()
    unmatched_reads = total_reads - matched_reads
    percent_matched_by_reads = matched_reads / total_reads if total_reads > 0 else 0

    # Append summary statistics to the output file
    with open(os.path.join(output_dir, f"{base}_family_with_cellsize.csv"), 'a') as f:
        f.write('\n# Summary Statistics\n')
        f.write(f"# Total unique families: {num_total_families}\n")
        f.write(f"# Matched families: {num_matched_families}\n")
        f.write(f"# Unmatched families: {num_unmatched_families}\n")
        f.write(f"# Percentage matched families: {percent_matched:.2%}\n")
        f.write(f"# Total reads: {total_reads}\n")
        f.write(f"# Reads in matched families: {matched_reads}\n")
        f.write(f"# Reads in unmatched families: {unmatched_reads}\n")
        f.write(f"# Percentage matched by read count: {percent_matched_by_reads:.2%}\n")
