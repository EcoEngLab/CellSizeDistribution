"""
This script processes amplicon sequencing data and merges it with bacteria cell size data at the genus level.

- Loads all host species TSV files from the input directory
- Pools read counts at the genus level for each host species
- Calculates median cell size values for each genus from the cell size database
- Merges genus-level read counts with genus-level median cell size values
- Outputs merged tables as CSV files for Bacteria only
- Appends summary statistics to each output file
"""

import pandas as pd
import numpy as np
import os
import glob
import re

# --- CONFIGURATION ---
input_dir = "qiime_taxon_tables_split (1)"  # Directory with host species files
cell_size_file = "data\VolumeOutputBacteria.csv"  # Genus-level cell size data
output_dir = "nopool_genus_level_output_bacteria"
os.makedirs(output_dir, exist_ok=True)

# Read cell size data and prepare genus-level medians
cell_size_df = pd.read_csv(cell_size_file)
taxonomy_cols = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']
cell_size_df[taxonomy_cols] = cell_size_df[taxonomy_cols].replace('NA', pd.NA)
cell_size_cols = ['LengthMin', 'LengthMax', 'WidthMin', 'WidthMax', 
                  'VolumeMin', 'VolumeMax', 'AvgVolume', 'GeoMeanVolume', 'SurfaceArea']
for col in cell_size_cols:
    cell_size_df[col] = pd.to_numeric(cell_size_df[col], errors='coerce')
genus_stats = cell_size_df.groupby('Genus').agg({
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
genus_stats = genus_stats.rename(columns={col: f"{col}_genus_median" for col in genus_stats.columns if col != "Genus"})

# Process each host species file
for file_path in glob.glob(os.path.join(input_dir, "*.tsv")):
    df = pd.read_csv(file_path, sep='\t')
    # Ensure 'genus' column exists
    if 'genus' not in df.columns:
        continue

    # Remove contaminants
    contaminant_keywords = ['Mitochondria', 'Chloroplast']
    tax_cols = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for col in tax_cols:
        if col in df.columns:
            df = df[~df[col].astype(str).str.contains('|'.join(contaminant_keywords), case=False, na=False)]

    # Filter for Bacteria only
    if 'domain' in df.columns:
        df = df[df['domain'].astype(str).str.lower() == 'bacteria']

    # Flexible genus matching
    def match_genus(abundance_genus, trait_genera):
        if not isinstance(abundance_genus, str) or not abundance_genus:
            return None
        # Try exact match first
        if abundance_genus in trait_genera:
            return abundance_genus
        for trait_genus in trait_genera:
            if not isinstance(trait_genus, str) or not trait_genus:
                continue
            pattern = r"^" + re.escape(trait_genus) + r"($|_|[(\s\dA-Z])"
            try:
                if re.match(pattern, abundance_genus):
                    return trait_genus
            except re.error:
                continue
        return None

    trait_genera = set(genus_stats['Genus'])
    df['genus_matched'] = df['genus'].apply(lambda g: match_genus(g, trait_genera))

    # Identify sample columns: only those starting with 'SRR', 'ERR', or 'DRR'
    sample_cols = [col for col in df.columns if col.startswith(('SRR', 'ERR', 'DRR'))]

    # Merge with genus-level cell size data
    merged = pd.merge(df, genus_stats, left_on='genus_matched', right_on='Genus', how='left')
    if 'Genus' in merged.columns:
        merged = merged.drop(columns=['Genus'])

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

    # Drop 'genus_matched' and 'genus_matched_rel_abund' columns if they exist
    cols_to_drop = [col for col in merged.columns if col == "genus_matched" or col.startswith("genus_matched_rel_abund")]
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
    group_cols = ['genus']
    cellsize_cols = [col for col in merged.columns if col.endswith('_genus_median') or col == 'Cell shape_genus_median']
    # Only keep one row per genus, summing sample columns and keeping median cell size columns
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
    grouped['LogGeoMeanVolume'] = np.log10(grouped['GeoMeanVolume_genus_median'])

    # Output merged file
    base = os.path.splitext(os.path.basename(file_path))[0]
    output_path = os.path.join(output_dir, f"{base}_genus_with_cellsize.csv")
    grouped.to_csv(output_path, index=False)

    # Calculate summary statistics
    num_total_genera = grouped['genus'].nunique()
    num_matched_genera = grouped[grouped['GeoMeanVolume_genus_median'].notna()]['genus'].nunique()
    num_unmatched_genera = num_total_genera - num_matched_genera
    percent_matched = num_matched_genera / num_total_genera if num_total_genera > 0 else 0

    # For total reads, sum across all sample columns
    total_reads = grouped[sample_cols].sum().sum()
    matched_reads = grouped[grouped['GeoMeanVolume_genus_median'].notna()][sample_cols].sum().sum()
    unmatched_reads = total_reads - matched_reads
    percent_matched_by_reads = matched_reads / total_reads if total_reads > 0 else 0

    # Append summary statistics to the output file
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
