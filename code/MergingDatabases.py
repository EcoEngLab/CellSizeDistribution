"""
This script processes amplicon sequencing data and merges it with a bacteria cell size database.

- Loads an ASV (amplicon sequence variant) table and a bacteria cell size database.
- Filters ASVs to only those in the Bacteria domain and with valid taxonomy.
- Calculates relative abundances for each sample.
- Merges ASV data with cell size data at the species level (matching on genus and species), and falls back to genus-level averages if no species match is found.
- Outputs a merged table and summary statistics about the matching process.
"""

import pandas as pd
import numpy as np

# Load ASV data
try:
    asv_df = pd.read_csv('Aepyceros_melampus_L7.tsv', sep='\t')
    print(f"Initial ASV data shape: {asv_df.shape}")
except FileNotFoundError:
    raise FileNotFoundError("'Aepyceros_melampus_L7.tsv' not found. Please ensure the file is in the correct directory.")

# Load bacteria cell size data
try:
    bacteria_df = pd.read_csv('OutputBacteria.csv')
    print(f"Initial Bacteria data shape: {bacteria_df.shape}")
except FileNotFoundError:
    raise FileNotFoundError("Required file 'OutputBacteria.csv' is missing.")

# Filter for Bacteria domain
asv_df_filtered = asv_df[asv_df['domain'].astype(str) == 'Bacteria'].copy()
print(f"ASV data shape after filtering domain to 'Bacteria': {asv_df_filtered.shape}")

# Remove rows with all NA across taxonomy levels
taxonomy_cols = ['phylum', 'class', 'order', 'family', 'genus', 'species']
asv_df_filtered[taxonomy_cols] = asv_df_filtered[taxonomy_cols].replace('NA', np.nan)
asv_df_filtered = asv_df_filtered.dropna(subset=taxonomy_cols, how='all')
print(f"ASV data shape after removing rows with all NA taxonomy: {asv_df_filtered.shape}")

# Identify sample columns
sample_cols = [col for col in asv_df_filtered.columns if col.startswith('SRR')]
if not sample_cols:
    print("Warning: No sample columns (starting with 'SRR') found in ASV data. Cannot calculate relative abundance.")
else:
    # Calculate total counts per sample and relative abundance
    for col in sample_cols:
        total_counts = asv_df_filtered[col].sum()
        asv_df_filtered[f'Relative_{col}'] = asv_df_filtered[col] / total_counts if total_counts > 0 else 0
    print("ASV DataFrame with calculated Relative Abundances:")
    print(asv_df_filtered.head())

# Clean up 'NA' values in bacteria df
bacteria_tax_cols = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']
bacteria_df[bacteria_tax_cols] = bacteria_df[bacteria_tax_cols].replace('NA', np.nan)
for col in ['AvgVolume', 'GeoMeanVolume', 'SurfaceArea']:
    bacteria_df[col] = pd.to_numeric(bacteria_df[col], errors='coerce')

# clean up 'NA' values in bacteria df
taxonomy_cols_bacteria = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']
bacteria_df[taxonomy_cols_bacteria] = bacteria_df[taxonomy_cols_bacteria].replace('NA', np.nan)
# ensure cell volume is numeric
bacteria_df['AvgVolume'] = pd.to_numeric(bacteria_df['AvgVolume'], errors='coerce')
bacteria_df['GeoMeanVolume'] = pd.to_numeric(bacteria_df['GeoMeanVolume'], errors='coerce')
bacteria_df['SurfaceArea'] = pd.to_numeric(bacteria_df['SurfaceArea'], errors='coerce')

# create a genus-level summary
genus_summary_data = {
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
genus_summary_df = bacteria_df.groupby('Genus').agg(genus_summary_data).reset_index()
genus_summary_df = genus_summary_df.rename(columns={col: f'{col}_genus_avg' for col in genus_summary_data.keys()})

# Prepare species column for matching (replace underscores with spaces)
asv_df_filtered['species_for_merge'] = asv_df_filtered['species'].astype(str).str.replace('_', ' ')
bacteria_df['Species'] = bacteria_df['Species'].astype(str)

# Merge on species and genus
merged_df = pd.merge(
    asv_df_filtered, bacteria_df,
    left_on=['species_for_merge', 'genus'],
    right_on=['Species', 'Genus'],
    how='left', suffixes=('_amplicon', '_cell_data')
)
merged_df['MatchType'] = None
merged_df.loc[merged_df['GeoMeanVolume'].notna(), 'MatchType'] = 'Species'
print(f"Merged DataFrame shape (species-level): {merged_df.shape}")

# Genus-level fallback for unmatched
unmatched = merged_df[merged_df['GeoMeanVolume'].isnull()].copy()
if not unmatched.empty:
    unmatched_merged = pd.merge(unmatched, genus_summary_df, left_on='genus', right_on='Genus', how='left')
    final_merged_df = merged_df.copy()
    for idx, row in unmatched_merged.iterrows():
        target_row = row.name
        for col_prefix in genus_summary_data.keys():
            orig_col = f"{col_prefix}"
            genus_col = f"{col_prefix}_genus_avg"
            if orig_col in final_merged_df.columns and pd.isnull(final_merged_df.loc[target_row, orig_col]):
                if genus_col in row.index:
                    final_merged_df.loc[target_row, orig_col] = row[genus_col]
                    if pd.notna(row[genus_col]):
                        final_merged_df.loc[target_row, 'MatchType'] = 'Genus'
    print(f"Final merged DataFrame shape (with genus fallback): {final_merged_df.shape}")
else:
    final_merged_df = merged_df.copy()
    print("All entries matched at the species level; no genus-level fallback needed.")

# Remove unwanted columns before saving
cols_to_drop = [
    'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain', 'species_for_merge'
]
final_merged_df = final_merged_df.drop(columns=[col for col in cols_to_drop if col in final_merged_df.columns])

# Calculate summary statistics
num_total = final_merged_df.shape[0]
num_matched = final_merged_df['GeoMeanVolume'].notnull().sum()
num_unmatched = final_merged_df['GeoMeanVolume'].isnull().sum()
matched_df = final_merged_df[final_merged_df['GeoMeanVolume'].notnull()]
num_species = (matched_df['MatchType'] == 'Species').sum()
num_genus = (matched_df['MatchType'] == 'Genus').sum()
prop_species = num_species / num_matched if num_matched > 0 else 0
prop_genus = num_genus / num_matched if num_matched > 0 else 0
percent_matched = num_matched / num_total if num_total > 0 else 0

# Save main data and append summary statistics
output_file_name = 'Aepyceros_melampus_community_cellsize_data.csv'
try:
    final_merged_df.to_csv(output_file_name, index=False)
    with open(output_file_name, 'a') as f:
        f.write('\n# Summary Statistics\n')
        f.write(f"# Total ASVs: {num_total}\n")
        f.write(f"# Matched ASVs: {num_matched}\n")
        f.write(f"# Unmatched ASVs: {num_unmatched}\n")
        f.write(f"# Percentage matched ASVs: {percent_matched:.2%}\n")
        f.write(f"# Proportion matched at species level: {prop_species:.2%}\n")
        f.write(f"# Proportion matched at genus level: {prop_genus:.2%}\n")
    print(f"Successfully saved results and summary to {output_file_name}")
except Exception as e:
    print(f"Error saving DataFrame to CSV: {e}")