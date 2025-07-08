import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set the folder containing the TSV files
folder_path = 'qiime_taxon_tables_split (1)'  

# Find all .tsv files in the folder
tsv_files = glob.glob(os.path.join(folder_path, '*.tsv'))
print(f"Found {len(tsv_files)} TSV files.")

# Lists to collect all sample-level and species-level proportions
all_sample_proportions = []
all_species_proportions = []

# Terms considered as ambiguous/unclassified for genus
ambiguous_terms = {'uncultured', 'gut_microbiome', 'nan', 'Chloroplast'}

def is_unclassified_below_family(genus):
    # Check for empty or null
    if pd.isnull(genus) or str(genus).strip() == '':
        return True
    genus_str = str(genus).strip().lower()
    if genus_str in ambiguous_terms:
        return True
    # Check for digits, underscore, or hyphen in genus
    if any(char.isdigit() for char in genus_str) or '_' in genus_str or '-' in genus_str:
        return True
    return False

for idx, tsv_file in enumerate(tsv_files, 1):
    print(f"Processing file {idx}/{len(tsv_files)}: {os.path.basename(tsv_file)}")
    try:
        df = pd.read_csv(tsv_file, sep='\t')
    except Exception as e:
        print(f"  Could not read {tsv_file}: {e}")
        continue

    # Filter to only include rows where domain is 'Bacteria'
    if 'domain' in df.columns:
        df = df[df['domain'].astype(str).str.lower() == 'bacteria']
        if df.empty:
            print(f"  No 'Bacteria' domain rows in {tsv_file}. Skipping.")
            continue

    # Identify sample columns (start with 'SRR')
    sample_cols = [col for col in df.columns if col.startswith('SRR')]
    if not sample_cols:
        print(f"  No sample columns found in {tsv_file}. Skipping.")
        continue

    # Classify unclassified-below-family ASVs (only using genus)
    mask_family_filled = df['family'].notnull() & (df['family'].astype(str).str.strip() != '')
    mask_genus_unclassified = df['genus'].apply(is_unclassified_below_family)
    mask_unclassified = mask_family_filled & mask_genus_unclassified

    # For each sample, calculate sample-level proportion hits
    for sample_col in sample_cols:
        total_reads = df[sample_col].sum()
        unclassified_reads = df.loc[mask_unclassified, sample_col].sum()
        if total_reads > 0:
            proportion = unclassified_reads / total_reads
            all_sample_proportions.append(proportion)
        else:
            all_sample_proportions.append(np.nan)

    # For the host species as a whole, calculate species-level proportion hits
    total_reads_all = df[sample_cols].values.sum()
    unclassified_reads_all = df.loc[mask_unclassified, sample_cols].values.sum()
    if total_reads_all > 0:
        species_proportion = unclassified_reads_all / total_reads_all
        all_species_proportions.append(species_proportion)
    else:
        all_species_proportions.append(np.nan)

print("\nFinished processing all files.")

# Remove NaNs from the lists
all_sample_proportions = [x for x in all_sample_proportions if pd.notnull(x)]
all_species_proportions = [x for x in all_species_proportions if pd.notnull(x)]

# Plot histogram for sample-level proportions
plt.figure(figsize=(10, 6))
plt.hist(all_sample_proportions, bins=30, color='skyblue', edgecolor='black')
plt.xlabel('Sample-level Proportion Hits (Unclassified Below Family)')
plt.ylabel('Number of Samples')
plt.title('Distribution of Sample-level Proportion Hits Across All Samples')
plt.tight_layout()
plt.show()

# Plot histogram for species-level proportions
plt.figure(figsize=(10, 6))
plt.hist(all_species_proportions, bins=30, color='salmon', edgecolor='black')
plt.xlabel('Species-level Proportion Hits (Unclassified Below Family)')
plt.ylabel('Number of Host Species')
plt.title('Distribution of Species-level Proportion Hits Across All Host Species')
plt.tight_layout()
plt.show()
