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

# Lists for family-level classified/unclassified proportions
all_sample_family_classified_props = []
all_sample_family_unclassified_props = []
# Add lists for species-level (host-level) family classified/unclassified
all_species_family_classified_props = []
all_species_family_unclassified_props = []

# Terms considered as ambiguous/unclassified for genus/family
ambiguous_terms = {'uncultured', 'gut_microbiome', 'nan', 'Chloroplast'}

def is_unclassified_below_family(taxon):
    # Check for empty or null
    if pd.isnull(taxon) or str(taxon).strip() == '':
        return True
    taxon_str = str(taxon).strip().lower()
    if taxon_str in ambiguous_terms:
        return True
    # Check for digits, underscore, or hyphen in taxon
    if any(char.isdigit() for char in taxon_str) or '_' in taxon_str or '-' in taxon_str:
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

    # Classify unclassified-below-family ASVs (using genus for original, family for new)
    mask_family_filled = df['family'].notnull() & (df['family'].astype(str).str.strip() != '')
    mask_genus_unclassified = df['genus'].apply(is_unclassified_below_family)
    mask_unclassified = mask_family_filled & mask_genus_unclassified

    # For each sample, calculate sample-level proportion hits (original, genus-based)
    for sample_col in sample_cols:
        total_reads = df[sample_col].sum()
        unclassified_reads = df.loc[mask_unclassified, sample_col].sum()
        if total_reads > 0:
            proportion = unclassified_reads / total_reads
            all_sample_proportions.append(proportion)
        else:
            all_sample_proportions.append(np.nan)

    # For the host species as a whole, calculate species-level proportion hits (original, genus-based)
    total_reads_all = df[sample_cols].values.sum()
    unclassified_reads_all = df.loc[mask_unclassified, sample_cols].values.sum()
    if total_reads_all > 0:
        species_proportion = unclassified_reads_all / total_reads_all
        all_species_proportions.append(species_proportion)
    else:
        all_species_proportions.append(np.nan)

    # --- Family-level classified/unclassified proportions (new) ---
    mask_family_unclassified = df['family'].apply(is_unclassified_below_family)
    mask_family_classified = ~mask_family_unclassified
    for sample_col in sample_cols:
        total_reads = df[sample_col].sum()
        unclassified_reads = df.loc[mask_family_unclassified, sample_col].sum()
        classified_reads = df.loc[mask_family_classified, sample_col].sum()
        if total_reads > 0:
            all_sample_family_unclassified_props.append(unclassified_reads / total_reads)
            all_sample_family_classified_props.append(classified_reads / total_reads)
        else:
            all_sample_family_unclassified_props.append(np.nan)
            all_sample_family_classified_props.append(np.nan)
    # --- Host species-level (across all samples) family classified/unclassified ---
    total_reads_family = df[sample_cols].values.sum()
    unclassified_reads_family = df.loc[mask_family_unclassified, sample_cols].values.sum()
    classified_reads_family = df.loc[mask_family_classified, sample_cols].values.sum()
    if total_reads_family > 0:
        all_species_family_unclassified_props.append(unclassified_reads_family / total_reads_family)
        all_species_family_classified_props.append(classified_reads_family / total_reads_family)
    else:
        all_species_family_unclassified_props.append(np.nan)
        all_species_family_classified_props.append(np.nan)

print("\nFinished processing all files.")

# Remove NaNs from the lists
all_sample_proportions = [x for x in all_sample_proportions if pd.notnull(x)]
all_species_proportions = [x for x in all_species_proportions if pd.notnull(x)]
all_sample_family_unclassified_props = [x for x in all_sample_family_unclassified_props if pd.notnull(x)]
all_sample_family_classified_props = [x for x in all_sample_family_classified_props if pd.notnull(x)]
# Remove NaNs for species-level family props
all_species_family_unclassified_props = [x for x in all_species_family_unclassified_props if pd.notnull(x)]
all_species_family_classified_props = [x for x in all_species_family_classified_props if pd.notnull(x)]

# --- Plot: Family-level classified vs. unclassified histograms (sample level) ---
plt.figure(figsize=(10, 6))
plt.hist(all_sample_family_unclassified_props, bins=30, color='orange', edgecolor='black', alpha=0.7, label='Unclassified at Family')
plt.hist(all_sample_family_classified_props, bins=30, color='green', edgecolor='black', alpha=0.5, label='Classified at Family')
plt.xlabel('Proportion of Reads')
plt.ylabel('Number of Samples')
plt.title('Distribution of Classified vs. Unclassified Reads at Family Level (per Sample)')
plt.legend()
plt.tight_layout()
plt.show()

# --- Plot: Family-level classified vs. unclassified histograms (host species level) ---
plt.figure(figsize=(10, 6))
plt.hist(all_species_family_unclassified_props, bins=30, color='orange', edgecolor='black', alpha=0.7, label='Unclassified at Family')
plt.hist(all_species_family_classified_props, bins=30, color='green', edgecolor='black', alpha=0.5, label='Classified at Family')
plt.xlabel('Proportion of Reads')
plt.ylabel('Number of Host Species')
plt.title('Distribution of Classified vs. Unclassified Reads at Family Level (per Host Species)')
plt.legend()
plt.tight_layout()
plt.show()

# --- Plot: Sample-level proportions (original, genus-based) ---
plt.figure(figsize=(10, 6))
plt.hist(all_sample_proportions, bins=30, color='skyblue', edgecolor='black')
plt.xlabel('Sample-level Proportion Hits (Unclassified Below Family)')
plt.ylabel('Number of Samples')
plt.title('Distribution of Sample-level Proportion Hits Across All Samples')
plt.tight_layout()
plt.show()

# --- Plot: Species-level proportions (original, genus-based) ---
plt.figure(figsize=(10, 6))
plt.hist(all_species_proportions, bins=30, color='salmon', edgecolor='black')
plt.xlabel('Species-level Proportion Hits (Unclassified Below Family)')
plt.ylabel('Number of Host Species')
plt.title('Distribution of Species-level Proportion Hits Across All Host Species')
plt.tight_layout()
plt.show()
