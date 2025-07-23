import pandas as pd
import glob
import os

# Set the directory containing your merged output CSVs
input_dir = "family_level_output_archaea"  
# For family-level analysis:
taxon_col = "family"
median_col = "GeoMeanVolume_family_median"
# For genus-level analysis:
# input_dir = "genus_level_output_archaea"
# taxon_col = "genus"; median_col = "GeoMeanVolume_genus_median"

all_files = glob.glob(os.path.join(input_dir, "*.csv"))
unique_taxa = set()

for file in all_files:
    try:
        df = pd.read_csv(file)
        if median_col in df.columns and taxon_col in df.columns:
            matched = df[df[median_col].notna()]
            unique_taxa.update(matched[taxon_col].dropna().unique())
    except Exception:
        continue

print(f"Total number of unique microbial {taxon_col}s with cell size data: {len(unique_taxa)}")