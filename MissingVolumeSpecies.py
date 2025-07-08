import pandas as pd

# Replace with your actual CSV filename
csv_filename = 'OutputBacteria.csv'
output_txt_filename = 'bacteria_missing_volume.txt'

# Read the CSV file
try:
    df = pd.read_csv(csv_filename)
except FileNotFoundError:
    print(f"File not found: {csv_filename}")
    exit(1)

# Find rows where 'GeoMeanVolume' is missing (empty or NaN)
missing_volume = df[df['GeoMeanVolume'].isnull() | (df['GeoMeanVolume'].astype(str).str.strip() == '')]

# Get the species names as a list
species_missing = missing_volume['Species'].dropna().astype(str).tolist()

# Print the species names
print("Species with missing GeoMeanVolume:")
for species in species_missing:
    print(species)

# Write the species names to a text file
with open(output_txt_filename, 'w') as f:
    for species in species_missing:
        f.write(species + '\n')

print(f"\nList of species with missing GeoMeanVolume written to {output_txt_filename}")
