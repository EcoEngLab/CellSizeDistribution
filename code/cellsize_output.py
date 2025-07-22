print("Starting script...")
from bacdive import BacdiveClient
client = BacdiveClient("ching.wan24@imperial.ac.uk", "Cm102056!")
print("Client setup complete.")

import pandas as pd
import re 
import math
import numpy as npy

try:
    df_cleaned = pd.read_csv('cleaned.archaea.raw.csv')
    bacids = df_cleaned['ID'].tolist()
    print("Successfully loaded BacDive IDs.")
    print("First 10 BacDive IDs:", bacids[:10])
    print("Total number of BacDive IDs:", len(bacids))

except FileNotFoundError:
    print("Please make sure the CSV file is in the same directory as the script.")
    bacids = [] 
except Exception as e:
    print(f"An error occurred while reading the CSV file: {e}")
    bacids = []

# Clean units using regex (remove anything not number, dash, dot, space), separate min and max for cell length and width
def parse_range(value):
    if value == "NA":
        return ("NA", "NA")
    if isinstance(value, str):
        value = value.strip()
        value = re.sub(r"[^\d\–\-\. ]", "", value)
        value = re.sub(r"[--]+", "-", value)
        parts = value.split("-")

        try:
            if len(parts) == 2:
                return (float(parts[0].strip()), float(parts[1].strip()))
            elif len(parts) == 1:
                val = float(parts[0].strip())
                return (val,val)
        except ValueError:
            return ("NA", "NA")
        return ("NA", "NA")
    
    return ("NA", "NA")

# volume and surface area calculations from jolan
def volumecalc(thisrow):
    '''
    Purpose:
        Calculate the minimum, maximum and average volume of the current species of bacteria or archea being iterated

    Parameters:
        thisrow: this is the current row in the csv file

    Returns:
        minVolume: minimum volume of cell
        maxVolume: maximum volume of cell
        avgVolume: average of the two
    '''
    pi = math.pi # Define pi here
    shape = thisrow["Cell shape"]
    minVolume = -1
    maxVolume = -1
    surfaceArea = -1
    minlengthnocaps = 0
    maxlengthnocaps = 0
    try:
        # These are the column names the function expects
        minlengthtotal = thisrow["LengthMin"]
        maxlengthtotal = thisrow["LengthMax"]
        minwidth = thisrow["WidthMin"]
        maxwidth = thisrow["WidthMax"]

        # Convert to numeric, as the function expects floats/ints
        minlengthtotal = pd.to_numeric(minlengthtotal, errors='coerce')
        maxlengthtotal = pd.to_numeric(maxlengthtotal, errors='coerce')
        minwidth = pd.to_numeric(minwidth, errors='coerce')
        maxwidth = pd.to_numeric(maxwidth, errors='coerce')

        # Check if any conversions failed
        if pd.isna(minlengthtotal) or pd.isna(maxlengthtotal) or pd.isna(minwidth) or pd.isna(maxwidth):
             raise ValueError("Invalid number format")

    except Exception: # Catches missing keys or conversion errors
        # This will now return the expected Series with NaNs
        return pd.Series({"VolumeMin":npy.nan, "VolumeMax":npy.nan, "SurfaceArea":npy.nan})

    minradius = minwidth/2
    maxradius = maxwidth/2
    if (str(shape).lower() in ["rod-shaped", "filament-shaped", "vibrio-shaped", "helical-shaped", "spiral-shaped", "curved-shaped", "spore-shaped"]):
        minlengthnocaps = minlengthtotal - minwidth
        maxlengthnocaps = maxlengthtotal - maxwidth
        if (minlengthnocaps < 0):
            minlengthnocaps = minlengthnocaps*-1
            minwidth = minlengthtotal
        if (maxlengthnocaps < 0):
            maxlengthtotal = maxlengthtotal*-1
            maxwidth = maxlengthtotal
        if (minlengthtotal == minwidth):
            minlengthnocaps = 0
        if (maxlengthtotal == maxwidth):
            maxlengthnocaps = 0 # Corrected a small typo here from '==' to '='
        minVolume = ((pi*minradius*minradius)*(((4/3)*minradius)+minlengthnocaps))
        maxVolume = ((pi*maxradius*maxradius)*(((4/3)*maxradius)+maxlengthnocaps))
        surfaceArea = npy.sqrt(((2*pi*minradius)*((2*minradius) + minlengthnocaps))*((2*pi*maxradius)*((2*maxradius) + maxlengthnocaps)))
    elif (str(shape).lower() in ["ellipsoidal", "coccus-shaped", "ovoid-shaped", "oval-shaped", "sphere-shaped", "crescent-shaped"]):
        minVolume = ((4/3)*pi*minradius*minradius*(minlengthtotal/2))
        maxVolume = ((4/3)*pi*maxradius*maxradius*(maxlengthtotal/2))
        surfaceArea = npy.sqrt((4*pi*minradius*minradius)*(4*pi*maxradius*maxradius))
    if (surfaceArea > 0 and minVolume > 0 and maxVolume > 0):
        return pd.Series({"VolumeMin":minVolume, "VolumeMax":maxVolume, "SurfaceArea":surfaceArea})
    else:
        return pd.Series({"VolumeMin":npy.nan, "VolumeMax":npy.nan, "SurfaceArea":npy.nan})

# retrieve cell morphology data from bacdive
results = []
#bacids = [141032, 132482, 133021, 140795, 133245, 141136, 17820] 
for bacid in bacids:
    #print(f"BacDive ID = {bacid}")
    client.search(id=int(bacid))
    for record in client.retrieve():
        gen = record.get("Name and taxonomic classification", {})
        morph = record.get("Morphology", {})
        seq = record.get("Sequence information", {})

        tax_id = {}  
        s16_data = seq.get("16S sequences")
        if isinstance(s16_data, list) and len(s16_data) > 0:
            #tax_id = s16_data[0].get("NCBI tax ID", "NA")
            for i, taxid_entry in enumerate(s16_data):
                if "NCBI tax ID" in taxid_entry:
                    tax_id = s16_data[i]
                    break
                else:
                    tax_id = {}
        elif isinstance(s16_data, dict):
            tax_id = s16_data

        if tax_id == {}:
            genome_data = seq.get("Genome sequences")
            if isinstance(genome_data, list) and len(genome_data) > 0:
                #tax_id = genome_data[0].get("NCBI tax ID", "NA")
                for j, taxid_entry in enumerate(genome_data):
                    if "NCBI tax ID" in taxid_entry:
                        tax_id = genome_data[j]
                        break
                    else:
                        tax_id = {}
            elif isinstance(genome_data, dict):
                tax_id = genome_data

        cell_morph_data = morph.get("cell morphology")
        final_cell_morph = {}
        shape_final_cell_morph = {}

        if isinstance(cell_morph_data, list) and len(cell_morph_data) > 0:
            for i, morph_entry in enumerate(cell_morph_data):
                if "cell length" in morph_entry and "cell width" in morph_entry:
                    final_cell_morph = cell_morph_data[i]
                    break
                else:
                    final_cell_morph = {}
            for j, morph_entry in enumerate(cell_morph_data):
                if "cell shape" in morph_entry:
                    shape_final_cell_morph = cell_morph_data[j]
                    break
                else:
                    shape_final_cell_morph = {}

        elif isinstance(cell_morph_data, dict):
            final_cell_morph = cell_morph_data
            shape_final_cell_morph = cell_morph_data

        if final_cell_morph and isinstance(gen, dict):
            cell_species = gen.get("species", "NA")
            cell_genus = gen.get("genus", "NA")
            cell_family = gen.get("family", "NA")
            cell_order = gen.get("order", "NA")
            cell_class = gen.get("class", "NA")
            cell_phylum = gen.get("phylum", "NA")
            cell_domain = gen.get("domain", "NA")
            cell_length = final_cell_morph.get("cell length", "NA")
            cell_width = final_cell_morph.get("cell width", "NA")
            cell_shape = shape_final_cell_morph.get("cell shape", "NA")
            cell_taxonomyid = tax_id.get("NCBI tax ID", "NA")

        else:
            print("  No cell morphology data found")

        results.append({
            "BacDive ID": bacid,
            "Species": cell_species,
            "Genus": cell_genus,
            "Family": cell_family,
            "Order": cell_order,
            "Class": cell_class,
            "Phylum": cell_phylum,
            "Domain": cell_domain,
            "NCBI taxonomy ID": cell_taxonomyid,
            "Cell length (μm)": cell_length,
            "Cell width (μm)": cell_width,
            "Cell shape": cell_shape,
        })

df = pd.DataFrame(results)

df["Cell length min (μm)"], df["Cell length max (μm)"] = zip(*df["Cell length (μm)"].map(parse_range))
df["Cell width min (μm)"], df["Cell width max (μm)"] = zip(*df["Cell width (μm)"].map(parse_range))
df = df.drop(columns=["Cell length (μm)", "Cell width (μm)"]) # remove the original cell length and width columns

# Create a temporary DataFrame with column names that `volumecalc` expects.
temp_df = df.rename(columns={
    "Cell length min (μm)": "LengthMin",
    "Cell length max (μm)": "LengthMax",
    "Cell width min (μm)": "WidthMin",
    "Cell width max (μm)": "WidthMax",
})

# STEP 2: Apply the new, unchanged `volumecalc` function.
# This creates new columns: "VolumeMin", "VolumeMax", and "SurfaceArea".
calc_results = temp_df.apply(volumecalc, axis=1)

# STEP 3: Join the results back to the original DataFrame.
df = df.join(calc_results)

# STEP 4: Rename the new columns to the desired final, descriptive format.
df = df.rename(columns={
    "VolumeMin": "Cell volume min (μm³)",
    "VolumeMax": "Cell volume max (μm³)",
    "SurfaceArea": "Cell surface area (μm²)" 
})

# STEP 5: Add average and geometric mean calculations from the new code.
df["Average cell volume (μm³)"] = df[["Cell volume min (μm³)", "Cell volume max (μm³)"]].mean(axis=1)

# Calculate Geometric Mean only for valid rows
valid_rows = (df["Cell volume min (μm³)"] > 0) & (df["Cell volume max (μm³)"] > 0)
df.loc[valid_rows, "Geometric mean cell volume (μm³)"] = npy.sqrt(df.loc[valid_rows, "Cell volume min (μm³)"] * df.loc[valid_rows, "Cell volume max (μm³)"])
df["Geometric mean cell volume (μm³)"].fillna(npy.nan, inplace=True) # Fill the rest with NaN

print(df)
#print(df)

df.to_csv("archaea.VolumeAndSurafaceArea.csv", index=False)