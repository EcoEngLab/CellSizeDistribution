"""
This script retrieves cell morphology and taxonomic data for a list of bacterial species from the BacDive database using the BacdiveClient API.

- Reads BacDive IDs from 'cleaned.bacteria.raw.csv'.
- Fetches taxonomy and cell morphology (length, width, shape) for each ID.
- Processes and cleans the data, extracting min/max cell dimensions.
- Saves the results to 'bacteria.bacdive_cell_data.ver2.csv'.

BacDive login credentials must be provided via the BACDIVE_EMAIL and BACDIVE_PASSWORD environment variables for security.
"""

import os
import pandas as pd
import re 
import math

# Set up BacDive client using credentials from environment variables
from bacdive import BacdiveClient
bacdive_email = os.environ.get("BACDIVE_EMAIL")
bacdive_password = os.environ.get("BACDIVE_PASSWORD")
if not bacdive_email or not bacdive_password:
    raise ValueError("Please set BACDIVE_EMAIL and BACDIVE_PASSWORD environment variables.")
client = BacdiveClient(bacdive_email, bacdive_password)

# Load BacDive IDs from the cleaned CSV file
try:
    df_cleaned = pd.read_csv('cleaned.bacteria.raw.csv')
    bacids = df_cleaned['ID'].tolist()
    # print("Successfully loaded BacDive IDs.")
    # print("First 10 BacDive IDs:", bacids[:10])
    # print("Total number of BacDive IDs:", len(bacids))
except FileNotFoundError:
    # print("Please make sure the CSV file is in the same directory as the script.")
    bacids = [] 
except Exception as e:
    # print(f"An error occurred while reading the CSV file: {e}")
    bacids = []

# Helper function to parse cell length/width ranges from strings
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

# Retrieve cell morphology and taxonomy data from BacDive for each ID
results = []
#bacids = [141032, 132482, 133021, 140795, 133245, 141136, 17820] 
for bacid in bacids:
    client.search(id=int(bacid))
    for record in client.retrieve():
        gen = record.get("Name and taxonomic classification", {})
        morph = record.get("Morphology", {})
        seq = record.get("Sequence information", {})

        # Extract NCBI taxonomy ID from sequence data
        tax_id = {}  
        s16_data = seq.get("16S sequences")
        if isinstance(s16_data, list) and len(s16_data) > 0:
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
                for j, taxid_entry in enumerate(genome_data):
                    if "NCBI tax ID" in taxid_entry:
                        tax_id = genome_data[j]
                        break
                    else:
                        tax_id = {}
            elif isinstance(genome_data, dict):
                tax_id = genome_data

        # Extract cell morphology data
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

        # Collect all relevant data for this BacDive ID
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

# Convert results to DataFrame for further processing
df = pd.DataFrame(results)

# Parse and split cell length and width into min and max columns
df["Cell length min (μm)"], df["Cell length max (μm)"] = zip(*df["Cell length (μm)"].map(parse_range))
df["Cell width min (μm)"], df["Cell width max (μm)"] = zip(*df["Cell width (μm)"].map(parse_range))
df = df.drop(columns=["Cell length (μm)", "Cell width (μm)"]) # remove the original cell length and width columns

# Print the final DataFrame and save to CSV
print(df)
df.to_csv("bacteria.bacdive_cell_data.ver2.csv", index=False)
