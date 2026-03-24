import pandas as pd
import csv

# CONFIGURATION
input_file = "uniprotwithinteracts.tsv" 
output_file = "blacklist_uniprot_2.txt"

print(f"Processing {input_file}...")

# Robust cleaning function
def clean_uniprot_id(raw_id):
    # Case 1: "PRO_0000037548 [Q9WMX2]" -> Extract Q9WMX2
    if "[" in raw_id and "]" in raw_id:
        return raw_id.split("[")[1].split("]")[0].strip()
    
    # Case 2: "P12345-2" (Isoform) -> Extract P12345
    return raw_id.split("-")[0].strip()

df = pd.read_csv(input_file, sep='\t', quoting=csv.QUOTE_NONE, on_bad_lines='skip')
df_interactions = df[df["Interacts with"].notna()]
forbidden_pairs = set()

for index, row in df_interactions.iterrows():
    # Clean Protein A
    protein_a = clean_uniprot_id(row["Entry"])
    
    partners_list = str(row["Interacts with"]).split(";")
    
    for partner in partners_list:
        # Clean Protein B (This handles the PRO_... [ID] case)
        protein_b = clean_uniprot_id(partner.strip())
        
        if protein_a == protein_b:
            continue
        
        pair = tuple(sorted((protein_a, protein_b)))
        forbidden_pairs.add(pair)

print(f"Writing {len(forbidden_pairs)} clean pairs to {output_file}...")

with open(output_file, "w") as f:
    for p1, p2 in forbidden_pairs:
        f.write(f"{p1}\t{p2}\n")

print("Correction Complete.")
