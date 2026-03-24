import csv
import sys

# Increase CSV field size limit for massive BioGRID lines
csv.field_size_limit(sys.maxsize)

# CONFIGURATION
input_file = "BIOGRID-ORGANISM-Homo_sapiens-5.0.252.tab3.txt"
output_file = "../blacklist_biogrid.txt"

print(f"Processing {input_file}...")

forbidden_pairs = set()

# Columns to target (Copied exactly from your header)
COL_A = "SWISS-PROT Accessions Interactor A"
COL_B = "SWISS-PROT Accessions Interactor B"

with open(input_file, 'r') as f:
    # DictReader automatically uses the first row (Header) to map keys
    reader = csv.DictReader(f, delimiter='\t')
    
    for row in reader:
        # 1. Get the Raw Strings
        raw_a = row[COL_A]
        raw_b = row[COL_B]
        
        # 2. Skip if empty or just a dash "-"
        if raw_a == "-" or raw_b == "-":
            continue
            
        # 3. Handle Multiple IDs (Split by pipe '|')
        ids_a = raw_a.split("|")
        ids_b = raw_b.split("|")
        
        # 4. Generate All Combinations
        for protein_a in ids_a:
            p_a = protein_a.strip()
            
            for protein_b in ids_b:
                p_b = protein_b.strip()
                
                # Exclude self-interactions
                if p_a == p_b:
                    continue
                
                # Sort and Store
                pair = tuple(sorted((p_a, p_b)))
                forbidden_pairs.add(pair)

print(f"Writing {len(forbidden_pairs)} pairs to {output_file}...")

with open(output_file, "w") as f:
    for p1, p2 in forbidden_pairs:
        f.write(f"{p1}\t{p2}\n")

print("Step 2 Complete. BioGRID Blacklist created.")
