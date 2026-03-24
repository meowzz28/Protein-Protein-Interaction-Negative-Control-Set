import pandas as pd
import csv

# ================= CONFIGURATION =================
mapping_file = "string_to_uniprot.m8"
links_file_functional = "string_protein_links.tsv"
links_file_physical = "string_physical_links.tsv"
output_file = "blacklist_string.txt"

# Thresholds 
MIN_IDENTITY = 95.0
MIN_COVERAGE = 0.90
# =================================================

print("--- Step 1: Loading and Filtering Sequence Map ---")

# Dictionary: string_id -> set(uniprot_ids)
id_map = {}

with open(mapping_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:

        if not row: continue

        string_id = row[0]
        
        # Use the ID directly
        raw_target = row[1].strip() 
        
        pident = float(row[2])
        qcov = float(row[3])
        tcov = float(row[4])

        # 1. Apply  Filters
        if pident < MIN_IDENTITY:
            continue
        if qcov < MIN_COVERAGE or tcov < MIN_COVERAGE:
            continue

        # 2. Add to map
        if string_id not in id_map:
            id_map[string_id] = set()
        
        # If the ID has versioning (like P12345.1), strip it just in case. 
        clean_id = raw_target.split(".")[0]
        id_map[string_id].add(clean_id)

print(f"Mapped {len(id_map)} STRING proteins to UniProt IDs.")

# =================================================
print("--- Step 2: Processing Interaction Files ---")

forbidden_pairs = set()

def process_string_file(filepath):
    count = 0
    with open(filepath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        
        # Read first row to check for header
        first_row = next(reader)
        
        # If the first row starts with "protein", it's a header. Skip it.
        # If it looks like data (starts with '9606.'), process it.
        if first_row[0].startswith("protein") or first_row[0].startswith("#"):
            pass 
        else:
            # It was data, process this row first!
            process_row(first_row)
        
        # Process the rest
        for row in reader:
            process_row(row)

def process_row(row):
    if len(row) < 2: return

    s_id_a = row[0]
    s_id_b = row[1]
    
    # 1. Check if we have valid maps for BOTH proteins
    if s_id_a in id_map and s_id_b in id_map:
        
        # 2. Get all UniProt possibilities
        uni_a_list = id_map[s_id_a]
        uni_b_list = id_map[s_id_b]
        
        # 3. Create combinations
        for u_a in uni_a_list:
            for u_b in uni_b_list:
                if u_a == u_b: continue
                
                # Sort and Add
                pair = tuple(sorted((u_a, u_b)))
                forbidden_pairs.add(pair)

# Process Functional Links
print(f"Reading {links_file_functional}...")
process_string_file(links_file_functional)

# Process Physical Links
print(f"Reading {links_file_physical}...")
process_string_file(links_file_physical)

# =================================================
print(f"--- Step 3: Writing {len(forbidden_pairs)} Unique Forbidden Pairs ---")

with open(output_file, "w") as f:
    for p1, p2 in forbidden_pairs:
        f.write(f"{p1}\t{p2}\n")

print("Done. STRING Blacklist created.")
