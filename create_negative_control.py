import sys
import csv

# ================= CONFIGURATION =================
# Inputs
blacklist_files = [
    "blacklist_uniprot_2.txt",
    "blacklist_biogrid.txt",
    "blacklist_string.txt"
]
homology_file = "homology_self.m8"
target_pairs_file = "uniprot_protein_pair_full.tsv"

# Output
output_file = "final_negative_control_set.txt"

# Homology Threshold (Safety check)
MIN_HOMOLOGY_ID = 40.0 
# =================================================

print("--- Phase 1: Loading Blacklists ---")
forbidden_set = set()

for b_file in blacklist_files:
    print(f"Loading {b_file}...")
    with open(b_file, 'r') as f:
        # Use simple split for speed, assuming tab-separated
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                # Store sorted tuple for symmetry
                pair = tuple(sorted((parts[0], parts[1])))
                forbidden_set.add(pair)

print(f"Total Unique Forbidden Pairs: {len(forbidden_set)}")


print("--- Phase 2: Loading Homology Map ---")
# Dict: protein -> set of homologs (including itself)
homologs = {}

with open(homology_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        # row: query, target, pident
        p1, p2 = row[0], row[1]
        
        # Clean IDs if they look like "sp|P12345|..."
        if "|" in p1: p1 = p1.split("|")[1]
        if "|" in p2: p2 = p2.split("|")[1]
        
        # Initialize if needed
        if p1 not in homologs: homologs[p1] = set()
        if p2 not in homologs: homologs[p2] = set()
        
        # Add to sets
        homologs[p1].add(p2)
        homologs[p2].add(p1) # Symmetry

print(f"Loaded homologs for {len(homologs)} proteins.")


print("--- Phase 3: Filtering the 200M Pairs (This will take time) ---")
saved_count = 0
dropped_count = 0

with open(target_pairs_file, 'r') as fin, open(output_file, 'w') as fout:
    # Buffer writing for speed
    writer = csv.writer(fout, delimiter='\t')
    
    for line in fin:
        parts = line.strip().split('\t')
        if len(parts) < 2: continue
        
        prot_a = parts[0]
        prot_b = parts[1]
        
        # 1. Direct Check (Fastest)
        direct_pair = tuple(sorted((prot_a, prot_b)))
        if direct_pair in forbidden_set:
            dropped_count += 1
            continue
            
        # 2. Homology Check 
        # Get homologs (default to just self if not in map)
        h_list_a = homologs.get(prot_a, {prot_a})
        h_list_b = homologs.get(prot_b, {prot_b})
        
        is_homologous_interaction = False

        for h_a in h_list_a:
            for h_b in h_list_b:
                # Check if this cousin-pair is forbidden
                h_pair = tuple(sorted((h_a, h_b)))
                if h_pair in forbidden_set:
                    is_homologous_interaction = True
                    break
            if is_homologous_interaction:
                break
        
        if is_homologous_interaction:
            dropped_count += 1
        else:
            writer.writerow([prot_a, prot_b])
            saved_count += 1
            
        if (saved_count + dropped_count) % 1000000 == 0:
            print(f"Processed: {saved_count + dropped_count/1e6:.1f}M | Saved: {saved_count} | Dropped: {dropped_count}")

print(f"DONE. Final Set Size: {saved_count}")
print(f"Total Dropped: {dropped_count}")
