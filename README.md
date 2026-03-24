# Protein-Protein-Interaction-Negative-Control-Set

## Overview
This repository contains the data processing pipeline for generating a rigorously filtered, high-confidence negative control set of protein-protein interactions (PPIs). This dataset serves as a strict baseline for benchmarking the performance of the **PPLM-PPI** model against **RF2-PPI** in binary interaction prediction.

To accurately evaluate PPLM-PPI's ability to recover interactions lacking strong co-evolutionary signals, the negative control set must be absolutely devoid of true positives. This pipeline achieves this by aggressively filtering potential negative pairs against multiple comprehensive interaction databases.

## Data Sources (Prerequisites)
To run this pipeline, you must acquire the following raw datasets and place them in the working directory:
* **STRING Database (v12.0):** `9606.protein.links.v12.0.txt.gz`, `9606.protein.physical.links.v12.0.txt.gz`, `9606.protein.sequences.v12.0.fa`
* **UniProt:** `uniprot_protein_pair_full.tsv`, `uniprotwithinteracts.tsv`
* **BioGRID:** Known interactions downloaded from the BioGRID database.
* **Homology Data:** Pre-computed `.m8` files (e.g., `string_to_uniprot.m8`, `homology_self.m8`) for ID mapping and paralog filtering.

## Pipeline Architecture
The pipeline executes in distinct filtering stages to generate blacklists of known interacting pairs, ensuring they are excluded from the final negative set:

1.  **STRING Filtering (`generate_string_blacklist.py`):** Parses the massive STRING network and physical links databases to identify and blacklist all protein pairs with known functional or physical associations. Outputs `blacklist_string.txt`.
2.  **UniProt Filtering (`generate_blacklist_uniprot.py`):** Processes UniProt interaction annotations to capture experimentally verified pairs. Outputs `blacklist_uniprot_2.txt`.
3.  **BioGRID Filtering (`generate_blacklist_biogrid.py`):** Parses the local `BIOGRID` files to dynamically extract and blacklist additional curated interactions. Outputs `blacklist_biogrid.txt`.
4.  **Final Assembly (`create_negative_control.py`):** Aggregates all blacklists, applies homology mapping to prevent biased sampling, and generates the definitive `final_negative_control_set.txt` for benchmarking.

## Usage
The entire pipeline is designed to be executed on a high-performance computing cluster using the SLURM workload manager.

Submit the pipeline script to your cluster:
```bash
sbatch create_negative_control.slurm
