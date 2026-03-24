#!/bin/bash
#SBATCH --job-name=filter_negative_controls      
#SBATCH --output=filter_negative_controls.out  
#SBATCH --error=filter_negative_controls.err    
#SBATCH --time=24:00:00                         
#SBATCH -p Normal                                 
#SBATCH -N 1                                     
#SBATCH --ntasks-per-node=1                     
#SBATCH --cpus-per-task=4                        
#SBATCH --mem=16G                                

module load python

python generate_final_negative_control_set.py

