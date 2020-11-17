#!/bin/bash

# Generate by pyGeneticPipe/shell/ShellMaker.py

# Setup batch control for slurm
#SBATCH --job-name=split_bed_by_chromosome
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:05:00
#SBATCH --mem=100m

# Source: zx8754 - https://www.biostars.org/p/387132/ 

# This pyGeneticPipe job takes the following direct args:
Plink_Path=apps/plink/2.00
Load_File=EUR.ldpred

# Load shell modules required to run
module load $Plink_Path

# For each chromosome within the file specified create a new file
for chr in {1..23}; do \
plink2 --bfile $Load_File --chr $chr --make-bed --out ${Load_File}_${chr}; \
done
