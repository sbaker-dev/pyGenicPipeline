#!/bin/bash

# Generate by pyGeneticPipe/shell/ShellMaker.py

# Setup batch control for slurm
#SBATCH --job-name=convert_to_bgen
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:05:00
#SBATCH --mem=100m

# See other conversions from https://www.well.ox.ac.uk/~gav/qctool/documentation/examples/converting.html 

# This pyGeneticPipe job takes the following direct args:
QCTool_Path=apps/qctool/2.0rc4
Load_File=EUR.ldpred

# Load shell modules required to run
module load $QCTool_Path

# For each .bed chromosome file within this directory convert it to .bgen v1.2
for chr in {1..23}; do \
qctool -g ${Load_File}_${chr}.bed -og ${Load_File}_${chr}.bgen; \
done
