---
layout: default
title: Split Bed by Chromosome
parent: Shell
grand_parent: Operations
nav_order: 2
---

## Split Bed by Chromosome
Some genetic data may come single bed file, but the majority of this pipeline **requires** all genetic file structures
to be split by chromosome. This facilitates multi-core/processing to be facilitated, reduces the looping requirements
within the source code, and can act as a significant speed/memory saving measure.

### Required 

The following is a snippet of the yaml file was create via ArgMaker. These are the variables that will be required to 
be set. You can create this uses ArgMaker in python as described in operations, or download this yaml config file to pass
to the pipeline [here][YD].

```yaml
# Project or output directory
Working_Directory: null

# The operation to be run, can take a string if the method name or a dict of operations where only one element is true
Operation: split_bed_by_chromosome

# Some operations will require a path to a set file
Load_File: null

# Takes the values of slurm, qsub, and null
Job_Scheduler: null

# Partition to use when using slurm
Partition: null

# Path to plink on sever, has a default based on bc4 bristol
Plink_Path: null



# OPTIONAL ARGS - NOT STRICTLY REQUIRED TO RUN BUT ALTERS OPERATION ####################################################

# Number of nodes to use, defaults to 1
Nodes: 1

# Job time via Job_Scheduler
Job_Time: 00:05:00

# Amount of memory to allocate via Job_Memory, defaults to 100m
Job_Memory: 100m


```

### Shell

This will be passed to Shell, which will configure these arguments into a bash file for you. You can download this file 
[here][ShellD] and just alter the bash file directly

```console

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
```

### Process Flow

In this case all pyGenicPipe does is create the submission file for you, so there is no process flow to describe.


[QC]: https://www.well.ox.ac.uk/~gav/qctool_v2/
[BG]: https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk
[YD]: https://github.com/sbaker-dev/pyGenicPipeline/blob/main/SubmissionScripts/Yaml/split_bed_by_chromosome.yaml
[ShellD]: https://github.com/sbaker-dev/pyGenicPipeline/blob/main/SubmissionScripts/Bash/split_bed_by_chromosome.sh

