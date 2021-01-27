---
layout: default
title: Convert to Bgen
parent: Shell
grand_parent: Operations
nav_order: 1
---

## Convert to Bgen
Sometimes you may need to convert plink files to bgen file format. We can do with [QCTool][QC] and [BGENIX][BG] which 
you will need to have installed. 

### Required 

The following is a snippet of the yaml file was create via ArgMaker. These are the variables that will be required to 
be set. You can create this uses ArgMaker in python as described in operations, or download this yaml config file to 
pass to the pipeline [here][YD].

```yaml
# MANDATORY ARGS - IF ANY ARE NONE JOB WILL NOT RUN ####################################################################

# Project or output directory
Working_Directory: null

# The operation to be run, can take a string if the method name or a dict of operations where only one element is true
Operation: convert_to_bgen

# Some operations will require a path to a set file
Load_File: null

# Takes the values of slurm, qsub, and null
Job_Scheduler: null

# Partition to use when using slurm
Partition: null

# Path to qctoolv2 on sever, has a default based on bc4 bristol
QCTool_Path: null

# Path to the compiled bgenix folder in build/apps/bgenix
BGENIX_Path: null

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
```

### Process Flow

In this case all pyGenicPipe does is create the submission file for you, so there is no process flow to describe.


[QC]: https://www.well.ox.ac.uk/~gav/qctool_v2/
[BG]: https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk
[YD]: https://github.com/sbaker-dev/pyGenicPipeline/blob/main/SubmissionScripts/Yaml/convert_to_bgen.yaml
[ShellD]: https://github.com/sbaker-dev/pyGenicPipeline/blob/main/SubmissionScripts/Bash/convert_to_bgen.sh

