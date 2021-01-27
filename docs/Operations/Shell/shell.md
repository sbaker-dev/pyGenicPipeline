---
layout: default
title: Shell
parent: Operations
has_children: true
nav_order: 1
---

## Shell Commands

These are the shell commands that may be useful or are required to format or structure data and output for this 
pipeline. Unlike many of the python commands, where bash is used as the submission tool alone, these scripts can also
be found in the github repository of SubmissionScripts/Slurm/Bash

### convert_to_bgen

The polygenic pipeline is now able to read bgen files, so strictly speaking this isn't required anymore. However, if you
are using other software that requires bgen then you can still use this

### split_bed_by_chromosome

Some genetic data may come single bed file, but the majority of this pipeline **requires** all genetic file structures
to be split by chromosome. This facilitates multi-core/processing to be facilitated, reduces the looping requirements
within the source code, and can act as a significant speed/memory saving measure.