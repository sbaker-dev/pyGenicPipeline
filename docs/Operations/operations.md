---
layout: default
title: Operations
nav_order: 3
has_children: true
---

# Operations

pyGenicPine is designed to combine as many packages and methods into a single pipeline for ease of use and standardised
documentation. Each of the sections below represent groups of operations that you can run, or a pipeline with steps such
as making a polygenic score.

## Shell

Whilst ever effort has been made to make files cross compatible with differing file structures, for example allowing
both plink and bgen files to work with polygenic scores, sometimes pipelines may require a set file type or structure.
This group of operations gives you the command you need to run common requirements for this pipeline or others.

## Polygenic Scores

This set of operations are used to construct polygenic scores . This pipeline allows for both plink and bgen files to
be used and can be broken down into multiple steps so you can use/validate output from any of the steps individually as
well as the final score itself.

# Configuring 

Preferences on how to use the pipeline will likely come from a combination of familiarity to similar processors you 
already use and your computing background. Where ever possible this pipeline is designed to allow non-computer 
scientists to use. Some processors may not be possible to run on a local machine, so a familiarity with bash will be 
helpful as you will likely have to run most of these processes via bash.

There are three core ways you can set Arguments. Natively within python, by editing a yaml file, and by using the GUI.
For these examples we will be running the [split_bed_by_chromosome](Shell/split%20bed%20by%20chromosome.md). In all 
cases we will be using data that was download from [this tutorial][exd], but you can download this which has been 
pre-split into chromosomes for you from this [dropbox directory][DBD].

## 1) Python
If you have a familiarity with python then you can use the packages ArgMaker class to initialise the require arguments
and then set them within python itself. The arguments take a form of a dict which will be used to write a yaml file, 
with yaml being chosen over json due to the ability to add comments. 

```python

from pyGenicPipeline import ArgMaker
from pathlib import Path

args_maker = ArgMaker(Path(r"Path/defaults.yaml"))

# Printing this will lead to output seen below
print(args_maker.convert_to_bgen)

output = {
    "Mandatory": 
    {
    "Working_Directory": None,
    "Load_File": None,
    "Operation": "convert_to_bgen",
    "Job_Scheduler": "slurm",
    "Partition": "test",
    'Plink_Path': 'apps/plink/2.00'
    }, 
    "Optional":
    {
    "Nodes": 1,
    "Job_Time": "00:05:00",
    "Job_Memory": "1000m"
    }
}
```

As you can see Jobs have mandatory args which you must set as well as optional ones which you may need or what to 
change. For our example case a 5 minute wall time is fine, but if your converting a 100TB file then its unlikely to be
sufficient. 

Whilst certain arguments are likely to be project specific such as the working directory many are likely variables you
only would want to set once. For example, here you will see that the the Job Scheduler has been set to slurm as this is
the system of the super computer i have access to. By providing a path to this file to the ArgMaker, these defaults will
automatically be filled in for you but this is optional and you do not need to provide any arguments to ArgMaker by 
default. 

You can also use this to save time. For example the Partition is set to test as for the majority of scripts that are 
generated they can be run in under 1 hour which is the maximum wall time for the test directory. However, lets say we 
needed this for 99% of the time but we need to update it for this process. Because the args are a dict, you can just
update it based on the keys.

```python

from pyGenicPipeline import ArgMaker
from pathlib import Path

args_maker = ArgMaker(Path(Path(__file__).parent, "defaults.yaml"))

args_dict = args_maker.split_bed_by_chromosome

args_dict["Mandatory"]["Working_Directory"] = "Working_Directory"

args_maker.write_yaml_args(args_dict, r"write_out_path")
```

Based on the defaults that where used, here you can see we just need to update the working directory before we can use 
it. At this point we now need to compile our args into a yaml config file. All you need to do is pass the dict you just
made as an the first argument, and the path to save it to as the second and then your finished constructing the 
arguments.

## 2) Yaml

If python isn't a language you are comfortable and just want to be able to have your arguments to submit to a sever then
you can also manually edit the yaml files. Within the SubmissionScripts folder of this github directory are all the 
bash submission files for Shell scripts such as split_bed_by_chromosome.sh and the yaml files for pyGenicPipeline to run
on the sever. These are the blank outputs of what ArgsMaker would produce in our previous example.

Yaml files act as dictionary and the syntax is quite simple. For example lets take the snippet of the
 split_bed_by_chromosome.yaml. The Orange Words, or words prefacing the colon, represent the variable name. You **should
not** edit these, if you do the system will crash. What you want to edit is the after the colon, so here you could 
copy the path to the working directory and paste it where the null is instead.

Config files for pyGenicPipeline have three types of arguments. Mandatory which should always be set, Optional which
can be left as null or at defaults but you may wish to change, and System which you can and should leave as null. Not
all processes need all the arguments, but during initialisation to determine what needs to be set from validating which
of the system arguments are validated as None, so if you delete them the system won't be able to run.

```yaml
# Project or output directory
Working_Directory: null
```

[exd]: https://choishingwan.github.io/PRS-Tutorial/target/
[DBD]: https://www.dropbox.com/sh/3uust1mcv10dj6l/AACkHUpjw9QgX_XqkzNZ7Wyva?dl=0

## 3) GUI

There is also a GUI you can use if you want a more 'software' like interface for editing yaml files. Each step is 
explained in more detail akin to the process Steps within the individual operation pages. See the setup instructions
within setup to install and use this on your PC.