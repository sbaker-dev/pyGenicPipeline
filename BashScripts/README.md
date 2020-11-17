# What this directory Contains

Within this directory are example scripts that have been produced at one point or another which can be edited and used 
should you wish.

If you are using slurm, then look within the slurm folder. This will be updated with qsub at some point.

Within the slurm folder you will find the bash files, but these where generate on in a windows environment. This means 
when you copy these across to your linux environment you need to run both chmod to allow it run but also dos2unix. 

```shell script
chmod +x ScriptName.sh
dos2unix ScriptName.sh
```

Most scripts are generate with direct and indirect variables. Indirect variable simply put the value that was assigned
in line. For example the amount of time to run in an sbatch will be shown as --time:0:05:00. Direct variables that are
to do with the actual program will be below a comment of 'This pyGeneticPipe Job takes...'. If you chang ethe variables
within there then, condition on your file being of equivalent size to the test data, they should run without error.

Most of these files take use files provided by [This tutorial][tutpath] for testing purposes. If you are running this on
actual UKBiobank or other full Genetic samples you will need to adjust your times, memory allocation and potential 
nodes accordingly. If these files have been tested they will be logged here as to the parameters that where used. 

UkBiobank: 


[tutpath]: https://choishingwan.github.io/PRS-Tutorial/target/