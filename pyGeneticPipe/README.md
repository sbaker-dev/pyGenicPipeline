### Code breakdown 

pyGeneticPipe has lots of separate code both original and re-worked from pre-existing repositories. This file should 
help explain what each file does so you can find the source code you are looking for and also provide links back to the
original sources should you want the unmodified versions.

### bgen

This contains a modified version of [pybgen][pbgen] that was altered slightly, not because there was anything wrong with
it, but so that it fit my own code style so it would be easier to debug and also add features. 

###### bgenObject.py

This contains the class of BgenObject that can parse in a file header of a bgen file and then use its properties.

### plink

This contains a modified version of [plinkio][pio] due to the linux binaries not compiling on one of the older super 
computers and it not being cross compatible.  

###### plinkObject

This file is used as a controlling for working with a set of plink files (bed, bim, fam).

###### supportObjects

Many files are set in rows, groups, sets or some for of object that needs to be reference. supportObjects contains 
objects for plinkObject to use, these are designed to the reference multiple times per file or operation and are 
unlikely to be needed by an end user. However, they are kept separate to avoid nest classes in case they are required.

supportObjects contains:

* ***BimLoci***: This contains a row a Bim file, which can be accessed via keys
* ***BimByChromosome***: This contains all the information from the bim, but formatted by chromosome


### PRS

This folder contains the code for the creation of polygenic scores and is a modified approach of [LDPred][ldp] 1 which
was coded in python to work with bgen as well as some other modifications.

###### Summary.py

This represents the class for managing cleaning of GWAS summary statistics.

###### prsConstructor.py

This represents the class that calls all other sub-classes to construct the PRS


### utils

This folder contains code that may be used across programs, non-class specific methods that would otherwise be static,
as well as a storage space for anything miscellaneous.

###### error_codes.py

A separate file for error codes which are called if assert condition X fails. 

###### Input.py

This class loads the inputted args and validates them when required via asserting they are applicable. The args are then
stored as class variables within the constructor which can be called via inheritance from other classes. 

###### misc.py

Contains a few methods that would otherwise be static within classes, if sufficient similarities between arguments of 
these methods are found they will be spun out into a new class within utils.  




[pbgen]: https://github.com/lemieuxl/pybgen
[ldp]: https://github.com/bvilhjal/ldpred
[pio]: https://github.com/mfranberg/libplinkio