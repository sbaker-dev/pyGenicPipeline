---
layout: default
title: Conda
parent: Setup
nav_order: 2
---

## Conda

This is going to assume your working on a linux server. We will be using a python environment, setup via conda. If you 
have not built one before, then follow these instructions. Otherwise you just need an activate conda environment with 
the pyGenicPipe module and the bgen-reader forked form sbaker-dev to allow for custom meta data locations.

First we need to add python anaconda to our bash environment. If you have anaconda on your server then you will need to 
check which version of anaconda you have. This can be done by typing 

```console
conda info
```

There are multi-places you can get this information, but if you look at the base environment you should see something
similar to:

```console
/cm/shared/languages/Anaconda3-2019.03  (read only)

or

/mnt/storage/software/languages/anaconda/Anaconda3-2018.12  (read only)
```

What we want here is the version number (ie 2019.03 or 2018.12). Now that we know which version of conda we need to use
we can add it to our environment via the code snippet below. Keep in mind if your uses another version of conda, to 
change the ending to be your version number.

```console
module add languages/python-anaconda3-2019.03
```

Once you have added conda to your environment, we now can create a virtual environment for python. To do this we give
the environment a name, and then set the python version (we want version 3) and then a list of modules. At the moment,
some modules are forked or otherwise not finalised so this package is not yet on pypi so just install pip, git, and 
setuptools in the base environment.

```console
conda create --name pyGenicPipe python=3 pip git setuptools
```

You will be asked if you want to Process ([y]/n)? Type y and then enter, then it should build a base environment which
we can then activate via 

```console
conda activate pyGenicPipe
```

Now we are working in this virtual environment we can install our other packages from github itself. First we want to
install the base package. This will install all the modules we need to run the solution

```console
pip install git+https://github.com/sbaker-dev/pyGenicPipeline
```

However, we need to do another update to a forked version of bgen-reader as the current version will try to right 
metadata files to the directory the file is in, which is likely to lead to issues in write permissions. To fix this, 
just download forked version via

```console
pip install git+https://github.com/sbaker-dev/bgen-reader-py
```

This should mean you now have a valid setup to work with!
