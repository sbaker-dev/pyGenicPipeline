---
layout: default
title: Setup
nav_order: 2
has_children: true
---

## Setup

### Server

If you are using this pipeline on a sever you will need to setup a virtual environment for python. If you haven't done
so before, check out the [conda setup](conda.md) instructions. These also include the instructions for setting up the 
python package as well.

### Python

To install the base package we just need to type the following command into the python console

```console
pip install git+https://github.com/sbaker-dev/pyGenicPipeline
```

However, we need to do another update to a forked version of bgen-reader as the current version will try to right 
metadata files to the directory the file is in, which is likely to lead to issues in write permissions. To fix this, 
just download forked version via

```console
pip install git+https://github.com/sbaker-dev/bgen-reader-py
```

### GUI

To install the GUI go the releases part of the [github page][GHP] and download the most up to date release. Once you 
have download this .... [TODO write this once complied]


[GHP]: https://github.com/sbaker-dev/pyGenicPipeline
