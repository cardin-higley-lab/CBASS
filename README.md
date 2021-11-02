## Overview


![Python Versions](https://img.shields.io/badge/python-3.6%20%7C%203.7-blue)
![Matlab Versions](https://img.shields.io/badge/MATLAB-2018%7C2019%7C2020-blue.svg?style=flat-square)

__CBASS__ stands for _**C**lustering **B**and-limited **A**ctivity by **S**tate and **S**pectro-temporal features_. It is a method designed to detect recurring spatio-temporal motifs in multi-channel time series. Motifs must have energy in a specified frequency band and their occurence must increase during specific epochs of the recording (i.e. state). The method was originally developed to analyze multichannel Local Field Potentials but can be applied to any time series in principle. A detailed description of the algorithm can be found in the [wiki](https://github.com/cardin-higley-lab/CBASS/wiki)


## Code organization
CBASS is implemented in Matlab and Python. Each implementation is contained in a dedicated folder. The code is organized similarly for both and is subdived in 2 main folders:
1. ***Demo*** contains script showing exemples of how to run CBASS. These scripts are meant to be edited and modified to the convenience of the user.
2. ***Pipeline*** contains the core functions implementing CBASS as well as a number of utilities. This part of the code is meant to be copied and edited by users wishing to modify the method or reuse some of its part.

The matlab implementation has an additional folder ***Figures*** containing scripts generating the figures shown in the [wiki](https://github.com/cardin-higley-lab/CBASS/wiki). These script operate on test data that can be downloaded [here](https://osf.io/3k7a5/?view_only=bbcb6ac653d041fab0bd1618301cab30).

## Link to download example data
The example data can be downloaded from OSF [here](https://osf.io/3k7a5/?view_only=bbcb6ac653d041fab0bd1618301cab30)

## Requirements
### Matlab 
The method has been developed primarily in Matlab 2018b. It should work on any posterior version. Matlab can be downloaded [here](https://www.mathworks.com/products/matlab.html).
### Python 

## Getting started
### Matalb
1. Make a copy of the CBASS directory and it to your Matlab path   `addpath(genpath( path_to_CBASS ))`.
2. The function *CBASS_Main_DetectEvents* encapsulate the whole pipeline. A Demonstration of how to run *CBASS_Main_DetectEvents* can be found in the demo script *CBASS_Call_Main__Template.m* in the [Demo](https://github.com/cardin-higley-lab/CBASS/edit/master/matlab/Demo) folder. Make a copy of this script and modify it to load your data and reflect your local path.
3. The source code for *CBASS_Main_DetectEvents* can be found in the [Pipeline](https://github.com/cardin-higley-lab/CBASS/edit/master/matlab/Pipeline) folder. The help section gives a detailed description its input, output and of the sufunctions implementing the different levels of the pipeline. Subfunctions also have a detailed help. Read and hack as needed.
### Python

## Setting up conda environment for optional plot
### Replicating our environment
Generating [PHATE](https://github.com/KrishnaswamyLab/PHATE) or [UMAP](https://umap-learn.readthedocs.io/en/latest/) plots of the intermediary steps of CBASS requires to set up [miniconda](https://docs.conda.io/en/latest/miniconda.html) together with python 3.7.6. Once miniconda is installed one need to set up an environment with the name gammaBouts_env. This can be done with the following steps:
1. The conda environment definition is in [gammaBouts_env.yml](gammaBouts_env.yml). To replicate this environment, do the following:
```
conda env create --name gammaBouts_env --file gammaBouts_env.yml
```
Further information about replicating the environment can be found at the [conda website](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually). 

2. Activate the environment:
+ Windows: `conda activate gammaBouts_env`
+ Linux/Mac: `source activate gammaBouts_env`

### Create your own environment
Replicating a conda environment can be challeinging. To work around these difficulties, a minimally functional conda environment can be created with the following steps:
1. Create an enviroment with python 3.7
```
conda create --name gammaBouts_env python=3.7.6
```

2. Activate your gammaBouts_env environment 
```
conda activate gammaBouts_env
```

3. Install key packages using the follwing commands
```
conda install -c conda-forge umap-learn
conda install -c anaconda matplotlib
conda install -c bioconda scprep
conda install -c bioconda phate
```

4. Needed for python version of CBASS
```
conda install h5py
conda install -c anaconda networkx
conda install pywavelets
conda install -c conda-forge statsmodels
conda install -c conda-forge python-louvain
```
