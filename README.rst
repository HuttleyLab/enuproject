###########################
Classification of mutations
###########################

This repository contains scripts used to generate the samples and perform the analyses reported in Statistical methods for *Classifying ENU induced mutations from spontaneous germline mutations in mouse with machine learning techniques* by Zhu, Ong and Huttley. 

*************
Initial Setup
*************

``ENU_project.yml`` file assists with setting up a virtual environment named ``mclass`` with all required dependencies installed within this virtual environment. Running the following commande to create the `mclass` environment using conda (please make sure the latest version of conda is installed).

```$ conda env create -f ENU_project.yml```

Activating the ``mclass`` environment by running:

```$ source activate mclass```

This environment has all required dependencies installed to run all scripts within this repository.


*************
Repo contents
*************

This repository contains scripts used to run two completely independent analyses. 

The `loglin` directory contains scripts and data to run the loglinear analysis. The log-linear analysis compares the neighbourhood effects between the ENU-induced and spontaneous mutations in mouse, in terms of the identity of the associated mutation motifs and their relative magnitude. Please refer to the README.rst file inside this directory for detailed explaination regarding the package installation and the analyses implementation.

The `classifier` directory contains scripts and data to perform the classification analysis. The classification analysis uses the neighbourhood effect discovered from the log-linear analysis to build classifiers, for predicting the mechanistic origin of individual point mutations. Again, for detailed installation and implementation guidelines, please refer to the README.rst file inside this directory.

The raw data used in this study are available at `Zenodo <http://zenodo.org/record/1204695>`_.


`manuscript_figs_tables.ipynb` is a Jupiter notenook used to produce all figures and tables used in the manuscripts. Some compounded figures are saved into the `compound_figures` directory in LATEX format.

The BSD 3-clause license is included in this repo as well, refer to `license.txt`
