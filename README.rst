###########################
Classification of mutations
###########################

This repository contains scripts used to generate the samples and perform the analyses reported in Statistical methods for *Classifying ENU induced mutations from spontaneous germline mutations in mouse with machine learning techniques* by Zhu, Ong and Huttley. 

*************
Initial Setup
*************

``ENU_project.yml`` file assists with initial settings, it set up all required dependencies to run scripts in this repository.

Running the following commande to create a new environment named `mclass` using conda, and to install all required dependencies. (please make sure the latest version of conda is installed.)

`conda env create -f ENU_project.yml`

Activating the ``mclass`` environment by running:

`source activate mclass`

*************
repo content
*************
This repository contains scripts used to run two completely independent analyses: loglin analysis and mutation classification analysis. 

DATA is up on zenodo
data for the loglin
data for classifier
jupiter note book is used to reproduc the figure and table prudoced in the manuscript.






The latter basically logs commands, file inputs and outputs, to assist with reproducible research.

%manuscript contains two completely independent analysis

*******
log-lin
*******

classifier: using what we discovered to build clf for predicting the sample of individual point muts
the raw data is avaliable on zenodo, <\link>
we need a copy of 3-clause bsd lisence, there is a licence, these code is released under that licence.
This repo consists of two analyses components: 

*******
loglin
*******

This directory contains the scripts and results of the log-linear analysis for examination of the neighbourhood effects on ENU-induced and spontaneous germline mutations in the mouse, respectively. Please refer to the README.rst file inside this directory for detailed explaination regarding the package installation and the analyses implementation.

**********
classifier
**********

This directory contains the scripts and results of applying machine learning techniques to classify variants as having derived from ENU treatment or occur spontaneously in the mouse genome. Again, for detailed installation and implementation guidelines, please refer to the README.rst file inside this directory.





