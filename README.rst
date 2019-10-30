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


`make_manuscript_figs.ipynb`, `make_manuscript_tables.ipynb` are Jupyter notebooks used to produce all figures and tables used in the manuscripts. Some compounded figures are saved into the `compound_figures` directory in LATEX format.

The BSD 3-clause license is included in this repo as well, refer to `license.txt`

*************
Data sampling
*************

Scripts are located in `sample_data/`.

Sampling the germline data
==========================

``read_mouse_germline.py`` reads files containing one-to-one ortholog alignments of mouse-rat-squirrel protein coding sequences, query the Ensembl variation database, and produce a summary table containing SNP symbols, locations, strands, effects, alleles, flanking sequences of required length (both 5' and 3', 250bp from each side), and relative coordinates of a SNP on mouse protein coding sequence.

The one-to-one ortholog alignment is obtained via the following steps:

1. Sampling homolog sequences by using HomologSampler, details please refer to `HomologSampler page <https://github.com/cogent3/homologsampler>`_.
2. Get one-to-one alignment by using Phyg-align. (This tool is defunct and has been replaced by ``cogent3``).

``sample_mouse_germline.py`` producing a summary table containing SNP ID, chromosome location, SNP strand, effects, alleles, ancestral base, variant base, flanking sequences of required lengths (250bp here) for mouse germline mutations, GC% and the response class, which is -1 for all germline muations. The results of this were saved into a .TXT file for later use.


Sampling the ENU data
=====================

Download ENU mutation files `SNVs20151101.txt <https://databases.apf.edu.au/mutations/>`_ from Australian Phenomics Facility database.

``get_ENU_variants.py`` reads ENU mutation files, according to chromosomes and coordinates given in the file, query Ensembl, obtain flanking sequences of required lengths (250bp), GC% and the response class, which is 1 for all ENU-induced muations. Then generate the data format consistent with that produced for the germline mutations.

Generating data for analysis
============================

In this example, we are working on SNP data on chromosome 1 only.

1. Acquire mouse-rat-squirrel protein coding sequences.
::

    $ homolog_sampler one2one --ref=Mouse --species=Mouse,Rat,Squirrel --release=88
        --outdir=mrs_homolog_seqs/ --force_overwrite

2. Produce CDS one-to-one alignments.
::

$ phyg-align -i mrs_homolog_seqs/ -o mrs_homolog_alns/ -m HKY85

.. note:: ``phyg-align`` is now replaced capabilities in a version of `cogent3`

3. Produce a summary table containing raw mouse germline SNP data::

    $ python read_mouse_germline.py -i ../path/to/mrs_homolog_alns/ -o ../path/to/mrs_homolog_vars/ -f 250

4. Produce a summary table containing required mouse germline SNP information. This infers the mutation direction from fitting a HKY85 model to the sequence alignments and inferring the most likely ancestral states for annotated SNOP locations.::

    $ python sample_mouse_germline.py -a ../path/to/mrs_homolog_alns/
        -v ../path/to/mrs_homolog_vars/
        -o ../path/to/germline_variants/germline_chrom1.txt
        -n Mouse -f 10 -min_len 5 -c 11

6. Produce a summary table containing required ENU-induced SNP information.::

    $ python get_ENU_variants.py -i ../path/to/SNVs20151101.txt
    -o ../path/to/ENU_variants/enu_chrom1.txt -f 250 -c 1

7. Run the ``variant_data/modify_variant_data.ipynb`` notebook to produce the standardised data formats for use in the analyses. Results are written to ``variant_data/Germline`` and ``variant_data/ENU`` as gzipped tsv files.

8. Run the ``variant_data/concat_data.ipynb`` notebook to merge the ENU and spontaneous germline varinats into single files for each chromosome. Results are written to ``variant_data/combined_data`` as gzipped tsv files. ENU variants are indicated by ``e`` value in the ``response`` column, while spontaneous germline are indicated by ``g``.

Using the ``mutation_origin`` command line tool
===============================================

We first note that ``mutation_origin`` is a rewrite of scripts authored by Yichneg Zhu. The rewrite was done to simplify inclusion of other classification algorithms. With hindsight of experience, optimisations for storage and performance were also included.

Note that all analyses done are logged using ``scitrack``. The generated log files are under the same directory and contain all run settings and md5 sums for the files used/produced.

For a full description of the command line options, see the ``mutation_origin`` `GitHub page <https://github.com/HuttleyLab/mutationorigin>`_.

Generating data for train and test
----------------------------------

::

    $ mutori_batch sample_data -ep variant_data/ENU/SNVs20151101_chrom1.tsv.gz -gp variant_data/Germline/mouse_germline_All_88_chrom1.tsv.gz -op classifier/chrom1_train/data -n 10 -N 3

Where ``-n`` is the number of replicates produced, ``-N`` the number of processors. This will generate balanced (equal numbers of randomly sampled ENU and Spontaneous germline) samples with total size of 1, 2, 4, 6, 8, and 16 thousand. The same samples are used for each classifier permutation.

Training classifiers, logistic regression as an example
-------------------------------------------------------

::

    $ mutori_batch lr_train -tp classifier/chrom1_train/data -op classifier/chrom1_train/lr/train -mr upto2 -N 20

This will trains a LR model with all possible terms up to 2-way interactions, for all data sets indicated by ``-tp`` and write the classifiers as python native serialised (``pickle`` formatted) files to matching paths indicated by ``-op``, using 20 processors.

Testing classifiers -- the prediction step
------------------------------------------

::

    $ mutori_batch predict -tp chrom1_train/data -cp chrom1_train/lr/train -op chrom1_train/lr/predict -N 3

Similar to above, it selects the matching files to those used for generating the classifier. For instance, for the classifier saved at ``chrom1_train/lr/train/1k/f0/train-0-classifier-lr.pkl`` will be applied to the testing data ``chrom1_train/data/1k/test-0.tsv.gz``. The result is a set of predictions for all the records in the testing set.

Evaluating performance
----------------------

::

    $ mutori_batch performance -tp chrom1_train/data -pp chrom1_train/lr/predict -op chrom1_train/lr/performance

Takes the results from the above and produces, for the performance statistics (typically AUC), the mean and standard deviation across cross-validation replicates.

Summarising performance across classifiers and sample sizes
-----------------------------------------------------------

::

    $ mutori_batch collate -bp chrom1_train/ -op chrom1_train/collated -O -ex genome

Takes all performance result files and combines into a single tsv. Excludes any files under the directory indicated by the ``-ex`` option. In this instance, this is where the whole genome prediction results are stored.

Predictions for the genome
==========================

Having chosen a classifier based on the last step, that classifier is applied to the entire genome, essentially recapping the steps from the prediction step through to the collate step. For example::

    $ mutori_batch predict -cp chrom1_train/lr/train/16k/f29d2p/train-1-classifier-lr.pkl -tp ../variant_data/combined_data/*.tsv.gz -op chrom1_train/genome/lr/predict -N 3

Where the value after ``-cp`` is the chosen LR classifier and ``-tp`` is the location of the genomic data.


