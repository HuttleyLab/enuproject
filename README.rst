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

*************
Data sampling
*************

Sampling the germline data
==========================

``read_mouse_germline.py`` reads files containing one-to-one ortholog alignments of mouse-rat-squirrel protein coding sequences, query the Ensembl variation database, and produce a summary table containing SNP symbols, locations, strands, effects, alleles, flanking sequences of required length (both 5' and 3', 250bp from each side), and relative coordinates of a SNP on mouse protein coding sequence.

The one-to-one ortholog alignment is obtained via the following steps:

1. Sampling homolog sequences by using Pycogent3 HomologSampler, details please refer to `HomologSampler Bitbucket page <https://bitbucket.org/pycogent3/homologsampler>`_.
2. Get one-to-one alignment by using Phyg-align, details please refer to `Phyg-align Bitbucket page <https://bitbucket.org/gavin.huttley/phyg>`_.

``sample_mouse_germline.py`` producing a summary table containing SNP ID, chromosome location, SNP strand, effects, alleles, ancestral base, variant base, flanking sequences of required lengths (250bp here) for mouse germline mutations, GC% and the response class, which is -1 for all germline muations. The results of this were saved into a .TXT file for later use.


Sampling the ENU data
=====================

Download ENU mutation files `SNVs20151101.txt <https://databases.apf.edu.au/mutations/>`_ from Australian Phenomics Facility database.

``get_ENU_variants.py`` reads ENU mutation files, according to chromosomes and coordinates given in the file, query Ensembl, obtain flanking sequences of required lengths (250bp), GC% and the response class, which is 1 for all ENU-induced muations. Then generate the data format consistent with that produced for the germline mutations.

.. ``sort_mut_dir.py`` categorise ENU and germline variant data according to their mutation directions, and save into different files.

Example workflow
================

In this example, we are working on SNP data on chromosome 1 only. 

1. Acquire mouse-rat-squirrel protein coding sequences. 
::

	$ homolog_sampler one2one --ref=Mouse --species=Mouse,Rat,Squirrel --release=88 
		--outdir=mrs_homolog_seqs/ --force_overwrite

2. Produce CDS one-to-one alignments. 
::
	
$ phyg-align -i mrs_homolog_seqs/ -o mrs_homolog_alns/ -m HKY85-GAMMA

3. Produce a summary table containing raw mouse germline SNP data
::
	
$ python read_mouse_germline.py -i ../path/to/mrs_homolog_alns/ -o ../path/to/mrs_homolog_vars/ -f 250

4. Produce a summary table containing required mouse germline SNP information.
::

	$ python sample_mouse_germline.py -a ../path/to/mrs_homolog_alns/ 
		-v ../path/to/mrs_homolog_vars/ 
		-o ../path/to/germline_variants/germline_chrom1.txt 
		-n Mouse -f 10 -min_len 5 -c 11

6. Produce a summary table containing required ENU-induced SNP information.
::

	$ python get_ENU_variants.py -i ../path/to/SNVs20151101.txt 
	-o ../path/to/ENU_variants/enu_chrom1.txt -f 250 -c 1


After implementing these commands, the resulted germline data and ENU data have a consistent format as follows:

+-------------+------------+--------+--------------------+--------------+-----------+----------+----------+-----------+----------+-------+-------------+----------+---------+----------+
| var_id      | chromosome | strand | effect             | allele freqs | alleles   | ref base | var base | 5' flank  | 3' flank | GC%   | pep_alleles | gene_loc | gene_id | response |
+-------------+------------+--------+--------------------+--------------+-----------+----------+----------+-----------+----------+-------+-------------+----------+---------+----------+
| germline_01 | 1          | -1     | synonymous_variant | None         | {'A','G'} | G        | A        | GAG...TGA | TC...GGT | 0.348 | None        | None     | None    | -1       |
+-------------+------------+--------+--------------------+--------------+-----------+----------+----------+-----------+----------+-------+-------------+----------+---------+----------+
| ...         | ...        | ...    | ...                | ...          | ...       | ...      | ...      | ...       | ...      | ...   | ...         | ...      | ...     |          |
+-------------+------------+--------+--------------------+--------------+-----------+----------+----------+-----------+----------+-------+-------------+----------+---------+----------+

+-------------+------------+--------+--------------------+--------------+-----------+----------+----------+-----------+----------+-------+-------------+----------+---------+----------+
| var_id      | chromosome | strand | effect             | allele freqs | alleles   | ref base | var base | 5' flank  | 3' flank | GC%   | pep_alleles | gene_loc | gene_id | response |
+-------------+------------+--------+--------------------+--------------+-----------+----------+----------+-----------+----------+-------+-------------+----------+---------+----------+
| ENU_01      | 1          | 1      | missense_variant   | None         | {'T','G'} | G        | T        | CAC...GTA | GC...GTG | 0.460 | None        | None     | None    | 1        |
+-------------+------------+--------+--------------------+--------------+-----------+----------+----------+-----------+----------+-------+-------------+----------+---------+----------+
| ...         | ...        | ...    | ...                | ...          | ...       | ...      | ...      | ...       | ...      | ...   | ...         | ...      | ...     |          |
+-------------+------------+--------+--------------------+--------------+-----------+----------+----------+-----------+----------+-------+-------------+----------+---------+----------+

