#######################
Mutation Classification
#######################

************
Installation
************


Before installing mutation_classifier, please make sure the latest version of conda is installed.

Use the provided conda environment file as

```bash
$ conda-env create -f mutation_classifier.yml
```


*************
Data sampling
*************


Sampling the germline data
==========================

``read_mouse_germline.py`` reads files containing one-to-one ortholog alignments of mouse-rat-squirrel protein coding sequences, query the Ensembl variation database, and produce a summary table containing SNP symbols, locations, strands, effects, alleles, flanking sequences of required length (both 5' and 3', 250bp from each side), and relative coordinates of a SNP on mouse protein coding sequence.

The one-to-one ortholog alignment is obtained via the following steps:

1. Sampling homolog sequences by using Pycogent3 HomologSampler, details please refer to <https://bitbucket.org/pycogent3/homologsampler>
2. Get one-to-one alignment by using Phyg-align, details please refer to <https://bitbucket.org/gavin.huttley/phyg>

``sample_mouse_germline.py`` producing a summary table containing SNP ID, chromosome location, SNP strand, effects, alleles, ancestral base, variant base, flanking sequences of required lengths (250bp here) for mouse germline mutations, GC% and the response class, which is -1 for all germline muations. The results of this were saved into a .TXT file for later use.


Sampling the ENU data
=====================

Download ENU mutation files `SNVs20151101.txt <https://databases.apf.edu.au/mutations/>`_ from Australian Phenomics Facility database.

``get_ENU_variants.py`` reads ENU mutation files, according to chromosomes and coordinates given in the file, query Ensembl, obtain flanking sequences of required lengths (250bp), GC% and the response class, which is 1 for all ENU-induced muations. Then generate the data format consistent with that produced for the germline mutations.

.. ``sort_mut_dir.py`` categorise ENU and germline variant data according to their mutation directions, and save into different files.

**Currently, we are working on variant data on chromosome 1 only, I have uploaded both germline and ENU chrom1 data onto this repo, and also included 2 example data file. For the code testing purpose, we can just use these data, and don't need to run the following command to download all these data again**:blush:

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


********************************
Sample training and testing data
********************************

``sample_train_test.py`` will merge ENU and germline data together, and producing one hot endoded training samples and testing samples with defined sample size. Each sample contains response (1 indicates an ENU SNPs, and -1 indicates a germline SNPs), ``one-hot`` encoded mutation direction features, and ``one-hot`` encoded position combination features of required dimension.

Example workflow
================

This example command is wrote to do the following:

1. To split the entire data file (i.e. ENU_variants + germline_variants) into trainging data (50%) and testing data (50%). In this example, we are using example data from ``../data/germline_variants/example_germline_chrom1.txt`` and ``../data/ENU_variants/example_enu_chrom1.txt``.

2. Encode the resulted training and testing data in ``one-hot`` scheme respectively, consider a 7bp window size (i.e. 3 flanking bases on each side of the mutation), and up to 1-way interactions between neighbouring position.

3. Repeat the whole sampling process 5 times so 5 different sample-sets will be produced and saved into the output directory.

are as follows:
::

	$ python scripts/sample_train_test.py -e data/ENU_variants/example_enu_chrom1 
	-g data/germline_variants/example_germline_chrom1 
	-o dir/path/to/save/output --train_size=0.5 --test_size=0.5 
	--directions All -d 1 -t 5 -f 3

After implementing ``sample_train_test.py``, the following files will be saved into specified output directory:

- ``training_i.csv``
- ``testing_i.csv``
- ``log/`` directory  

**********************************
Logistic regression classification
**********************************

``classification_analysis.py`` will implement the classification analysis, and produce the performance AUC score. When doing classification, three separate analyses will be done:


**1. Setting C and penalty value with cross validation**

For the logistic regression classification, hyperparameter C and penalty needs to be set. To do this, we split training data to actual training data and validation data. We train the classifier on actual training data, and set hyperparameters on validation data. Within each validation process, performances of algorithms with different C and penalty values were compared, and the hyperparameter generating the best performance was saved for further analyses. The default C options are set as '0.1,1,10,100', and the default penalty option is set as l1.


**2. Evaluating classification performance on training data**

After obtaining the classifier, we evaluate the performance of the classifier on training data.


**3. Evaluating classification performance on testing data**

Finally, we evaluate the performance of the classifier on testing data.


Example workflow
================
1. Find the best classifier by doing cross validation:
::

	$ python scripts/classification_analysis.py logreg_train_and_validate 
		--training_data path/to/training_i.csv -o output/dir -f 3 -gc

Remove ``-gc`` flag if GC% feature is not included in the analyses.

After implementing this command, the following files will be saved in to the specified output directory:

- ``cross_validation/cv_training.csv`` 
- ``cross_validation/cv_validation.csv``  
- ``logreg_classifier.pkl``
- ``cv_report.json``
- ``logs/`` directory
- ``betas.json``
- ``scaler.pkl`` if ``-gc`` flag is ON in the command line


2. Apply resulted classifier on training data:
::

	$ python scripts/classification_analysis.py logreg_on_train -o output/dir -gc 

Remove ``-gc`` flag if GC% feature is not included in the analyses.

After implementing this command, a ``clf_on_train.json`` file is saved into specified directory.

**Very important!!!** The ``--output_dir`` option defined for step 1 and step 2 should be exactly the same, because in step 2, we are applying the classifier (and scaler) obtained from the 1st step to the exact cv_training file produced from the first step, and then save results into the same output directory, therefore, please do not change the ``--output_dir`` option setting in this step.

3. Apply resulted classifier on testing data:
::

	$ python scripts/classification_analysis.py logreg_testing 
		--testing_data path/to/testing_i.csv -o output/dir 
		--clf_file path/to/classifier.pkl 
		--scaler_file path/to/scaler.pkl -gc

Remove ``-gc`` flag and ``scaler_file`` if GC% feature is not included in the analyses.

Here the output directory can be difined as any output directory.

After implementing this command, a ``classification_report.json`` file is saved into specified directory.


************************************
Bernoulli Naive Bayes classification
************************************

The overall NB classification analyses is very similar to the logistic regression classification analyses, it also contains the three-step analyses: setting hyperparameter with cross validation, evaluating classifier on training data, and evaluating classifier on testing data. The same training data files, and sample testing data files as used to do the logistic regression classification analysis are used here.


**1. Setting alpha value with cross validation**

For the Naive Bayes classification, hyperparameter alpha needs to be set. To do this, we split training data to actual training data and validation data. We train the classifier on actual training data, and set hyperparameters on validation data. Within each validation process, performances of algorithms with different alpha values were compared, and the alpha value generating the best performance was saved for further analyses. The default alpha options are set as '0.01,0.1,1,2,3'.


**2. Evaluating classification performance on training data**

After obtaining the classifier, we evaluate the performance of the classifier on training data.


**3. Evaluating classification performance on testing data**

Finally, we evaluate the performance of the classifier on testing data.


Example workflow
================
1. Find the best classifier by doing cross validation:
::

	$ python scripts/classification_analysis.py bernoullinb_train_and_validate 
	--training_data path/to/training_i.csv -o output/dir -gc

Remove ``-gc`` flag if GC% feature is not included in the analyses.

After implementing this command, the following files will be saved in to the specified output directory:

- ``cross_validation/cv_training.csv`` 
- ``cross_validation/cv_validation.csv``  
- ``nb_classifier.pkl``
- ``cv_report.json``
- ``logs/`` directory
- ``scaler.pkl`` if ``-gc`` flag is ON in the command line


2. Apply resulted classifier on training data:
::

	$ python scripts/classification_analysis.py bernoullinb_on_train 
		-o output/dir -gc

Remove ``-gc`` flag if GC% feature is not included in the analyses.

After implementing this command, a ``clf_on_train.json`` file is saved into specified directory.

**Very important!!!** The ``--output_dir`` option defined for step 1 and step 2 should be exactly the same, because in step 2, we are applying the classifier (and scaler) obtained from the 1st step to the exact cv_training file produced from the first step, and then save results into the same output directory, therefore, please do not change the ``--output_dir`` option setting in this step.

3. Apply resulted classifier on testing data:
::
	
	$ python scripts/classification_analysis.py bernoullinb_testing 
		--testing_data path/to/testing_i.csv -o output/dir 
		--clf_file path/to/classifier.pkl 
		--scaler_file path/to/scaler.pkl -gc

Remove ``-gc`` flag and ``--scaler_file`` if GC% feature is not included in the analyses.

After implementing this command, a ``classification_report.json`` file is saved into specified directory.
