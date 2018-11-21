#######################
Mutation Classification
#######################


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
