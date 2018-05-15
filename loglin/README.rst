##########################################
Analysing mutations with log-linear models
##########################################

This loglin repository contains scripts used to generate the samples and perform the analyses reported in *Statistical methods for identifying sequence motifs affecting point mutations* by Zhu, Neeman, Yap and Huttley. Running all these scripts requires the `MutationMotif Python library <https://bitbucket.org/gavin.huttley/mutationmotif>`_, and all of its dependencies. In addition, an install of Jupyter notebook is required. (Note, scripts ending in ``.ipy`` and ``.ipynb`` require Jupyter/Ipython.)

As the repository also includes the counts data from which all analyses were conducted, running all scripts listed under the "Analysis of neighbourhood effects" and the "Analysis of mutation spectra" sections should reproduce exactly the tables and figures reported in the manuscript.

***********************************
Sampling the mouse spontaneous data
***********************************

``read_mouse_germline.py`` reads files containing one-to-one ortholog alignments of mouse-rat-squirrel protein coding sequences, query the Ensembl variation database, and produce a summary table containing SNP symbols, locations, strands, effects, alleles, flanking sequences of required length (both 5' and 3', 250bp from each side), and relative coordinates of a SNP on mouse protein coding sequence.

The one-to-one ortholog alignment is obtained via the following steps:

1. Sampling homolog sequences by using Pycogent3 HomologSampler, details please refer to `HomologSampler Bitbucket page <https://bitbucket.org/pycogent3/homologsampler>`_.
2. Get one-to-one alignment by using Phyg-align, details please refer to `Phyg-align Bitbucket page <https://bitbucket.org/gavin.huttley/phyg>`_.

``sample_mouse_germline.py`` producing a summary table containing SNP ID, chromosome location, SNP strand, effects, alleles, ancestral base, variant base and flanking sequences of required lengths (250bp here) for mouse germline mutations. The results of this were saved into a .TXT file for later use.


*********************
Sampling the ENU data
*********************

Download ENU mutation files `SNVs20151101.txt <https://databases.apf.edu.au/mutations/>`_ from Australian Phenomics Facility database.

``get_ENU_variants.py`` reads ENU mutation files, according to chromosomes and coordinates given in the file, query Ensembl and obtain flanking sequences of required lengths (250bp). Then generate the data format consistent with that produced for the germline mutations.


**************************************
Convert sequence data into count table
**************************************

``generate_alignments.ipy`` calls the `snptables_to_aln.py`, and aligns sequences (sampled from previous steps) together with mutations locating in the middle, and save the alignment into FASTA format.

``generate_counts.ipy`` calls the ``aln_to_counts`` command in `MutationMotif Python <https://bitbucket.org/gavin.huttley/mutationmotif>`_, and generate a count table for neighbouring base counts with flank size of 2 bp on each side of a mutation.

``generate_long_flanks_counts.ipy`` calls the ``aln_to_counts`` command in `MutationMotif <https://bitbucket.org/gavin.huttley/mutationmotif>`_, and generate a count table for neighbouring base counts with longer flank size of 10 bp on each side of a mutation.


*********************************
Analysis of neighbourhood effects
*********************************

``do_nbr.ipy`` Statistical analysis of a single group for contributions of neighbouring bases to point mutations. Results from this written to
``../results/<mut_source>/<chrom>/directions/<direction>``. 

``do_nbr_long_flanks.ipy`` Statistical analysis of a single group for contributions of larger neighbouring bases to point mutations. Results written to ``../results/<mut_source>/long_flanks/<chrom>/directions/<direction>``

``do_nbr_compare.ipy`` Statistical analysis comparing neighbourhood effects between groups (e.g. between ENU-induced mutations and spontaneous germline mutations). Results are written to ``<group_vs_group>`` directories within ``../results/``, depending on the level of comparison, e.g. ``../results/ENU_vs_germline``.


****************************
Analysis of mutation spectra
****************************

Using the ``spectra`` command in `MutationMotif <https://bitbucket.org/gavin.huttley/mutationmotif>`_ to perform the spectra comparison between groups (e.g. between ENU-induced mutations and spontaneous germline mutations). Output files (``spectra_analysis.json``, ``spectra_analysis.log`` and ``spectra_summary.txt``) are written to the same directories as the corresponding neighbour analysis.





