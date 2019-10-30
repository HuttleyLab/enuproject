#!/usr/bin/env python3.5
"""code for annotating SNP positions on reference sequences in an alignment"""
import os, re, pickle, glob, gzip, time

import copy
from copy import deepcopy
from cogent3 import LoadSeqs, LoadTable, LoadTree, DNA
from cogent3.core.annotation import Feature
from cogent3.evolve.models import HKY85
from cogent3.core.alignment import Alignment

import click
from scitrack import CachingLogger

LOGGER = CachingLogger()

def get_files(input_path):
    """return a list of files with a specific glob pattern from a directory"""
    fns = []
    for file_extension in ['fasta.gz', 'fasta', 'fa', 'fa.gz']:
        fns += glob.glob(os.path.join(input_path, '*.' + file_extension))
    if not fns:
        raise RuntimeError('Error selscting files')
    
    return fns

def get_syntenic_alignment(aln, var_name, var_start, ref_name, aln_flank):
    """target the variant and slice out the syntenic region, with the
     variant allocate in the middle."""
    annot_synt_region = aln.get_seq(ref_name).add_annotation(Feature, var_name, 'syntenic_region', [(max(var_start - aln_flank, 0), var_start + 1 + aln_flank)])
    annotate_aln = list(aln.get_annotations_from_seq(ref_name, var_name))[0]
    sliced_aln = annotate_aln.get_slice()
    
    return sliced_aln
    
def align_checker(syn_aln, ref_name, aln_flank, min_length):
    """returns alignment with annotated SNP at central position of ref_name
    
    Returns None if:
     - the ref sequence length does not satisfy 2x+1,
     - the number of strictly nucleotide columns is < min_length
     - non ref_name sequences have missing data at the SNP position
     """
    ref_seq = syn_aln.get_seq(ref_name)
    ref_length = len(ref_seq)
    if ref_length != aln_flank * 2 + 1:
        print("The reference sequence does not satisfy the length requirement")
        return None
    
    no_gaps = syn_aln.no_degenerates()
    if no_gaps is None:
        print("The number of strictly nucleotide columns is less than the minimum requirement.")
        return None
        
    if no_gaps is not None:
        if len(no_gaps) <= min_length:
            print("The number of strictly nucleotide columns is less than the minimum requirement.")
            return None
    
    variant = syn_aln.get_seq(ref_name).add_feature('variant', 'variant', 
                           [(aln_flank, aln_flank + 1)])
    aln_vars = syn_aln.get_annotations_from_seq(ref_name, 'variant')
    
    aln_var_bases = list(aln_vars)[0].get_slice().todict()
    if not set(aln_var_bases.values()) <= set(DNA):
        print("The variant base must be a DNA nucleotide.")
        return None
    
    return syn_aln

def get_root_alns(aln, ref_name, aln_flank):
    """returns alns from selected files, build a tree and produce ancestor 
    sequence and then project a sequence position onto the alignment position,
    and get the ancestral state, which is the mutation start state"""
    ref_seq = aln.get_seq(ref_name)
    variant = aln.get_seq(ref_name).add_feature('variant', 'variant', 
                           [(aln_flank, aln_flank + 1)])
    
    tree = LoadTree(tip_names=aln.names)
    sm = HKY85()
    lf = sm.make_likelihood_function(tree, digits=3, space=2)
    lf.set_alignment(aln)
    lf.optimise(show_progress=False, local=True)
    ancestor = lf.likely_ancestral_seqs()
    
    seq = ancestor.take_seqs(['root'])
    alns = aln.add_seqs(seq)
    
    return alns

def get_start_state(aln, ref_name, aln_flank):
    """load aln from selected files, build a tree and produce ancestor 
    sequence and then project a sequence position onto the alignment position,
    and get the ancestral state, which is the mutation start state
    """
    alns = get_root_alns(aln, ref_name, aln_flank)
    variant = alns.get_annotations_from_seq(ref_name, 'variant')
    aln_bases = list(variant)[0].get_slice()
    rt = aln_bases.take_seqs(['root'])
    rt_base = list(rt.todict().values())[0]
    
    return rt_base

def get_end_state(start_base, allele):
    """get the ending state of SNP of interest"""
    if start_base in allele:
        end_base = allele - set(start_base)
        end_base = ''.join(end_base)
        if len(end_base) is 1:
            return end_base
        else:
            return None
    else:
        return None

def is_correct_chrom(chroms, chromosome):
    if chroms is 'All':
        return lambda x: True
    if chroms is not 'All':
        chroms = chroms.split(',')

        return chromosome in [c.upper() for c in chroms]

def is_correct_effect(var_effect, effects):
    if var_effect is 'All':
        return lambda x: True
    if var_effect is not 'All':
        return var_effect in effects

def get_gc(seq):
    """to calculate the total GC content of a sequence"""
    gc = sum(seq.count(b) for b in 'GC')
    total = len(seq)
    return gc / total

def get_var_data(aln, variant, ref_name, aln_flank, min_length, chroms):
    if variant == ['']:
        return None
        
    [var_name, var_chrom, exon_strand, var_effects, var_alleles, flank_5_seq, flank_3_seq, var_coord] = variant
    
    if not is_correct_chrom(chroms, var_chrom):
        return None 
    
    var_alleles = set(var_alleles.split('/'))
    
    var_start = int(var_coord.split(',')[0])
    #get  alignment
    syn_aln = get_syntenic_alignment(aln, var_name, var_start, ref_name, aln_flank)
    syn_aln = LoadSeqs(data=copy.deepcopy(syn_aln.todict()), moltype=DNA, array_align=False)
    #check alignments and only keep the alignment meet requirements
    checked_aln = align_checker(syn_aln, ref_name, aln_flank, min_length)
    if not checked_aln:
        return None
    
    start_base = get_start_state(checked_aln, ref_name, aln_flank)
    
    end_base = get_end_state(start_base, var_alleles)
    
    if not end_base:
        return None
    
    if end_base is '':
        return None
    
    nbr_seq = DNA.make_seq(flank_5_seq + flank_3_seq)
    gc_content = get_gc(nbr_seq)
    
    allele_freqs = pep_alleles = gene_loc = gene_id = 'None'
    
    response = '-1'
    
    return (var_name, var_chrom, exon_strand, var_effects, allele_freqs, str(var_alleles), str(start_base), str(end_base), str(flank_5_seq), str(flank_3_seq), str(gc_content), pep_alleles, gene_loc, gene_id, response)

@click.command()
@click.option('-a','--aln_directory', help='Directory containing one-to-one ortholog alignments.')
@click.option('-v','--variant_records', help='Directory containing files of variant records')
@click.option('-o','--output_datafile', help='Path to write data.')
@click.option('-n','--ref_name', help='Reference species')
@click.option('-f','--aln_flank', type=int, help='Required flank length in an alignment on each side of a variant')
@click.option('-min_len','--min_length', type=int, help='Minimum number of strictly nucleotide columns in syntenic sub-alignment with variant allocating in the middle')
@click.option('-e','--var_effect', default='All',
    type=click.Choice(['All', 'missense_variant', 'synonymous_variant']), help='Variant effect')
@click.option('-c','--chroms', default='All', 
help="Chromosome location of interest. Examples of choices of the chroms option are:\
     - 'All' for all chromosomes, \
     - '1' for chromosome 1, or \
     - '1,2,3' for multiple chromosomes 1, 2 and 3, chromosomes are separated by comma")


def main(aln_directory, variant_records, output_datafile, aln_flank, min_length, ref_name, var_effect, chroms):
    args = locals()
    
    output_dir = os.path.dirname(output_datafile)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logfile_path = os.path.join(output_dir, "logs/sample_germline.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    
    input_path = os.path.abspath(aln_directory)
    aln_files = get_files(input_path)
    
    num = 0
    with open(output_datafile, mode='w') as output:
        LOGGER.output_file(output_datafile)
        start_time = time.time()
        
        for aln_f in aln_files:
            LOGGER.input_file(aln_f)
            gene_id = os.path.basename(aln_f).split('.')[0]
            
            print("get syntenic alignment for protein coding sequence in gene %s"%gene_id)
            
            var_record_file = os.path.join(variant_records, '%s.txt'%gene_id)
            if not os.path.isfile(var_record_file):
                print("Variant file for gene %s does not exist"%gene_id)
                continue
            
            aln = LoadSeqs(aln_f, moltype=DNA, array_align=False)
            ref_name = ref_name
            
            with open(var_record_file, 'r') as input_file:
                for line in input_file:
                    variant = line.strip().split('\t')
                    print('Recording variant %s'%variant[0])
                    
                    if is_correct_effect(var_effect, variant[4]):
                        record = get_var_data(aln, variant, ref_name, aln_flank, min_length, chroms)
                        
                        if record is None:
                            continue
                            
                        output.write('\t'.join(record)+'\n')
                        num += 1
    
    print ('DONE!')
    print ('Number of variants recorded: %s' % num)
    duration = time.time() - start_time
    LOGGER.log_message("%.2f" % (duration / 60.), label="run duration (minutes)")
    LOGGER.log_message("%s"%num, label="Recorded snp number")
                    
   
if __name__ == "__main__":
    main()
