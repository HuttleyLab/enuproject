import os, re
import time, glob
from ensembldb3 import Genome, HostAccount

import click
from scitrack import CachingLogger

LOGGER = CachingLogger()

def get_files(input_dir):
    """return a list of .TXT files from a directory"""
    fns = []
    fns += glob.glob(os.path.join(input_dir, '*.txt'))
    
    if not fns:
        raise RuntimeError('Error selscting files')
    
    return fns

def is_correct_chrom(chroms, chromosome):
    if chroms is 'All':
        return lambda x: True
    if chroms is not 'All':
        chroms = chroms.split(',')

        return chromosome in [c.upper() for c in chroms]

def get_ENU_data(records, chroms, coord_range):
    """only sample data of inerested, 
    - variants with exonic effect only.
    - variants on a specific chromosome, 
    - or variants within a specific range in a genome."""
    var_id, chromosome, coordinate, ref_base, var_base, effect = \
    records[0], records[2], int(records[3]), records[15], records[16], records[19]
    
    if effect not in 'Non-synonymous, Synonymous':
        print  ("Not non-synonymous nor synonymous variant")
        return None

    if not is_correct_chrom(chroms, chromosome):
        return None
        
    if len(ref_base) == 0:
        return None
        
    if len(var_base) == 0:
        return None
    
    if not coord_range:
        #print ('there is no genemic location requirement')
        return (var_id, chromosome, int(coordinate), ref_base, var_base, effect)
        
    if coord_range:
        min_val, max_val = coord_range.split(',')
        min_val, max_val = int(min_val), int(max_val)
        
        if coordinate not in range(min_val, max_val + 1):
            return None
        else:
            return var_id, chromosome, int(coordinate), ref_base, var_base, effect

def get_flanks(snp_region, flank_size):
    """return the flanking regions for a location in a chromosome"""
    region_seqs = snp_region.seq
    rawflank5 = region_seqs[0: flank_size]
    rawflank3 = region_seqs[flank_size+1: flank_size*2+1]
    ensembl_ref = region_seqs[flank_size: flank_size+1]
    return rawflank5, ensembl_ref, rawflank3

def get_snp_region(chromosome, coordinate, flank_size, genome):
    """return snp region from Ensembl database, according to the "Chromosome" and "Coordinate (Assembly version: GRCm38)" data from the given ENU mutagenesis data table"""
    start = coordinate - flank_size
    end = coordinate + flank_size
    
    region = genome.get_region(coord_name=chromosome, start=start, end=end, ensembl_coord=True)
    
    return region

def get_snp_data(var_id, chromosome, coordinate, ref_base, var_base, effect, genome, flank_size):
    snp_region = get_snp_region(chromosome, coordinate, flank_size, genome)
    
    strand = str(snp_region.location.strand)
    
    f5, ensemble_base, f3 = get_flanks(snp_region, flank_size)
    
    if str(ensemble_base) != ref_base:
        print ("snp %s is not recorded, because the ensembl base %s is different from the given ref base %s."%(var_id, ensemble_base, ref_base))
        return None
    
    alleles = str(set((str(ensemble_base), var_base)))
    
    pep_alleles = gene_loc = gene_id = allele_freqs = 'None'
    
    return (var_id, chromosome, strand, effect, allele_freqs, alleles,
                str(ref_base), str(f5), str(f3), pep_alleles, gene_id, gene_loc)

@click.command()
@click.option('-i','--input_dir', help='Input file containing mutation data.')
@click.option('-o','--output_datafile', help='File to write output data.')
@click.option('-f','--flank_size', help='Number of flanking bases on each side of a variant.')
@click.option('-c','--chroms', default='All', 
help="Chromosome location of interest. Examples of choices of the chroms option are:\
     - 'All' for all chromosomes, \
     - '1' for chromosome 1, or \
     - '1,2,3' for multiple chromosomes 1, 2 and 3, chromosomes are separated by comma")
@click.option('-g','--coord_range', default = None, help='Genomic location range of interest. The minimum value and the maximum value to be seperated by comma.')

def main(input_dir, output_datafile, flank_size, chroms, coord_range):
    args = locals()
    
    output_dir = os.path.dirname(output_datafile)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logfile_path = os.path.join(output_dir, "logs/sample_ENU.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    
    account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
    mouse = Genome('mouse', release=88, account=account, pool_recycle=10000)
    
    input_dir = os.path.dirname(input_dir)
    file_paths = get_files(input_dir)
    
    start_time = time.time()
    num = 0 
    for fn in file_paths:
        with open(fn, mode='rb') as f_in:
            LOGGER.input_file(fn)
            with open(output_datafile, mode='w') as f_out:
                LOGGER.output_file(output_datafile)
                first_line = f_in.readline()
                for line in f_in:
                    line = line.decode("utf-8")
                    records = line.split('|')
                    if not get_ENU_data(records, chroms, coord_range):###
                        continue
                    
                    var_id, chromosome, coordinate, ref_base, var_base, effect = get_ENU_data(records, chroms, coord_range)###
                    
                    record = get_snp_data(var_id, chromosome, coordinate, ref_base, var_base, effect, mouse, int(flank_size))
                    
                    if not record:
                        continue
                    
                    f_out.write('\t'.join(record)+'\n')
                    num += 1
        
    LOGGER.log_message("%s" % num, label="number of variants written")
    f_out.close()
    
    #determine runtime
    duration = time.time() - start_time
    LOGGER.log_message("%.2f" % (duration / 60.), label="run duration (minutes)")        
        
if __name__ == "__main__":
    main()