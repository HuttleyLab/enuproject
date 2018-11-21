import os, re, pickle, time, glob

from cogent3 import LoadSeqs, LoadTable, DNA
from cogent3.core.annotation import Feature
from ensembldb3 import Genome, HostAccount

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

def var_match_exon(exon_strand, var_strand, alleles, flank_5, flank_3):
    """return the matched variant data on exon strand, including snp base, 
    alleles, and flanking sequences.
    if the variant has the same strand as the allocated exon, return the 
    original variant data;
    if the variant has different strand to the allocated exon, return the 
    complement data"""
    #annotate the variant with its var_symbol and allele bases at the coordinate position in an exon sequence  
    if var_strand == exon_strand:
        new_alleles = alleles
        new_flank_5 = flank_5
        new_flank_3 = flank_3

    #get the complement allele bases if the variant strand and exon strand are different.
    else:
        orig_alleles = DNA.make_seq(alleles.split('/'))
        rc_alleles = '/'.join(orig_alleles.complement())
        new_alleles = rc_alleles
        new_flank_5 = flank_3.rc()
        new_flank_3 = flank_5.rc()
    
    return new_alleles, new_flank_5, new_flank_3
  
def get_coord(exon_strand, exon_start, exon_end, var_start, var_end):
    """return variant coordinates according to the exonic strand"""
    if exon_strand == 1:
        var_coord = (int(var_start) - int(exon_start), int(var_end) - int(exon_start))
        
    if exon_strand == -1:
        var_coord = (int(exon_end) - int(var_end), int(exon_end) - int(var_start))
        
    return var_coord

def is_effect(required_effects, var_effects):
    #for snps with more than one effect types
    if type(var_effects) is list:
        if any(effect in required_effects for effect in var_effects) is False:
            print ("variant is not a synsynonymous or missense variant.")
            return None
        
        if 'splice_region_variant' in var_effects:
            print("The splice region variant is ignored")
            return None
    
    #for snps with only one effect type        
    if type(var_effects) is str:
        if var_effects not in required_effects:
            print ("variant is not a synsynonymous or missense variant.")
            return None
    
    return True

def get_var_info(gene, gene_id, flank_size):
    """given a gene, return translated exon sequences with variants positions annotated"""
    
    cds = gene.canonical_transcript.cds
    translated_exons = gene.canonical_transcript.translated_exons
    
    seq_records = DNA.make_seq('')
    for exon in translated_exons:
        exon_seq = exon.seq
        
        if exon.variants == ():
            seq_records += exon_seq
        
        else:
            exon_chrom = exon.location.coord_name
            exon_start = exon.location.start
            exon_end = exon.location.end
            exon_strand = exon.location.strand

            variants = exon.variants
            for var in variants:
                var_name = var.symbol
                var_alleles = var.alleles
                
                print ('Recording variant %s' % var_name)
                
                if var.somatic is True:
                    print("Only include germline variants, but %s is a somatic variant"%var_name)
                    continue
                
                if var.num_alleles != 2:
                    print("The number of SNP alleles should be 2, but variant %s has alleles %s"%(var_name, var_alleles))
                    continue
                
                var_effects = var.effect
                required_effects = ['synonymous_variant', 'missense_variant']
                if not is_effect(required_effects, var_effects):
                    continue
                
                var_strand = var.location.strand
                var_chrom = var.location.coord_name
                var_start = var.location.start
                var_end = var.location.end
                
                if var_start == var_end:
                    print("The variant %s is not allocated in side an exon"%var_name)
                    continue
                
                flank_seqs = var.flanking_seq
                if not flank_seqs:
                    continue
                
                flank_5 = flank_seqs[0][-int(flank_size):]
                flank_3 = flank_seqs[1][:int(flank_size)]
                adj_alleles, adj_flank_5, adj_flank_3 = var_match_exon(exon_strand, var_strand, var_alleles, flank_5, flank_3)
                var_coord = get_coord(exon_strand, exon_start, exon_end, var_start, var_end)

                #add annotation to each translated exon sequence 
                exon_seq.add_annotation(Feature, 'variant', '%s*%s*%s*%s*%s*%s*%s'%(var_name, var_chrom, exon_strand, var_effects, adj_alleles, adj_flank_5, adj_flank_3), [var_coord])
                
            seq_records += exon_seq
    
    assert (cds == seq_records),'for gene %s, the translated exon and CDs are different.'%gene_id
    
    ###put all annotations together and acquire variants' relative location on the whole CDs.
    var_records = seq_records.get_annotations_matching('variant')
    
    return var_records

def var_records(variant):
    """for each of the variant in the variant records, extract variant ID, variant location, variant alleles, variant effects, and the according location of the variant in the proteing coding sequence"""
    annotation = variant.split('"')[1]
    
    var_name = annotation.split('*')[0]
    var_chrom = annotation.split('*')[1]
    exon_strand = annotation.split('*')[2]
    var_effects = annotation.split('*')[3]
    var_alleles = annotation.split('*')[4]
    flank_5_seq = annotation.split('*')[5]
    flank_3_seq = annotation.split('*')[6]
    
    var_coord = re.compile('[at,\s/]+').split(str(variant).split('"')[2])[1].strip('[]').replace(":", ",")
    
    return (var_name, var_chrom, exon_strand, var_effects, var_alleles, flank_5_seq, flank_3_seq, var_coord)

@click.command()
@click.option('-i','--input_directory', help='Input directory containing one-to-one ortholog alignments.')
@click.option('-o','--output', help='Path to write variant information.')
@click.option('-f','--flank_size', help='Length of flank sequence from each side of a variant.')

def main(input_directory, output, flank_size):
    args = locals()
    
    if not os.path.exists(output):
        os.makedirs(output)
        
    logfile_path = os.path.join(output, "mouse_germline.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    
    account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
    mouse = Genome('mouse', release=88, account=account, pool_recycle=10000)
    
    input_path = os.path.abspath(input_directory)
    file_paths = get_files(input_path)
    
    for fn in file_paths:
        LOGGER.input_file(fn)
        
        gene_id = os.path.basename(fn).split('.')[0]
        
        print ("Acquiring variants from gene %s"%gene_id)
        gene = mouse.get_gene_by_stableid(stableid=gene_id)
        
        
        output_file = os.path.join(output, '%s.txt' % gene_id)
        
        start_time = time.time()
        num = 0
        with open(output_file, mode='w') as out_file:
            LOGGER.output_file(output_file)
            try:
                variants = get_var_info(gene, gene_id, flank_size)
                
                for var in variants:
                    record = var_records(str(var))
                    out_file.write('\t'.join(record) + '\n')
                    num += 1
                    
            except AssertionError:
                print('for gene %s, the translated exon and CDs are different.' % gene_id)
                os.remove(output_file)
            
            print ("finish getting variants on gene %s" % gene_id)
            LOGGER.log_message("%s"%num, label="Number of SNPs recorded")
            print ()
            
    print ()
    print ('Done')
    
    #determine runtime
    duration = time.time() - start_time
    LOGGER.log_message("%.2f" % (duration / 60.), label="run duration (minutes)")
    
   
if __name__ == "__main__":
    main()
