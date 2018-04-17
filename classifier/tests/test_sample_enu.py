import os, sys
sys.path.append('..')

from ensembldb3 import Genome, HostAccount
from cogent3 import LoadTable, DNA
from cogent3.util.unit_test import TestCase, main
from scripts.get_ENU_variants import get_files, get_gc, is_correct_chrom, get_ENU_data,\
    get_flanks, get_snp_region, get_snp_data

account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
mouse = Genome('mouse', release=88, account=account, pool_recycle=10000)

class TestSampleENU(TestCase):
    def test_get_ENU_data(self):
        data = [['774776', 'IGL00164', '9', '119119678', 'T->P', '', '0.99', 'Probably damaging', '0.01', 'deleterious', 'MGI:2443671', 'Dlec1', 'deleted in lung and esophageal cancer 1 [Source:MGI Symbol;Acc:MGI:2443671]', 'Heterozygous', '', 'A', 'C', '617', '33.0', 'Non-synonymous', 'Cryopreserved', 'APF-G1'], 
        ['774777', 'IGL00164', '2', '87656238', 'M->V', '', '0.95', 'Possibly damaging', '0.02', 'deleterious', 'MGI:3030968', 'Olfr1134', 'olfactory receptor 1134 [Source:MGI Symbol;Acc:MGI:3030968]', 'Heterozygous', '', 'T', 'C', '325', '39.0', 'Splice', 'Cryopreserved', 'APF-G1'], 
        ['774779', 'IGL00164', '5', '31299576', 'V->M', '', '1.0', 'Probably damaging', '0.02', 'deleterious', 'MGI:1096345', 'Gckr', 'glucokinase regulatory protein [Source:MGI Symbol;Acc:MGI:1096345]', 'Heterozygous', 'Homozygous inactivation of this gene leads to reduced glucokinase protein levels and activity in the liver and altered glucose homeostasis.', 'G', 'A', '274', '39.0', 'Non-synonymous', 'Cryopreserved', 'APF-G1'], 
        ['774780', 'IGL00164', '9', '22073014', 'G->D', '', '0.12', 'Benign', '0.35', 'tolerated', 'MGI:1349469', 'Ecsit', 'ECSIT homolog (Drosophila) [Source:MGI Symbol;Acc:MGI:1349469]', 'Heterozygous', 'Homozygous mutant mice die around the stage of gastrulation showing abnormal epiblast patterning.', 'C', '', '231', '38.0', 'Non-synonymous', 'Cryopreserved', 'APF-G1'], 
        ['774781', 'IGL00164', '14', '55065026', 'S->P', '', '', 'Unknown', '', 'deleterious', 'MGI:2686934', 'Zfhx2', 'zinc finger homeobox 2 [Source:MGI Symbol;Acc:MGI:2686934]', 'Heterozygous', '', '', 'G', '185', '36.0', 'Non-synonymous', 'Cryopreserved', 'APF-G1']]
        
        self.assertIsNone(get_ENU_data(data[1], chroms='All', coord_range=None))
        self.assertIsNone(get_ENU_data(data[2], chroms='9', coord_range=None))
        self.assertIsNone(get_ENU_data(data[3], chroms='All', coord_range=None))
        self.assertIsNone(get_ENU_data(data[4], chroms='All', coord_range=None))
        
        var_id, chromosome, coordinate, ref_base, var_base, effect = get_ENU_data(data[0], chroms='9', coord_range=None)
        records_got = [var_id, chromosome, coordinate, ref_base, var_base, effect]
        records_expect = ['774776', '9', 119119678, 'A', 'C', 'Non-synonymous']
        self.assertEqual(records_got, records_expect)
        
    def test_get_files(self):
        """return a list of .TXT files from a directory"""
        input_dir = 'data/sample_enu/'
        file_paths_got = get_files(input_dir)
        file_paths_expect = ['data/sample_enu/ENU_test_sample_4.txt',
                             'data/sample_enu/ENU_test_sample_5.txt',
                             'data/sample_enu/ENU_test_sample_1.txt',
                             'data/sample_enu/ENU_test_sample_2.txt',
                             'data/sample_enu/ENU_test_sample_3.txt']
        self.assertEqual(file_paths_got, file_paths_expect)
    
    def test_get_gc(self):
        """return correct GC%"""
        seq1 = 'TATAGATGGGCATGTTGTCC'
        seq2 = 'CACAGGGATTCGGTAGCCGC'
        seq3 = 'AACCGCGCAGACGCAAAGAT'
        seq4 = 'AGGGCAGTGAGGAGGCAGGC'
        seq5 = 'TAGTTTAACAGCACATACAA'
        
        gc_got = [0.45, 0.65, 0.55, 0.7, 0.3]
        gc_expect = [get_gc(seq1), get_gc(seq2), get_gc(seq3), get_gc(seq4), get_gc(seq5)]
                
        self.assertEqual(gc_got, gc_expect)
                
    def test_is_correct_chrom(self):
        """return true if chromosome is in chroms"""
        chroms_1 = 'All'
        chroms_2 = '1'
        chroms_3 = 'x'
        chroms_4 = '1,4,x'
        chromosomes = ['1', '12', 'X']
        
        self.assertTrue(is_correct_chrom(chroms_1, chromosomes[0]))
        self.assertTrue(is_correct_chrom(chroms_2, chromosomes[0]))
        self.assertFalse(is_correct_chrom(chroms_3, chromosomes[0]))
        self.assertTrue(is_correct_chrom(chroms_4, chromosomes[0]))
        
        self.assertTrue(is_correct_chrom(chroms_1, chromosomes[1]))
        self.assertFalse(is_correct_chrom(chroms_2, chromosomes[1]))
        self.assertFalse(is_correct_chrom(chroms_3, chromosomes[1]))
        self.assertFalse(is_correct_chrom(chroms_4, chromosomes[1]))
        
        self.assertTrue(is_correct_chrom(chroms_1, chromosomes[2]))
        self.assertFalse(is_correct_chrom(chroms_2, chromosomes[2]))
        self.assertTrue(is_correct_chrom(chroms_3, chromosomes[2]))
        self.assertTrue(is_correct_chrom(chroms_4, chromosomes[2]))
    
    def test_get_flanks(self):
        """return the correct flanking regions for a location in a chromosome"""
        snp_region = get_snp_region(chromosome='9', coordinate=119119678, flank_size=2, genome=mouse)
        rawflank5_got, ensembl_ref_got, rawflank3_got = get_flanks(snp_region, flank_size=2)
        rawflank5_exp, ensembl_ref_exp, rawflank3_exp = 'CC', 'A', 'CC'
        self.assertEqual(rawflank5_got, rawflank5_exp)
        self.assertEqual(ensembl_ref_got, ensembl_ref_exp)
        self.assertEqual(rawflank3_got, rawflank3_exp)
    
    def test_get_snp_region(self):
        chromosome = '9'
        coordinate = 119119678
        flank_size = 2
        
        region = get_snp_region(chromosome, coordinate, flank_size, mouse)
        records_got = [region.genome.species, region.location.coord_name, region.location.start, region.location.end, len(region.seq)]
        records_expect = ['Mus musculus', chromosome, coordinate - flank_size - 1, coordinate + flank_size, 2 * flank_size + 1]
        self.assertEqual(records_got, records_expect)

if __name__ == '__main__':
    main()
