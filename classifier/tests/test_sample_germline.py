import os, sys
sys.path.append('..')

from ensembldb3 import Genome, HostAccount

from cogent3 import LoadSeqs, DNA
from cogent3.util.unit_test import TestCase, main
from scripts.read_mouse_germline import var_match_exon, get_coord, is_effect, get_var_info, var_records
from scripts.sample_mouse_germline import get_syntenic_alignment, align_checker, get_root_alns, get_start_state, get_end_state, get_var_data


account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
mouse = Genome('mouse', release=88, account=account, pool_recycle=10000)

class TestSampleGermline(TestCase):
    def test_var_match_exon(self):
        """return the matched variant data on exon strand, including snp base, 
        alleles, and flanking sequences.
        if the variant has the same strand as the allocated exon, return the 
        original variant data;
        if the variant has different strand to the allocated exon, return the 
        complement data"""
        alleles = 'A/C'
        flank_5 = DNA.make_seq('AG')
        flank_3 = DNA.make_seq('TG')
        var_strand = 1
        
        exon_strand = 1
        alleles_got, flank_5_got, flank_3_got = var_match_exon(exon_strand, var_strand, alleles, flank_5, flank_3)
        alleles_exp, flank_5_exp, flank_3_exp = 'A/C', 'AG', 'TG'
        self.assertEqual(alleles_got, alleles_exp)
        self.assertEqual(flank_5_got, flank_5_exp)
        self.assertEqual(flank_3_got, flank_3_exp)
        
        exon_strand = -1
        alleles_got, flank_5_got, flank_3_got = var_match_exon(exon_strand, var_strand, alleles, flank_5, flank_3)
        alleles_exp, flank_5_exp, flank_3_exp = 'T/G', 'CA', 'CT'
        self.assertEqual(alleles_got, alleles_exp)
        self.assertEqual(flank_5_got, flank_5_exp)
        self.assertEqual(flank_3_got, flank_3_exp)
    
    def test_get_coord(self):
        """return variant coordinates according to the exonic strand"""
        exon_start = 20
        exon_end = 30
        var_start = 23 
        var_end = 24
        
        exon_strand = 1
        coord_got = get_coord(exon_strand, exon_start, exon_end, var_start, var_end)
        coord_exp = (3, 4)
        self.assertEqual(coord_got, coord_exp)
        
        exon_strand = -1
        coord_got = get_coord(exon_strand, exon_start, exon_end, var_start, var_end)
        coord_exp = (6, 7)
        self.assertEqual(coord_got, coord_exp)
    
    def test_is_effect(self):
        """return True if snp effect is synonymous or missense, otherwise, return False"""
        required_effects = ['synonymous_variant', 'missense_variant']
        self.assertTrue(is_effect(required_effects, var_effects='synonymous_variant'))
        self.assertTrue(is_effect(required_effects, var_effects='missense_variant'))
        self.assertIsNone(is_effect(required_effects, var_effects=['splice_region_variant', 'missense_variant']))
        self.assertTrue(is_effect(required_effects, var_effects=['other_variant', 'missense_variant']))
    
    def test_get_syntenic_alignment(self):
        """should return correct syntenic alignment"""
        aln = LoadSeqs(data=[('Rat', 'ATGTTCTTAGGACTCGTATCTCTTTTCTATTGCAGGATGAATAAGAGATACTTACAGAAAGCAACACAAGGGAAGCTTCTGATCATTATTTTTATAGTGACCTTGTGGGGGAAAGCTGTTTCCAGCGCCAACCATCACAAAGCTCACCATGTTAGAACTGGGACTTGCGAGGTGGTGGCGCTGCACAGATGCTGTAATAAGAACAAGATAGAAGAACGGTCCCAAACGGTCAAGTGCTCCTGCTTCCCTGGGCAGGTGGCAGGCACTACCCGAGCTGCCCCATCTTGTGTGGATGCATCCATAGTGGAACAAAAATGGTGGTGTCATATGCAGCCATGTCTGGAGGGAGAGGAATGTAAAGTCCTTCCAGATCGCAAAGGATGGAGCTGTTCCTCTGGAAACAAAGTAAAAACAACTAGGGTAACTCAT'),
                ('Mouse', 'ATG--------------ATCAC----------CAAGATGAATAAGAGATACTTGCAGAAAGCAACACAAGGAAAGCTTCTGATAATTATTTTTATAGTGACCTTGTGGGGGAAAGCCGTTTCCAGCGCCAACCATCACAAAGCTCACCATGTTAGAACTGGGACTTGCGAGGTTGTGGCGCTGCACAGATGCTGTAATAAGAACAAGATAGAAGAACGGTCCCAAACGGTCAAGTGCTCCTGCTTCCCTGGGCAGGTGGCAGGCACTACCCGAGCTGCTCCGTCTTGTGTGGATGCATCCATAGTGGAACAAAAGTGGTGGTGTCATATGCAGCCATGCCTGGAGGGAGAGGAATGTAAAGTCCTTCCAGATCGCAAAGGATGGAGCTGTTCCTCTGGAAACAAAGTAAAAACAACTAGGGTAACCCAT'),
                ('Squirrel', '---------------------------------------------------------------------------------------------------------------------------------------------GCTCACCATGTTAAAACGGGAACTTGCGAGGTGGTGGCACTCCACAGATGCTGTAATAAGAACAAGATAGAAGAACGATCACAGACAGTCAAGTGCTCCTGCTTCCCCGGGCAGGTGGCAGGCACCACACGAGCGGCTCCGTCTTGTGTGGATGCATCTATAGTGGAACAGAAATGGTGGTGTCATATGCAGCCATGTTTGGAGGGAGAGGAGTGTAAAGTTCTTCCAGATCGGAAAGGATGGAGCTGTTCTTCTGGGAATAAAGTGAAAACAACAAGGGTAAGTGGA')], moltype=DNA, array_align=False)
        var_name = 'rs265299216'
        var_start = 19
        ref_name = 'Mouse'
        flank_size = 20
        
        syn_aln_got = get_syntenic_alignment(aln, var_name, var_start, ref_name, flank_size)
        syn_aln_ex = LoadSeqs(
            data=[('Rat', 'ATGATCTCCAGGATGAATAAGAGATACTTACAGAAAGCAA'),
                ('Mouse', 'ATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA'),
             ('Squirrel', '----------------------------------------')],
        moltype=DNA, array_align=False)
        self.assertEqual(syn_aln_got, syn_aln_ex)
    
    def test_align_checker(self):
        """return filtered alignment"""
        ref_name = 'Mouse'
        flank_size = 20
        min_length = 10
        
        aln1 = LoadSeqs(
            data=[('Rat', 'ATGATCTCCAGGATGAATAAGAGATACTTACAGAAAGCAA'),
                ('Mouse', 'ATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA'),
             ('Squirrel', 'ATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA')],
        moltype=DNA, array_align=False)
        self.assertIsNone(align_checker(aln1, ref_name, flank_size, min_length))
        
        
        aln2 = LoadSeqs(
            data=[('Rat', 'AATGATCTCCAGGATGAATANGAGATACTTACAGAAAGCAA'),
                ('Mouse', 'AATGATCACCAAGATGAATANGAGATACTTGCAGAAAGCAA'),
             ('Squirrel', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA')],
        moltype=DNA, array_align=False)
        self.assertIsNone(align_checker(aln2, ref_name, flank_size, min_length))
        
        
        aln3 = LoadSeqs(
            data=[('Rat', 'AATGATCTCCAGGATGAATAAGAGATACTTACAGAAAGCAA'),
                ('Mouse', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA'),
             ('Squirrel', '-----------------------------------------')],
        moltype=DNA, array_align=False)
        self.assertIsNone(align_checker(aln3, ref_name, flank_size, min_length))
        
        aln4 = LoadSeqs(
            data=[('Rat', 'AATGATCTCCAGGATGAATAAGAGATACTTACAGAAAGCAA'),
                ('Mouse', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA'),
             ('Squirrel', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA')],
        moltype=DNA, array_align=False)
        
        syn_aln_ex = LoadSeqs(
            data=[('Rat', 'AATGATCTCCAGGATGAATAAGAGATACTTACAGAAAGCAA'),
                ('Mouse', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA'),
             ('Squirrel', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA')],
        moltype=DNA, array_align=False)
        self.assertEqual(aln4, syn_aln_ex)
    
    def test_get_root_alns(self):
        """return syntenic alignment and the inferring ancestral sequence"""
        aln = LoadSeqs(
            data=[('Rat', 'AATGATCTCCAGGATGAATAAGAGATACTTACAGAAAGCAA'),
                ('Mouse', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA'),
             ('Squirrel', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA')],
        moltype=DNA, array_align=False)
        ref_name = 'Mouse'
        r_aln_got = get_root_alns(aln, ref_name, aln_flank=20)
        r_aln_ex = LoadSeqs(
            data=[('Rat', 'AATGATCTCCAGGATGAATAAGAGATACTTACAGAAAGCAA'),
                ('Mouse', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA'),
             ('Squirrel', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA'),
                 ('root', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA')],
        moltype=DNA, array_align=False)
        
        self.assertEqual(r_aln_got, r_aln_ex)
    
    def test_get_start_state(self):
        """return the correct ancestral base for a variant"""
        aln = LoadSeqs(
            data=[('Rat', 'AATGATCTCCAGGATGAATAAGAGATACTTACAGAAAGCAA'),
                ('Mouse', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA'),
             ('Squirrel', 'AATGATCACCAAGATGAATAAGAGATACTTGCAGAAAGCAA')],
        moltype=DNA, array_align=False)
        ref_name = 'Mouse'
        
        r_base_got = get_start_state(aln, ref_name, aln_flank=20)
        r_base_ex = 'A'
        self.assertEqual(r_base_got, r_base_ex)
    
    def test_get_end_state(self):
        """return the correct ending base for a variant"""
        start_base = 'A'
        allele = set('G/A'.split('/'))
        
        e_base_got = get_end_state(start_base, allele)
        e_base_ex = 'G'
        self.assertEqual(e_base_got, e_base_ex)
        

if __name__ == '__main__':
    main()
