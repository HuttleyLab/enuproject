import os, sys, numpy as np
sys.path.append('..')

from cogent3 import LoadTable, DNA
from cogent3.util.unit_test import TestCase, main
from scripts.sample_train_test import str_to_numerical, get_learning_data,\
 nck, get_nbr_combs_features, sort_col_by_dim, feature_selector, onehot_encoder

class TestTrainTestSample(TestCase):
    def test_str_to_numerical(self):
        """correctly convert number string to float or int type"""
        num_got = str_to_numerical(num_str='0.38')
        num_exp = float(0.38)
        self.assertEqual(num_got, num_exp)
        
        num_got = str_to_numerical(num_str='1')
        num_exp = int(1)
        self.assertEqual(num_got, num_exp)
        
        num_got = str_to_numerical(num_str='1.0')
        num_exp = float(1.0)
        self.assertEqual(num_got, num_exp)
    
    def test_get_learning_data(self):
        """read through a data file, correctly return var_id, response, 
        neighboring seqs and mutation direction """
        data_file = "data/ENU_data.txt"
        with open(data_file, mode='r') as d_file:
            data = [line.rstrip() for line in d_file]
            records_got = get_learning_data(data, flank_size=3)
            records_exp = [['774776', '1', '0.45', 'AtoC', 'GGGCAT'], 
                           ['774777', '1', '0.65', 'TtoC', 'ATTCGG'], 
                           ['774779', '1', '0.55', 'GtoA', 'CAGACG'], 
                           ['774780', '1', '0.7', 'CtoT', 'TGAGGA'], 
                           ['774781', '1', '0.3', 'AtoG', 'ACAGCA']]
            
            self.assertEqual(records_got, records_exp)
    
    def test_nck(self):
        """return the number of possible n choose k combinations, where n>k"""
        n_choice_got = nck(4, 2)
        n_choice_exp = 6
        
        self.assertEqual(n_choice_got, n_choice_exp)
    
    def test_get_nbr_combs_features(self):
        """read through a feature matrix row by row, transform neighboring sequences column in feature matrix into columns of neighboring base combinations of required dimension, and create a new matrix containing newly generated neighbhood combination columns
    
        For example, given a matrix row "CAGA", we have the neighouring sequence 5'-CA[mut]GA-3'. If upto 2-way dependent effects are considered, then transform this row into:
        "C, A, G, A, CA, CG, CA, AG, AA, GA"
        """
        features_X = np.array(['CAGA'])
        feature_comb_got = get_nbr_combs_features(features_X, feature_dim=2)
        feature_comb_exp = np.array([['C', 'A', 'G', 'A', 'CA', 'CG', 'CA', 'AG', 'AA', 'GA']])
        self.assertEqual(feature_comb_got, feature_comb_exp)
        
    def test_sort_col_by_dim(self):
        """given a dimension value, return columns in a feature matrix with
         defined dimensional neighborhood"""
        feature_X = np.array([['CtoA', 'C', 'A', 'G', 'A', 'CA', 'CG', 'CA', 'AG', 'AA', 'GA']])
        flank_window_size = 4
        single_pos_got = sort_col_by_dim(feature_X, flank_window_size, feature_dim=1)
        single_pos_exp = np.array([['C', 'A', 'G', 'A']])
        self.assertEqual(single_pos_got, single_pos_exp)
        
        two_way_pos_got = sort_col_by_dim(feature_X, flank_window_size, feature_dim=2)
        two_way_pos_exp = np.array([['CA', 'CG', 'CA', 'AG', 'AA', 'GA']])
        self.assertEqual(two_way_pos_got, two_way_pos_exp)
    
    def test_feature_selector(self):
        """depending on features specified, correctly transform the data 
        matrix into a feature matrix of selected features and response array"""
        data_matrix = np.array([['1', '0.45', 'AtoC', 'CCCC'], 
                                ['1', '0.65', 'TtoC', 'CACC'], 
                                ['1', '0.55', 'GtoA', 'TGTG'], 
                                ['1', '0.7', 'CtoT', 'AACA'], 
                                ['1', '0.3', 'AtoG', 'AGCA']])
                                
        features_got_1, responses_got_1 = feature_selector(data_matrix, feature_dim=2, directions='All')
        features_exp_1 = np.array([
        ['0.45', 'AtoC', 'C', 'C', 'C', 'C', 'CC', 'CC', 'CC', 'CC', 'CC', 'CC'],
        ['0.65', 'TtoC', 'C', 'A', 'C', 'C', 'CA', 'CC', 'CC', 'AC', 'AC', 'CC'],
        ['0.55', 'GtoA', 'T', 'G', 'T', 'G', 'TG', 'TT', 'TG', 'GT', 'GG', 'TG'],
        ['0.7', 'CtoT', 'A', 'A', 'C', 'A', 'AA', 'AC', 'AA', 'AC', 'AA', 'CA'],
        ['0.3', 'AtoG', 'A', 'G', 'C', 'A', 'AG', 'AC', 'AA', 'GC', 'GA', 'CA']])
        responses_exp_1 = np.array(['1', '1', '1','1', '1'])
        self.assertEqual(features_got_1, features_exp_1)
        self.assertEqual(responses_got_1, responses_exp_1)
        
        features_got_2, responses_got_2 = feature_selector(data_matrix, feature_dim=2, directions='None')
        features_exp_2 = np.array([
            ['0.45', 'C', 'C', 'C', 'C', 'CC', 'CC', 'CC', 'CC', 'CC', 'CC'],
            ['0.65', 'C', 'A', 'C', 'C', 'CA', 'CC', 'CC', 'AC', 'AC', 'CC'],
            ['0.55', 'T', 'G', 'T', 'G', 'TG', 'TT', 'TG', 'GT', 'GG', 'TG'],
            ['0.7', 'A', 'A', 'C', 'A', 'AA', 'AC', 'AA', 'AC', 'AA', 'CA'],
            ['0.3', 'A', 'G', 'C', 'A', 'AG', 'AC', 'AA', 'GC', 'GA', 'CA']])
        responses_exp_2 = np.array(['1', '1', '1', '1', '1'])
        self.assertEqual(features_got_2, features_exp_2)
        self.assertEqual(responses_got_2, responses_exp_2)
        
    def test_onehot_encoder(self):
        """Encode categorical features using one-hot scheme, where TRUE 
        features are encoded in 1, and FALSE features are encoded in -1."""
        data_matrix = np.array([['1', '0.45', 'AtoC', 'CCCC'], 
                                ['1', '0.65', 'TtoC', 'CACC'], 
                                ['1', '0.55', 'GtoA', 'TGTG'], 
                                ['1', '0.7', 'CtoT', 'AACA'], 
                                ['1', '0.3', 'AtoG', 'AGCA']])
                                
        enc_all_dir_got = onehot_encoder(data_matrix, flank_size=2, feature_dim=1, directions='All').tolist()
        enc_all_dir_exp = np.array([[1.0, 0.45, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0], 
        [1.0, 0.65, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0], 
        [1.0, 0.55, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0], 
        [1.0, 0.7, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0], 
        [1.0, 0.3, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0]])
        self.assertEqual(enc_all_dir_got, enc_all_dir_exp)
        
        enc_no_dir_got = onehot_encoder(data_matrix, flank_size=2, feature_dim=1, directions='None').tolist()
        enc_no_dir_exp = np.array([[1.0, 0.45, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0], 
        [1.0, 0.65, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0], 
        [1.0, 0.55, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0], 
        [1.0, 0.7, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0], 
        [1.0, 0.3, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0]])
        self.assertEqual(enc_no_dir_got, enc_no_dir_exp)
        
if __name__ == '__main__':
    main()
