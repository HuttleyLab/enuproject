import os, time, re, glob, random, itertools, click
import numpy as np
import operator as op
from fractions import Fraction
from functools import reduce
from random import seed as set_seed
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from scitrack import CachingLogger

LOGGER = CachingLogger()

def str_to_numerical(num_str):
    """convert number string to float or int type"""
    if '.' in num_str:
        num = float(num_str)
    else:
        num = int(num_str)
    return num

def get_learning_data(data, flank_size):
    """get var_id, response, neighboring seqs and mutation direction info from 
    enu data or germline data"""
    records = []
    for d in data:
        var_id, location, strand, effect, allele_freqs, alleles, ref_base, var_base, f5, f3, gc_content, pep_alleles, gene_loc, gene_id, response = d.split('\t')
        seqs = f5[-flank_size:] + f3[:flank_size]
        mut_dir = ref_base + 'to' + var_base
        records.append([var_id, response, gc_content, mut_dir, seqs])
    return records

def nck(n,k):
    """return the number of possible n choose k combinations, where n>k"""
    return int(reduce(op.mul, (Fraction(n-i, i+1) for i in range(k)), 1))

def get_train_test_split(data_file, flank_size, test_size, train_size, seed, enu_ratio):
    """split data into traing set and testing set of defined sizes"""
    with open(data_file, mode='r') as df:
        data = [line.rstrip() for line in df]
        records = np.array(get_learning_data(data, flank_size))
        respose = records[:, 1]
        feature = records[:, 2:]
        train_size = str_to_numerical(train_size)
        
        if test_size:
            test_size = str_to_numerical(test_size)
        else:
            if type(train_size) is float:
                test_size = 1 - train_size
            else:
                test_size = len(data) - train_size 
        
        test_size = test_size * enu_ratio
        train_size = train_size * enu_ratio
        
        if train_size == 0.0:
            train = None
            test = records[:, 1:]
        else:
            train_features, test_features, train_responses, test_responses = train_test_split(feature, respose, test_size=test_size, train_size=train_size, random_state=seed)
            train_responses = np.array([train_responses.tolist()])
            test_responses = np.array([test_responses.tolist()])
        
            train = np.concatenate((train_responses.T, train_features), axis=1)
            test = np.concatenate((test_responses.T, test_features), axis=1)
        
        return train, test

def get_nbr_combs_features(seq_array, feature_dim):
    """read through a feature matrix row by row, transform neighboring sequences column in feature matrix into columns of neighboring base combinations of required feature_dim, and create a new matrix containing newly generated neighbhood combination columns
    
    For example, given a matrix row "CAGA", we have the neighouring sequence 5'-CA[mut]GA-3'. If upto 2-way dependent effects are considered, then transform this row into:
    "C, A, G, A, CA, CG, CA, AG, AA, GA"
    """
    all_nbr_combs = []
    for seq in seq_array:
        nbr_bases = np.array([base for base in seq])
        indices = np.arange(0, nbr_bases.shape[0])
        nbr_combs = []
        for n in range(feature_dim):
            idx_combs = [i for i in itertools.combinations(indices, n+1)]
            for idx in idx_combs:
                nbr_comb = ''.join(nbr_bases.take(list(idx), axis=0).tolist())
                nbr_combs.append(nbr_comb)
        all_nbr_combs.append(nbr_combs)
    
    all_nbr_combs = np.array(all_nbr_combs)
    return all_nbr_combs

def sort_col_by_dim(features, flank_window_size, feature_dim):
    """given a feature_dim value, return columns in a feature matrix with defined feature_dim neighborhood"""
    start_idx = sum(list(map(lambda x : nck(flank_window_size, x), [i for i in range(feature_dim)])))
    end_idx = sum(list(map(lambda x : nck(flank_window_size, x), [i + 1 for i in range(feature_dim)]))) + 1
    return features[:, start_idx: end_idx]

def feature_selector(data_matrix, feature_dim, directions):
    """depending on features considered in building a classifier, transform 
    the data matrix into a feature matrix of selected features and response
    array"""
    if directions == 'None' or 'to' in directions:
        assert feature_dim > 0, 'Dimension of the neighborhood effect has to be greater than 0.'
        response_array = data_matrix[:, 0]
        gc_array = data_matrix[:, 1]
        nbr_array = data_matrix[:, 3]
        nbr_X = get_nbr_combs_features(nbr_array, feature_dim)
        new_feature_X = np.column_stack((gc_array, nbr_X))
    
    else:
        response_array = data_matrix[:, 0]
        gc_array = data_matrix[:, 1]
        direction_array = data_matrix[:, 2]
        nbr_array = data_matrix[:, 3]
        nbr_X = get_nbr_combs_features(nbr_array, feature_dim)
        new_feature_X = np.column_stack((gc_array, direction_array, nbr_X))

    return new_feature_X, response_array
    
def label_encoder(features, flank_window_size, feature_dim, directions):
    """Encode categorical features labels with value between 0 and n_classes-1"""
    le = preprocessing.LabelEncoder()
    le_feature_list = []
    n_values = []
    
    if directions == 'All':
        mut_dir_choices = ['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT',
                           'GtoA', 'GtoC', 'GtoT', 'TtoA', 'TtoC', 'TtoG']
        n_values.append(len(mut_dir_choices))
        
        le_mut_dir = le.fit(np.array(mut_dir_choices))
        mut_dir_X = le_mut_dir.transform(features[:, 0]).tolist()
        le_feature_list.append(mut_dir_X)
    else:
        features = np.insert(features, 0, 'NA', axis=1)

    for dim in range(1, feature_dim + 1):
        base_choices = [''.join(list(i)) for i in list(itertools.product(*[['A', 'C', 'G', 'T']] * dim))]
        sorted_pos_X = sort_col_by_dim(features, flank_window_size, dim)
        
        for n in range(sorted_pos_X.shape[1]):
            le_nbr_comb = le.fit(np.array(base_choices))
            nbr_combs_X = le_nbr_comb.transform(sorted_pos_X[:, n]).tolist()
            le_feature_list.append(nbr_combs_X)
            n_values.append(len(base_choices))
            
    le_feature_X = np.array(le_feature_list).T
    n_values = np.array(n_values)
    return le_feature_X, n_values
    
def onehot_encoder(matrix, flank_size, feature_dim, directions):
    """Encode categorical features using one-hot scheme, where TRUE features are encoded in 1, and FALSE features are encoded in -1."""
    flank_window_size = flank_size * 2
    feature_X, response = feature_selector(matrix, feature_dim, directions)
    le_feature_X, n_values = label_encoder(feature_X[:, 1:], flank_window_size, feature_dim, directions)
    onehot = OneHotEncoder(n_values=n_values)
    onehot_features = onehot.fit_transform(le_feature_X).toarray()
    onehot_features[onehot_features == 0] = -1
    
    onehot_enc_X = np.concatenate((np.column_stack((response, feature_X[:, 0])), onehot_features), axis=1).astype(np.float64)
    return onehot_enc_X

@click.command()
@click.option('-S', '--seed', type=int, default=None,
              help='Seed for random number generator (e.g. 17, or 2017-04-11).'
              ' Defaults to system time.')
@click.option('--directions',
    type=click.Choice(['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT',
                       'GtoA', 'GtoC', 'GtoT', 'TtoA', 'TtoC', 'TtoG',
                       'All', 'None']), default='All',
    help='Examples of choices of mut_directions option are:\
     - All for all directions, \
     - specific direction (e.g. AtoC) for a specified direction, or \
     - None for no direction type spedified when building a classifier')
@click.option('-d','--feature_dim', type=click.IntRange(0, 21),
    default=1,
    help='Max dependent interaction order/dimension considered when constructing position features.')
@click.option('-e','--enu_prefix', required=True, help='Prefix of input file path containing ENU mutation data.')
@click.option('-g','--germline_prefix', required=True, help='Prefix of input file path containing germline mutation data.')
@click.option('-o','--ouput_path', required=True, help='Path to write output.')
@click.option('-f', '--flank_size', type=int, required=True, help='flank size considered when query the data file')
@click.option('--train_size', required=True, help='float or int. '
              'If float, should be between 0.0 and 1.0 and represent the proportion of the dataset to include in the train split. '
              'If int, represents the absolute number of train samples.')
@click.option('--test_size', default=None,
              help=' float, int, or None (default is None). '
              'If float, should be between 0.0 and 1.0 and represent the proportion of the dataset to include in the test split. '
              'If int, represents the absolute number of test samples. '
              'If None, the value is automatically set to the complement of the train size.')
@click.option('-r','--enu_ratio', type=click.Choice([1, 10, 100]), default=1, help='Ratio of ENU to germline.')
@click.option('-t','--run_time', required=True, help='Number of times to run the splitting process for.')

def main(enu_prefix, germline_prefix, ouput_path, seed, flank_size, train_size, test_size, enu_ratio, feature_dim, directions, run_time):
    args = locals()
    start_time = time.time()
    
    if 'to' in directions:
        enu_file = enu_prefix + '_%s.txt' % (directions)
        germline_file = germline_prefix + '_%s.txt' % (directions)
    else:
        enu_file = enu_prefix + '.txt'
        germline_file = germline_prefix + '.txt'
    
    for t in range(int(run_time)):
        if not seed:
            seed = int(time.time() + t)
        if seed:
            seed = int(seed + t)
        
        assert enu_ratio >= 1, "ENU ratio has to be greater than 1"
        enu_training, enu_testing = get_train_test_split(enu_file, flank_size, test_size, train_size, seed, enu_ratio)
        germline_training, germline_testing = get_train_test_split(germline_file, flank_size, test_size, train_size, seed, enu_ratio=1)

        n_train_sample = 0
        if enu_training is not None and germline_training is not None:
            training_data = np.concatenate((enu_training, germline_training), axis=0)
            try:
                onehot_training = onehot_encoder(training_data, flank_size, feature_dim, directions)
            except ValueError:
                print ('Training data has no %s data' % directions)
                quit()
            n_train_sample += onehot_training.shape[0]
        
        testing_data = np.concatenate((enu_testing, germline_testing), axis=0)

        try:
            onehot_testing = onehot_encoder(testing_data, flank_size, feature_dim, directions)
        except ValueError:
            print ('Testing data has no %s data' % directions)
            quit()

        if directions:
            output_dir = os.path.join(ouput_path, 'direction_%s/%s_way_samples/train_size_%s/' % (directions, str(feature_dim), n_train_sample), 'sample_%s/'%(t+1))
        
        else:
            output_dir = os.path.join(ouput_path, 'direction_None/%s_way_samples/train_size_%s/' % (str(feature_dim), n_train_sample), 'sample_%s/'%(t+1))
            
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        logfile_path = os.path.join(output_dir, "logs/data_sampling.log")
        LOGGER.log_file_path = logfile_path
        LOGGER.log_message(str(args), label="vars")
        LOGGER.log_message(str(seed), label="random_number_seed_%s"%(t+1))
        LOGGER.log_message(enu_file, label="ENU_input_file")
        LOGGER.log_message(germline_file, label="germline_input_file")
        
        if enu_training is not None and germline_training is not None:
            training_file = os.path.join(output_dir, 'training_%s.csv.gz'%(t+1))
            LOGGER.log_message(training_file, label="training_file_path")
            np.savetxt(training_file, onehot_training, fmt='%.3f', delimiter=",")
            print ('sample_%s/training_%s.csv.gz is saved to %s' % (t+1, t+1, output_dir))
        
        testing_file = os.path.join(output_dir, 'testing_%s.csv.gz'%(t+1))
        LOGGER.log_message(testing_file, label="testing_file_path")
        np.savetxt(testing_file, onehot_testing, fmt='%.3f', delimiter=",")
        
        print('sample_%s/testing_%s.csv.gz is saved to %s' % (t+1, t+1, output_dir))
        print()
    
    duration = time.time() - start_time
    LOGGER.log_message("%.2f" % (duration/60.), label="run duration (minutes)")
                
if __name__ == "__main__":
    main()