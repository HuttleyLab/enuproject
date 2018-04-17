import os, glob, gzip, pickle, time, click, csv, json, collections, itertools
import pandas as pd
import numpy as np
from random import seed as set_seed
from sklearn import metrics, linear_model
from sklearn.metrics import roc_auc_score, accuracy_score, brier_score_loss
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import BernoulliNB
from sklearn.model_selection import ShuffleSplit, GridSearchCV
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.neighbors.kde import KernelDensity
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from cogent3.util.misc import open_

from scitrack import CachingLogger

LOGGER = CachingLogger()

def get_positions(flank_size):
    """return list containing all positions according to flank size"""
    pos_num = flank_size * 2
    positions = []
    for i in range(pos_num):
        positions.append('Pos_%s' % i)
    return positions

def get_scaled_data(training_data, output_dir):
    feature = training_data[:,1:]
    response = training_data[:,0].astype(np.float)
        
    scaler = StandardScaler()
    scaler.fit(feature)
    scaler_file = os.path.join(output_dir, "scaler.pkl")
    with open(scaler_file, 'wb') as norm:
        LOGGER.log_message(scaler_file, label="Pickled scaler")
        pickle.dump(scaler, norm)
    print ('scaler.pkl is saved to output directory!')
        
    scaled_feature = scaler.transform(feature)
    training_sample = np.insert(scaled_feature, 0, response, axis=1)
    
    return training_sample

def get_feature_X(training_data, gc_feature=None):
    """return correct feature matrix
    if with gc_feature:
        include gc column
    if without gc_feature:
        do not include gc column
    """
    if gc_feature:
        print('Include GC feature')
        feature_matrix = training_data[:,1:]
    if not gc_feature:
        print('No GC feature')
        feature_matrix = training_data[:,2:]
    return feature_matrix
        
def make_reader(infile_path, chunk_size=1000, sep=","):
    records = []
    with open_(infile_path) as infile:
        for line in infile:
            line = line.strip()
            if not line:  # handle empty lines
                break

            line = line.split(sep)
            records.append(line)

            if len(records) == chunk_size:
                yield np.array(records, dtype=float)

                records = []  # then reset records to be empty

    if records:  # catch the last group
        yield np.array(records, dtype=float)

def make_prediction(data, output_dir, pd_func, score_func, gc_feature=None, scaler_file=None):
    if not gc_feature:
        feature_matrix = data[:, 2:]
        response = data[:,0].astype(np.float)
            
        expected = response
        predicted = pd_func(feature_matrix)
        predicted_score = score_func(feature_matrix)
     
    if gc_feature:
        if scaler_file:
            scaler = pickle.load(open(scaler_file, 'rb'))
            LOGGER.log_message(scaler_file, label="Scaler")
            feature_matrix = data[:,1:]
            scaled_feature_matrix = scaler.transform(feature_matrix)
        if not scaler_file:
            scaled_feature_matrix = data[:,1:]
        response = data[:,0].astype(np.float)
        
        expected = response
        predicted = pd_func(scaled_feature_matrix)
        predicted_score = score_func(scaled_feature_matrix)
    
    return expected, predicted, predicted_score

def report_to_json(report):
    """convert string type report to json"""
    items = ['']
    for i in ''.join(report).split('\n'):
        if i is '':
            continue
        for j in i.split('  '):
            if j is '':
                continue
            items.append(j.strip())
    data = np.reshape(items, (4, 5))
    df = pd.DataFrame(data=data[1:,1:], index=data[1:,0], columns=data[0,1:])
    json_report = json.loads(df.to_json())
    return json_report

@click.group()
def main():
    pass

_training_file = click.option('--training_file', help='Input file containing training data.')
_testing_file = click.option('--testing_file', help='Input file containing testing data.')
_output_dir = click.option('-o','--output_dir', help='Directory to write output data.')
_clf_file = click.option('--clf_file', help='Pickled file containing saved classifier.')
_scaler_file = click.option('--scaler_file', default=None, help='Pickled file containing saved scaler.')
_gc_feature = click.option('-gc','--gc_feature', is_flag=True, help='Indlude GC content as a feature.')
_flank_size = click.option('-f', '--flank_size', type=int, required=True, help='flank size considered when query the data file')
_c_options = click.option('-C', '--c_options', default='0.1,1,10,100', help='C values choosed for model, for example: 0.1,1,10,100')
_penalty_options = click.option('-P', '--penalty_options', default='l1', help="penalty parameter choosed for model, Examples of choices of options are: 'l1','l2', or 'l1,l2'")
_alpha_options = click.option('-a', '--alpha_options', default="0.01,0.1,1,2,3", help='Alpha values choosed for model, for example: 0.01,0.1,1,2,3')
_seed = click.option('-S', '--seed', type=int, default=None,
              help='Seed for random number generator (e.g. 17, or 2017-04-11).'
              ' Defaults to system time.')
                
@main.command()
@_seed
@_training_file
@_output_dir
@_gc_feature
@_c_options
@_penalty_options
@_flank_size
def logreg_train_and_validate(training_file, output_dir, gc_feature, seed, c_options, penalty_options, flank_size):
    """Doing nested cross validation on training data, and choose parameter C and penalty in the algorithm that return the best performance"""
    args = locals()
    if not seed:
        seed = int(time.time())
    set_seed(seed)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logfile_path = os.path.join(output_dir, "logs/logreg-cv.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    LOGGER.log_message(str(seed), label="random_number_seed")
    LOGGER.log_message(penalty_options, label="Penalty options")
    LOGGER.log_message(str(c_options), label="C value options")
    LOGGER.log_message(training_file, label="Training data file")
    
    training_data = np.loadtxt(training_file, delimiter=",")
    if gc_feature:
        training_data = get_scaled_data(training_data, output_dir)
    feature_matrix = get_feature_X(training_data, gc_feature=gc_feature)
    response = training_data[:,0].astype(np.float)
    
    rs = ShuffleSplit(n_splits=1, test_size=0.5, random_state=seed)
    cv_dir = os.path.join(output_dir, "cross_validation")
    if not os.path.exists(cv_dir):
        os.makedirs(cv_dir)
    for train_index, validate_index in rs.split(training_data):
        cv_training = training_data.take(train_index, axis=0)
        cv_training_file = os.path.join(cv_dir, 'cv_training.csv')
        LOGGER.log_message(cv_training_file, label="cv_training_file_path")
        np.savetxt(cv_training_file, cv_training, fmt='%.3f', delimiter=",")
        print ('cross_validation/cv_training.csv is saved to output directory!')
        
        cv_validation = training_data.take(validate_index, axis=0)
        validation_file = os.path.join(cv_dir, 'cv_validation.csv')
        LOGGER.log_message(validation_file, label="cv_validation_file_path")
        np.savetxt(validation_file, cv_validation, fmt='%.3f', delimiter=",")
        print ('cross_validation/cv_validation.csv is saved to output directory!')
    
    penalty_list = penalty_options.split(',')
    c_list = [float(c) for c in c_options.split(',')]
    param_grid = {'C': c_list, 'penalty': penalty_list}
    
    log_reg = LogisticRegression(class_weight='balanced')#imbalanced sampling
    logreg_classifier = GridSearchCV(estimator=log_reg, param_grid=param_grid, scoring='roc_auc', cv=rs)
    logreg_classifier.fit(feature_matrix, response)
    cv_results = pd.DataFrame.from_dict(logreg_classifier.cv_results_)
    best_logreg = logreg_classifier.best_estimator_
    
    model_file = os.path.join(output_dir, 'logreg_classifier.pkl')
    with open(model_file, 'wb') as clf_file:
        LOGGER.log_message(model_file, label="Pickled logistic regression model")
        pickle.dump(best_logreg, clf_file)
    print ('logreg_classifier.pkl is saved to output directory!')
    
    cv_results_dict = {}
    cv_rank = json.loads(cv_results[['param_C', 'param_penalty', 'params', 'rank_test_score', 'split0_test_score']].to_json())
    cv_results_dict.update(cv_rank)
    best_c = {'best_c_choice' : logreg_classifier.best_estimator_.C}
    cv_results_dict.update(best_c)
    cv_report_file = os.path.join(output_dir, 'cv_report.json')
    with open(cv_report_file, 'w') as report_file:
        LOGGER.log_message(cv_report_file, label="Cross validation report file")
        report_file.write(json.dumps(cv_results_dict))
    print ('cv_report.json is saved to output directory!')
    
    dimension, _, _ = training_file.split('/')[-4].split('_')
    _, mut_dir_type = training_file.split('/')[-5].split('_')
    
    feature_col_nm = []
    if gc_feature:
        feature_col_nm.append('GC_feature')
    
    if mut_dir_type == 'All':
        mut_dir_choices = ['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT',
                           'GtoA', 'GtoC', 'GtoT', 'TtoA', 'TtoC', 'TtoG']
        feature_col_nm.extend(mut_dir_choices)
    
    positions = get_positions(int(flank_size))
    for dim in range(1, int(dimension) + 1):
        position_combs = [':'.join(list(i)) for i in list(itertools.combinations(positions, dim))]
        base_choices = [''.join(list(i)) for i in list(itertools.product(*[['A', 'C', 'G', 'T']] * dim))]
        pos_nbr_comb = ['_'.join(i) for i in list(itertools.product(position_combs, base_choices))]
        feature_col_nm.extend(pos_nbr_comb)    
    
    beta_list = best_logreg.coef_.tolist()[0]
    beta_dict = dict(zip(feature_col_nm, beta_list))
    
    coef_dict_file = os.path.join(output_dir, 'betas.json')
    with open(coef_dict_file, 'w') as dict_file:
        LOGGER.log_message(coef_dict_file, label="Betas saved to file")
        dict_file.write(json.dumps(beta_dict))
    print ('betas.json is saved to output directory!')
    print ('=' * 50)

@main.command()
@_output_dir
@_gc_feature
def logreg_on_train(output_dir, gc_feature):
    """apply the logistic regression model got after the cross validation 
       back to the training sample used in the cross validation 
       process
    """
    args = locals()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logfile_path = os.path.join(output_dir, "logs/logreg-on-train.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    
    logreg_clf = os.path.join(output_dir, 'logreg_classifier.pkl')
    LOGGER.log_message(logreg_clf, label="Logistic regression model")
    with open(logreg_clf, 'rb') as clf_file:
        logreg_classifier = pickle.load(clf_file)
        
    cv_training = os.path.join(output_dir, "cross_validation/cv_training.csv")
    LOGGER.log_message(cv_training, label="CV training data file")
    data_reader = make_reader(cv_training, chunk_size=1000)
    predict_func = logreg_classifier.predict
    score_func = logreg_classifier.decision_function
    
    expects = np.array([])
    predictions = np.array([])
    predicted_scores = np.array([])
    for testing_sample in data_reader:
        exp, pred, scr = make_prediction(testing_sample, output_dir, predict_func, score_func, gc_feature=gc_feature, scaler_file=None)
        expects = np.concatenate((expects, exp), axis=0)
        predictions = np.concatenate((predictions, pred), axis=0)
        predicted_scores = np.concatenate((predicted_scores, scr), axis=0)

    result_dict = {}
    response_classes = ['Germline', 'ENU']
    report = metrics.classification_report(expects, predictions, target_names=response_classes)
    report = report_to_json(report)
    result_dict.update(report)
    
    confusion_matrix = pd.DataFrame(metrics.confusion_matrix(expects, predictions), index=[i for i in response_classes], columns=[i for i in response_classes]).rename_axis('actual / predicted', axis=1)
    confusion_matrix = json.loads(confusion_matrix.to_json())
    result_dict.update(confusion_matrix)
    au_roc_score = {'auROC' : roc_auc_score(expects, predicted_scores)}
    result_dict.update(au_roc_score)
    
    output_file = os.path.join(output_dir, 'clf_on_train.json')
    with open(output_file, 'w') as report_file:
        LOGGER.log_message(output_file, label="Result report file")
        report_file.write(json.dumps(result_dict))
        
    print ()
    print ('clf_on_train.json is saved to output directory!')
    print ('=' * 50)
 
@main.command()
@_testing_file
@_clf_file
@_scaler_file
@_output_dir
@_gc_feature
def logreg_testing(testing_file, clf_file, scaler_file, output_dir, gc_feature):
    """tests a logistic regression model"""
    args = locals()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logfile_path = os.path.join(output_dir, "logs/logreg-test.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    
    LOGGER.log_message(clf_file, label="Logistic regression model")
    with open(clf_file, 'rb') as clf:
        logreg_classifier = pickle.load(clf)
    
    data_reader = make_reader(testing_file, chunk_size=1000) 
    LOGGER.log_message(testing_file, label="Testing data file")
    predict_func = logreg_classifier.predict
    score_func = logreg_classifier.decision_function
    
    expects = np.array([])
    predictions = np.array([])
    predicted_scores = np.array([])
    for testing_sample in data_reader:
        exp, pred, scr = make_prediction(testing_sample, output_dir, predict_func, score_func, gc_feature=gc_feature, scaler_file=scaler_file)
        expects = np.concatenate((expects, exp), axis=0)
        predictions = np.concatenate((predictions, pred), axis=0)
        predicted_scores = np.concatenate((predicted_scores, scr), axis=0)
    
    result_dict = {}
    response_classes = ['Germline', 'ENU']
    report = metrics.classification_report(expects, predictions, target_names=response_classes)
    report = report_to_json(report)
    result_dict.update(report)
    
    confusion_matrix = pd.DataFrame(metrics.confusion_matrix(expects, predictions), index=[i for i in response_classes], columns=[i for i in response_classes]).rename_axis('actual / predicted', axis=1)
    confusion_matrix = json.loads(confusion_matrix.to_json())
    result_dict.update(confusion_matrix)
    au_roc_score = {'auROC' : roc_auc_score(expects, predicted_scores)}
    result_dict.update(au_roc_score)
    
    output_file = os.path.join(output_dir, 'classification_report.json')
    with open(output_file, 'w') as report_file:
        LOGGER.log_message(output_file, label="Result report file")
        report_file.write(json.dumps(result_dict))
        
    print ()
    print ('classification_report.json is saved to output directory!')
    print ('=' * 50)

################################################################################

@main.command()
@_seed
@_training_file
@_output_dir
@_gc_feature
@_alpha_options
def bernoullinb_train_and_validate(training_file, output_dir, gc_feature, seed, alpha_options):
    """Doing nested cross validation on training data, and choose the best alpha value that return the best classification performance"""
    args = locals()
    if not seed:
        seed = int(time.time())
    set_seed(seed)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logfile_path = os.path.join(output_dir, "logs/nb-cv.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    LOGGER.log_message(str(seed), label="random_number_seed")
    LOGGER.log_message(str(alpha_options), label="Alpha value options")
    LOGGER.log_message(training_file, label="Training data file")
    
    training_data = np.loadtxt(training_file, delimiter=",")
    if gc_feature:
        training_data = get_scaled_data(training_data, output_dir)
    feature_matrix = get_feature_X(training_data, gc_feature=gc_feature)
    response = training_data[:,0].astype(np.float)
    
    rs = ShuffleSplit(n_splits=1, test_size=0.5, random_state=seed)
    cv_dir = os.path.join(output_dir, "cross_validation")
    if not os.path.exists(cv_dir):
        os.makedirs(cv_dir)
    for train_index, validate_index in rs.split(training_data):
        cv_training = training_data.take(train_index, axis=0)
        cv_training_file = os.path.join(cv_dir, 'cv_training.csv')
        LOGGER.log_message(cv_training_file, label="cv_training_file_path")
        np.savetxt(cv_training_file, cv_training, fmt='%.3f', delimiter=",")
        print ('cross_validation/cv_training.csv is saved to output directory!')
        
        cv_validation = training_data.take(validate_index, axis=0)
        validation_file = os.path.join(cv_dir, 'cv_validation.csv')
        LOGGER.log_message(validation_file, label="cv_validation_file_path")
        np.savetxt(validation_file, cv_validation, fmt='%.3f', delimiter=",")
        print ('cross_validation/cv_validation.csv is saved to output directory!')
    
    a_list = [float(a) for a in alpha_options.split(',')]
    param_grid = {'alpha': a_list}
    
    nb = BernoulliNB()
    nb_classifier = GridSearchCV(estimator=nb, param_grid=param_grid, scoring='roc_auc', cv=rs)
    nb_classifier.fit(feature_matrix, response)
    cv_results = pd.DataFrame.from_dict(nb_classifier.cv_results_)
    best_nb = nb_classifier.best_estimator_
    
    model_file = os.path.join(output_dir, 'nb_classifier.pkl')
    with open(model_file, 'wb') as clf_file:
        LOGGER.log_message(model_file, label="Pickled NB model")
        pickle.dump(best_nb, clf_file)
    print ('nb_classifier.pkl is saved to output directory!')
    
    cv_results_dict = {}
    cv_rank = json.loads(cv_results[['param_alpha', 'params', 'rank_test_score', 'split0_test_score']].to_json())
    cv_results_dict.update(cv_rank)
    best_a = {'best_a_choice' : nb_classifier.best_estimator_.alpha}
    cv_results_dict.update(best_a)
    
    cv_report_file = os.path.join(output_dir, 'cv_report.json')
    with open(cv_report_file, 'w') as report_file:
        LOGGER.log_message(cv_report_file, label="Cross validation report file")
        report_file.write(json.dumps(cv_results_dict))
    print ('cv_report.json is saved to output directory!')

@main.command()
@_output_dir
@_gc_feature
def bernoullinb_on_train(output_dir, gc_feature):
    """apply the naive bayes model got after the cross validation 
       back to the training sample used in the cross validation 
       process"""
    args = locals()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logfile_path = os.path.join(output_dir, "logs/nb-on-train.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    
    nb_clf = os.path.join(output_dir, "nb_classifier.pkl")
    LOGGER.log_message(nb_clf, label="Naive Bayes model")
    with open(nb_clf, 'rb') as clf_file:
        nb_classifier = pickle.load(clf_file)
        
    cv_training = os.path.join(output_dir, "cross_validation/cv_training.csv")
    LOGGER.log_message(cv_training, label="CV training data file")
    
    data_reader = make_reader(cv_training, chunk_size=1000)
    predict_func = nb_classifier.predict
    score_func = nb_classifier.predict_proba
    
    expects = np.array([])
    predictions = np.array([])
    predicted_scores = np.array([])
    for testing_sample in data_reader:
        exp, pred, scr = make_prediction(testing_sample, output_dir, predict_func, score_func, gc_feature=gc_feature, scaler_file=None)
        scr = scr[:, 1]
        expects = np.concatenate((expects, exp), axis=0)
        predictions = np.concatenate((predictions, pred), axis=0)
        predicted_scores = np.concatenate((predicted_scores, scr), axis=0)
    
    result_dict = {}
    response_classes = ['Germline', 'ENU']
    report = metrics.classification_report(expects, predictions, target_names=response_classes)
    report = report_to_json(report)
    result_dict.update(report)
    
    confusion_matrix = pd.DataFrame(metrics.confusion_matrix(expects, predictions), index=[i for i in response_classes], columns=[i for i in response_classes]).rename_axis('actual / predicted', axis=1)
    confusion_matrix = json.loads(confusion_matrix.to_json())
    result_dict.update(confusion_matrix)
    au_roc_score = {'auROC' : roc_auc_score(expects, predicted_scores)}
    result_dict.update(au_roc_score)
    
    output_file = os.path.join(output_dir, 'clf_on_train.json')
    with open(output_file, 'w') as report_file:
        LOGGER.log_message(output_file, label="Result report file")
        report_file.write(json.dumps(result_dict))
        
    print ()
    print ('clf_on_train.json is saved to output directory!')
    print ('=' * 50)

@main.command()  
@_testing_file
@_clf_file
@_scaler_file
@_output_dir
@_gc_feature
def bernoullinb_testing(testing_file, clf_file, scaler_file, output_dir, gc_feature):
    """tests a naive bayes model"""
    args = locals()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logfile_path = os.path.join(output_dir, "logs/nb-test.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    
    LOGGER.log_message(clf_file, label="Naive Bayes model")
    with open(clf_file, 'rb') as clf:
        nb_classifier = pickle.load(clf)
    
    data_reader = make_reader(testing_file, chunk_size=1000)    
    LOGGER.log_message(testing_file, label="Testing data file")
    predict_func = nb_classifier.predict
    score_func = nb_classifier.predict_proba
    
    expects = np.array([])
    predictions = np.array([])
    predicted_scores = np.array([])
    for testing_sample in data_reader:
        exp, pred, scr = make_prediction(testing_sample, output_dir, predict_func, score_func, gc_feature=gc_feature, scaler_file=None)
        scr = scr[:, 1]
        expects = np.concatenate((expects, exp), axis=0)
        predictions = np.concatenate((predictions, pred), axis=0)
        predicted_scores = np.concatenate((predicted_scores, scr), axis=0)
    
    result_dict = {}
    response_classes = ['Germline', 'ENU']
    report = metrics.classification_report(expects, predictions, target_names=response_classes)
    report = report_to_json(report)
    result_dict.update(report)
    
    confusion_matrix = pd.DataFrame(metrics.confusion_matrix(expects, predictions), index=[i for i in response_classes], columns=[i for i in response_classes]).rename_axis('actual / predicted', axis=1)
    confusion_matrix = json.loads(confusion_matrix.to_json())
    result_dict.update(confusion_matrix)
    au_roc_score = {'auROC' : roc_auc_score(expects, predicted_scores)}
    result_dict.update(au_roc_score)
    
    output_file = os.path.join(output_dir, 'classification_report.json')
    with open(output_file, 'w') as report_file:
        LOGGER.log_message(output_file, label="Result report file")
        report_file.write(json.dumps(result_dict))
        
    print ()
    print ('classification_report.json is saved to output directory!')
    print ('=' * 50)

if __name__ == "__main__":
    main()    