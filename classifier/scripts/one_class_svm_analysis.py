import os, glob, gzip, pickle, time, click, csv, json
import pandas as pd
import numpy as np
from sklearn import svm, metrics
from sklearn.metrics import roc_auc_score
from cogent3.util.misc import open_

from scitrack import CachingLogger

LOGGER = CachingLogger()
        
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
                 
@main.command()
@_training_file
@_output_dir
def ocs_training(training_file, output_dir):
    """tests a logistic regression model"""
    args = locals()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logfile_path = os.path.join(output_dir, "logs/ocs_train.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    LOGGER.log_message(training_file, label="Training data file")
    
    training_data = np.loadtxt(training_file, delimiter=",")
    # germline(reference) is classified as class -1 in our original data, but for one-class, it should be labeled as 1, so we need to change the labelling here, I did it by original_label * -1.
    germline_feature = training_data[training_data[:,0] == -1][:, 2:]
    germline_response = training_data[training_data[:,0] == -1][:, 0] * -1
    
    # fit the model
    clf_ocs = svm.OneClassSVM(nu=0.3, kernel="linear")
    clf_ocs.fit(germline_feature)
    model_file = os.path.join(output_dir, 'one_class_svm.pkl')
    with open(model_file, 'wb') as clf_file:
        LOGGER.log_message(model_file, label="Pickled one class svm model")
        pickle.dump(clf_ocs, clf_file)
    print ('one_class_svm.pkl is saved to output directory!')

 
@main.command()
@_clf_file
@_testing_file
@_output_dir
def ocs_testing(testing_file, clf_file, output_dir):    
    args = locals()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logfile_path = os.path.join(output_dir, "logs/ocs_test.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    
    LOGGER.log_message(clf_file, label="Logistic regression model")
    with open(clf_file, 'rb') as clf:
        clf_ocs = pickle.load(clf)
        
    data_reader = make_reader(testing_file, chunk_size=1000) 
    LOGGER.log_message(testing_file, label="Testing data file")
    predict_func = clf_ocs.predict
    score_func = clf_ocs.decision_function
    
    expects = np.array([])
    predictions = np.array([])
    predicted_scores = np.array([])
    for testing_sample in data_reader:
        feature_testing = testing_sample[:, 2:]
        #change the labelling here
        response_testing = testing_sample[:,0].astype(np.float) * -1
        
        exp = response_testing
        pred = clf_ocs.predict(feature_testing)
        scr = clf_ocs.decision_function(feature_testing).flatten()

        expects = np.concatenate((expects, exp), axis=0)
        predictions = np.concatenate((predictions, pred), axis=0)
        predicted_scores = np.concatenate((predicted_scores, scr), axis=0)
        
    result_dict = {}
    response_classes = ['Germline', 'ENU']
    report = metrics.classification_report(expects, predictions, target_names=response_classes)
    print (report)
    report = report_to_json(report)
    result_dict.update(report)
    
    confusion_matrix = pd.DataFrame(metrics.confusion_matrix(expects, predictions), index=[i for i in response_classes], columns=[i for i in response_classes]).rename_axis('actual / predicted', axis=1)
    print (confusion_matrix)
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