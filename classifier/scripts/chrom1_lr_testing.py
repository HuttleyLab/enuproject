#!/usr/bin/env python
from warnings import filterwarnings
filterwarnings("ignore", ".+not using MPI.+")

import click
import os
import sys
import subprocess
from itertools import product
os.environ['DONT_USE_MPI'] = "1"
from cogent3.util import parallel


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0

    Parameters
    ----------

    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""
    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        msg = err
        sys.stderr.writelines("FAILED: %s\n%s" % (cmnd, msg))
        exit(proc.returncode)
    if out is not None:
        r = out.decode('utf8')
    else:
        r = None
    return r

def make_param_groups(flank_sizes, train_sizes, sample_idxs, clf_names):
    groups = []
    for flank_size in flank_sizes:
        #get model dimension according to flank_sizes
        max_dim = 2 * flank_size + 1
        feature_dims = list(range(max_dim))
        features = product(feature_dims, train_sizes, sample_idxs, clf_names)
        groups.extend([(flank_size, fd, train_sz, idx, nm) for fd, train_sz, idx, nm in features])
    return groups

#make_command(sample, args_template, resultdir, md, fsz, fd, train_sz, idx)
def make_command(exec_cmnd, args_template, flank_size, feature_dim, train_size, sample_id, clf_nm):
    window_size = flank_size * 2 + 1
    testing_path = "../results/chrom1/%d-mer/data_ratio_1/direction_All/%d_way_samples/train_size_%d/sample_%d/testing_%d.csv.gz" % (window_size, feature_dim, train_size, sample_id, sample_id)
    clf_path = "../results/chrom1/%d-mer/data_ratio_1/direction_All/%d_way_samples/train_size_%d/sample_%d/%s/logreg_classifier.pkl" % (window_size, feature_dim, train_size, sample_id, clf_nm)
    scaler_path = "../results/chrom1/%d-mer/data_ratio_1/direction_All/%d_way_samples/train_size_%d/sample_%d/%s/scaler.pkl" % (window_size, feature_dim, train_size, sample_id, clf_nm)
    outdir = "../results/chrom1/%d-mer/data_ratio_1/direction_All/%d_way_samples/train_size_%d/sample_%d/%s" % (window_size, feature_dim, train_size, sample_id, clf_nm)

    vals = dict(testing_file=testing_path, clf_file=clf_path, scaler_file=scaler_path, outdir=outdir)
                
    cmnd = " ".join([exec_cmnd, args_template % vals])
    return cmnd


@click.command()
@click.option('-n', '--numprocs', type=int, default=1,
              help='number of procs if using multiprocess')
@click.option('-t', '--test', is_flag=True, help='test run')
@click.option('-gc','--gc_feature', is_flag=True, help='Indlude GC content as a feature.')
def main(numprocs, test, gc_feature):
    script_path = os.path.join(os.path.dirname(__file__),
                               "classification_analysis.py")
    assert os.path.exists(script_path), "classification_analysis.py not found"
    sample = "python %s" % script_path
    if gc_feature:
        args = ["logreg_testing",
                "--testing_file %(testing_file)s",
                "--clf_file %(clf_file)s",
                "--scaler_file %(scaler_file)s",
                "-o %(outdir)s",
                "-gc"]
        args_template = " ".join(args)
    
    if not gc_feature:
        args = ["logreg_testing",
                "--testing_file %(testing_file)s",
                "--clf_file %(clf_file)s",
                "-o %(outdir)s"]
        args_template = " ".join(args)
    
    _flank_sizes = (1, 2, 3)
    _train_sizes = (1009, 2050, 4101, 10255, 16408)
    _sample_idxs = (1, 2, 3, 4, 5)
    if gc_feature:
        _clf_names = ['LR_GC']
    if not gc_feature:
        _clf_names = ['LR']

    group = make_param_groups(flank_sizes=_flank_sizes, train_sizes=_train_sizes, sample_idxs=_sample_idxs, clf_names=_clf_names)
    
    commands = [make_command(sample, args_template, fsz, fd, train_sz, idx, nm) for fsz, fd, train_sz, idx, nm in make_param_groups(flank_sizes=_flank_sizes, train_sizes=_train_sizes, sample_idxs=_sample_idxs, clf_names=_clf_names)]
    
    if test and numprocs > 1:
        commands = commands[:numprocs]
    elif test and parallel.size > 1:
        commands = commands[:parallel.size]
        
    if test:
        print(commands)

    if numprocs > 1:
        parallel.use_multiprocessing(numprocs)

    for r in parallel.imap(exec_command, commands):
        pass


if __name__ == "__main__":
    main()
