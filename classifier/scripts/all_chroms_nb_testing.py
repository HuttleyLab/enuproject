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

def make_param_groups(flank_sizes, chroms, clf_names):
    groups = []
    for flank_size in flank_sizes:
        #get model dimension according to flank_sizes
        max_dim = 2 * flank_size + 1
        feature_dims = list(range(max_dim))
        features = product(feature_dims, chroms, clf_names)
        groups.extend([(chrom, flank_size, fd, nm) for fd, chrom, nm in features])
    return groups

#make_command(sample, args_template, resultdir, md, fsz, fd, train_sz, idx)
def make_command(exec_cmnd, args_template, chrom, flank_size, feature_dim, clf_nm):
    window_size = flank_size * 2 + 1
    testing_path = "../results/chrom%s/%d-mer/data_ratio_1/direction_All/%d_way_samples/train_size_0/sample_1/testing_1.csv.gz" % (chrom, window_size, feature_dim)
    clf_path = "../results/chrom1/%d-mer/data_ratio_1/direction_All/%d_way_samples/train_size_10255/sample_1/%s/nb_classifier.pkl" % (window_size, feature_dim, clf_nm)
    scaler_path = "../results/chrom1/%d-mer/data_ratio_1/direction_All/%d_way_samples/train_size_10255/sample_1/%s/scaler.pkl" % (window_size, feature_dim, clf_nm)
    outdir = "../results/chrom%s/%d-mer/data_ratio_1/direction_All/%d_way_samples/train_size_0/sample_1/%s" % (chrom, window_size, feature_dim, clf_nm)

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
        args = ["bernoullinb_testing",
                "--testing_file %(testing_file)s",
                "--clf_file %(clf_file)s",
                "--scaler_file %(scaler_file)s",
                "-o %(outdir)s",
                "-gc"]
        args_template = " ".join(args)
    
    if not gc_feature:
        args = ["bernoullinb_testing",
                "--testing_file %(testing_file)s",
                "--clf_file %(clf_file)s",
                "-o %(outdir)s"]
        args_template = " ".join(args)
    
    _chroms = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 'XY']
    _flank_sizes = (1, 2, 3)
    
    if gc_feature:
        _clf_names = ['NB_GC']
    if not gc_feature:
        _clf_names = ['NB']

    group = make_param_groups(flank_sizes=_flank_sizes, chroms=_chroms, clf_names=_clf_names)
    
    commands = [make_command(sample, args_template, chrom, fsz, fd, nm) for chrom, fsz, fd, nm in make_param_groups(flank_sizes=_flank_sizes, chroms=_chroms, clf_names=_clf_names)]
    
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
