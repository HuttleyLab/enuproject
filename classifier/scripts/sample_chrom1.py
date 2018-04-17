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

def make_param_groups(flank_sizes, train_sizes):
    groups = []
    for flank_size in flank_sizes:
        #get model dimension according to flank_sizes
        feature_dimensions = list(range(2 * flank_size + 1))
        feature_by_tsz = product(feature_dimensions, train_sizes)
        groups.extend([(flank_size, fd, tsz) for fd, tsz in feature_by_tsz])
    return groups

def make_command(exec_cmnd, args_template, basedir, flank_size, feature_dim, train_size):
    window_size = flank_size * 2 + 1
    outdir = os.path.join(basedir, "chrom1/%d-mer/data_ratio_1/" % (window_size))
    # these filenames actually have a .txt
    enu_path = "../data/ENU_variants/SNVs20151101_chrom1"
    germ_path = "../data/germline_variants/mouse_germline_All_88_chrom1"
    
    vals = dict(enu_data=enu_path, germ_data=germ_path, seed=1493876701,
                num_reps=5, flank_size=flank_size, feature_dim=feature_dim,
                train_size=train_size, outdir=outdir)
                
    cmnd = " ".join([exec_cmnd, args_template % vals])
    return cmnd


@click.command()
@click.option('-r', '--resultdir', default="../results",
              help='parent dir to write results')
@click.option('-n', '--numprocs', type=int, default=1,
              help='number of procs if using multiprocess')
@click.option('-t', '--test', is_flag=True, help='test run')
def main(resultdir, numprocs, test):
    script_path = os.path.join(os.path.dirname(__file__),
                               "sample_train_test.py")
    assert os.path.exists(script_path), "sample_train_test.py not found"
    sample = "python %s" % script_path
    
    _flank_sizes = (1, 2, 3)
    _train_sizes = [0.029, 0.0589, 0.1178, 0.2945, 0.4712]

    args = ["-e %(enu_data)s",
            "-g %(germ_data)s",
            "-o %(outdir)s",
            "-S %(seed)s",
            "--train_size=%(train_size)s",
            "--feature_dim %(feature_dim)d",
            "-t %(num_reps)s",
            "--flank_size %(flank_size)d"]
    args_template = " ".join(args)
    
    commands = [make_command(sample, args_template, resultdir, fsz, fd, tsz)
                for fsz, fd, tsz in
                make_param_groups(flank_sizes=_flank_sizes, train_sizes=_train_sizes)]
                
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
