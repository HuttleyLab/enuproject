import os, re, glob, click
import mutation_motif
from mutation_motif.util import load_table_from_delimited_file, open_
from scitrack import CachingLogger

LOGGER = CachingLogger()

def get_files(input_dir):
    """return a list of .TXT files from a directory"""
    fns = []
    fns += glob.glob(os.path.join(input_dir, '*/*/*.txt'))
    
    if not fns:
        raise RuntimeError('Error selscting files')
    
    return fns

@click.command()
@click.option('-i','--input_dir', help='Input file containing counts data.')
@click.option('-o','--output_datafile', help='File to write output data.')
    
def main(input_dir, output_datafile):
    args = locals()
    
    output_dir = os.path.dirname(output_datafile)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    logfile_path = os.path.join(output_dir, "logs/combine_counts.log")
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    
    input_dir = os.path.dirname(input_dir)
    file_paths = get_files(input_dir)
    
    with open(output_datafile, mode='w') as out_file:
        LOGGER.output_file(output_datafile)
    
        num = 0
        for i in range(len(file_paths)):
            print ()
            print (file_paths[i])
            print ()
            with open (file_paths[i], 'r') as f:
                first_line = f.readline()
                for line in f:
                    out_file.write(line)
        
            num += 1
            if num == len(file_paths):
                break

    
    
if __name__ == "__main__":
    main()