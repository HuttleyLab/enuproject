import os, click
from scitrack import CachingLogger

LOGGER = CachingLogger()

@click.command()
@click.option('-i','--input_file', help='Input file containing variants data')
@click.option('--direction', help='Specified mutation direction')

def main(input_file, direction):
    args = locals()
    
    output_dir = os.path.dirname(input_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
      
    logfile_path = os.path.join(output_dir, "logs/sample_chromXY_%s.log" % direction)
    LOGGER.log_file_path = logfile_path
    LOGGER.log_message(str(args), label="vars")
    
    outfile_name = os.path.basename(input_file).split('.')[0] + '_%s.txt' % direction
    output_file = os.path.join(output_dir, outfile_name)
    
    with open(input_file, mode='r') as df:
        LOGGER.input_file(input_file)
        with open(output_file, mode='w') as output:
            LOGGER.output_file(output_file)
            
            data = [line.rstrip() for line in df]
            num = 0
            for d in data:
                var_id, location, strand, effect, allele_freqs, alleles, ref_base, var_base, f5, f3, pep_alleles, gene_loc, gene_id, response = d.split('\t')
                mut_dir = ref_base + 'to' + var_base
                if mut_dir == direction:
                    output.write(d + '\n')
                    num += 1
                    
            LOGGER.log_message('%s' % num, label="%s variant count" % direction)
                    
        
        print("num written", num)
        output.close()


if __name__ == "__main__":
    main()