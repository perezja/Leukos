import os
import sys
import glob
import json
import logging
import pandas as pd
import argparse

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)

def parse_arguments(debug=False):

    parser = argparse.ArgumentParser(prog='Get pair mates from FASTQ directory.', description='')
    parser.add_argument('fastq_dir', type=str)
    parser.add_argument('sample_lookup', type=str, help='TSV file linking sample name to tab-separated list of run files.') 
    parser.add_argument('--prefix', type=str, default='input', help='TSV file linking sample name to tab-separated list of run files.') 
    parser.add_argument('--not_paired', action='store_true', help='') 
    parser.add_argument('--subset', type=str, default=None, help='List of run ids to process.') 
    parser.add_argument('--output_dir', '-o', default='.', type=str, help='Output directory.')
    args = parser.parse_args()

    return args

def get_sample_map(sample_lookup, subset=None):

    # sample_lookup: run id -> sample id
    sample_dict = { i : list() for i in sample_lookup.iloc[:,0] }
    for n, x in enumerate(sample_lookup.index.tolist()):

        sample_id = sample_lookup.iloc[n,0]

        if subset and (x not in subset):
            try:
                sample_dict.pop(sample_id)
            except KeyError:
                continue
        else:
            try:
                sample_dict[sample_id].append(x)
            except KeyError:
                continue

    return(sample_dict)

def main():

    args = parse_arguments()
    
    if args.subset:
        with open(args.subset) as fp:
            args.subset = fp.read().strip().split('\n')

    df = pd.read_csv(args.sample_lookup, index_col=0, sep='\t')

    sample_to_run = get_sample_map(df, args.subset)

    fqs = os.listdir(args.fastq_dir)
    cur_run_ids = {i.split('_')[0] for n, i in enumerate(fqs) if ((i.endswith('.fq') or i.endswith('.fastq') or i.endswith('.gz')) and (i.split('_')[0] not in fqs[:n])) } 
    
    sample_to_fastqs = dict()
    for sample_id, all_run_ids in sample_to_run.items():

        sample_to_fastqs[sample_id] = list()

        for rid in all_run_ids:

            assert rid in cur_run_ids 

            run_files = glob.glob(os.path.join(args.fastq_dir, rid + '*'))
            run_files.sort()

            if args.not_paired:
                if not len(run_files) == 1: raise ValueError
            else:
                if not len(run_files) == 2: raise ValueError
                
            sample_to_fastqs[sample_id].append(run_files)

    outfile = os.path.join(args.output_dir, args.prefix + '.fastq.')  
    if args.not_paired:
        outfile += 'se.json'
    else:
        outfile += 'pe.json'

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(outfile, 'w+') as fp:
        json.dump(sample_to_fastqs, fp, sort_keys=True, indent=4, separators=(',', ': '))

    print('wrote to: {}'.format(outfile)) 
            
if __name__ == '__main__':
    main()
