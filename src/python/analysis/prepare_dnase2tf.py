import argparse
import subprocess
import pandas as pd
import numpy as np
from datetime import datetime
from common import bsub
import time
import tempfile
import shutil
import json
import os

def run_bam2split(bam, mem=25, gtmp=8, docker_image="apollodorus:dnase2tf:1.0.1"):  

    sample_id = os.path.split(bam)[1].split('_')[0]
    fn = sample_id + 'b.counts'

    outfile = os.path.join(os.path.abspath(tmpdir), fn)
    
    cmd = 'featureCounts -T 5' \
      + ' -a ' + saf \
      + ' -F SAF' \
      + ' -s 0' \
      + ' -o ' + outfile + ' ' \
      + bam 

    cmd = bsub(cmd, mem, gtmp, docker_image, 'featureCounts')
    po = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    return((sample_id, outfile), po)

def combine_counts(id_count_files):

    sample_ids = np.array([x[0] for x in id_count_files]) 
    paths = np.array([x[1] for x in id_count_files]) 

    i = np.argsort(sample_ids)
    sample_ids = sample_ids[i]
    paths = paths[i]

    colnames = ['peak_id', 'chr', 'start', 'end', 'strand', 'length', sample_ids[0]]
    df = pd.read_csv(paths[0], sep='\t', skiprows=2, header=None, index_col=0, names=colnames)

    pct_df = pd.DataFrame(0, index=df.index, columns=list(sample_ids), dtype=np.int32)
    pct_df[sample_ids[0]] = df[sample_ids[0]].astype(np.int32)

    for k, (i, p) in enumerate(zip(sample_ids[1:], paths[1:])):

        print('\rProcessing {}/{}'.format(k+2, len(paths)), end='', flush=True)

        df = pd.read_csv(p, sep='\t', skiprows=2, header=None, usecols=[0,6], index_col=0, names=['peak_id', i], dtype={'peak_id':str, i:np.int32})
        pct_df[i] = df[i]

    # filter peaks features with total counts under threshold 

    idx = (pct_df.sum(axis=1) >= args.count_threshold)
    print('\nfiltered {} peaks with low fragment counts.'.format(str(sum([not i for i in idx]))))

    pct_df = pct_df[idx]

    outfile = os.path.join(os.path.abspath(args.output_dir), args.prefix + '.pct')
    pct_df.to_csv(outfile, sep='\t') 

    
    return(outfile) 

parser = argparse.ArgumentParser(prog='Make read count matrix for peak feature set.')
parser.add_argument('json', type=str, help='Path to input json listing .bam and .narrowPeak input files.')
parser.add_argument('peak_features', type=str, help='Path to peak features in .saf format for featureCounts.')
parser.add_argument('prefix', type=str, help='Prefix for outfile.')
parser.add_argument('-t', '--count_threshold', type=int, default=5, help='Total fragment count threshold for containing a peak feature.')
parser.add_argument('-o', '--output_dir', default='.', help='Window size for merging intra and inter replicate peaks.')
parser.add_argument('--debug', action='store_true', help='')
args = parser.parse_args()

def main():
    
    with open(args.json) as fp:
        peak_data = json.load(fp)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)

    count_files = list()
    processes = list()

    for group, paths in peak_data.items():

        bam_files = paths.get('bams')

        for bam in bam_files:

            id_counts_tuple, po = run_feature_counts(bam, args.peak_features, tmpdir) 
            count_files.append(id_counts_tuple)
            processes.append(po)

    print("[ {} ] Running featureCounts.".format(datetime.now().strftime("%b %d %H:%M:%S")))

    while processes:
        po = processes.pop()
   
        while po.poll() is None:
            time.sleep(0.5)

    print("[ {} ] Combining feature counts into matrix.".format(datetime.now().strftime("%b %d %H:%M:%S")))

    outfile = combine_counts(count_files)

    print("[ {} ] Done.".format(datetime.now().strftime("%b %d %H:%M:%S")))
    print("wrote to: {}".format(outfile)) 
    
    shutil.rmtree(tmpdir)
            
if __name__ == '__main__':
    main()

