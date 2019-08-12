import argparse
import subprocess
import pandas as pd
import numpy as np
from collections import OrderedDict
from common import bsub
import json
import time
import re
import os

parser = argparse.ArgumentParser(prog='')
parser.add_argument('summit_coverage_files', type=str, help='List of summit coverage files')
parser.add_argument('-o', '--output_dir', default='.', help='')
args = parser.parse_args()

def main():

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.summit_coverage_files) as fp:
        sample_files = [i for i in fp.read().strip().split('\n')]

    pat = r'\w+(?=[0-9]$)' # regex for groups from sample ids

    def sid(fp):
        return(os.path.split(fp)[1].split('.')[0])

    sample_files.sort(key=lambda x: sid(x)) 
    sample_data = [(sid(x), x) for x in sample_files] 

    df = pd.read_csv(sample_data[0][1], sep='\t', header=None, names=['chr', 'start', 'end', 'peak_id', 'pos', sample_data[0][0]], index_col=['peak_id', 'pos'])

    group_summit_df = pd.DataFrame(0, index=df.index, columns = [s[0] for s in sample_data])
    group_summit_df[sample_data[0][0]] = df[sample_data[0][0]]

    for i,p in sample_data[1:]:

        df = pd.read_csv(p, sep='\t', header=None, names=['chr', 'start', 'end', 'peak_id', 'pos', i], usecols=['peak_id', 'pos', i], index_col=['peak_id', 'pos'])
        group_summit_df[i] = df[i] 

    group_summit_df.dropna(axis=1, inplace=True)
    group_summit_df.astype(np.int32)

    groups = OrderedDict()
    for i,sample_id in enumerate(list(group_summit_df.columns)):
        groups.setdefault(re.search(pat, sample_id).group(0), []).append(i)

    mu_summit_df = pd.DataFrame(0, index=group_summit_df.index, columns=groups.keys())
    for grp, idc in groups.items():
        group_cols = group_summit_df.columns[idc]
        mu_summit_df[grp] = group_summit_df[group_cols].apply(np.mean, axis=1)
        
    # log2 read density
    mu_summit_df = mu_summit_df.applymap(lambda x: np.log2(x+1))
    mu_summit_df.reset_index(inplace=True)

    for grp in groups.keys():
        outfile = os.path.join(args.output_dir, grp+'_peakSummitCoverage.txt')
        mu_summit_df[['peak_id', 'pos', grp]].to_csv(outfile, index=False, header=False, sep='\t')

    print('wrote to: {}'.format(args.output_dir))

if __name__ == '__main__':
    main()

