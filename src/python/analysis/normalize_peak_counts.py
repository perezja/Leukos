import argparse
import subprocess
import pandas as pd
import numpy as np
from rnaseqnorm import normalize_quantiles
from collections import OrderedDict
import re
import os

parser = argparse.ArgumentParser(prog='Quantile normalize peak count matrix and take mean across replicates.')
parser.add_argument('count_matrix', type=str, help='Path to count matrix.')
parser.add_argument('groups', type=str, help='List of group levels for each sample column.')
parser.add_argument('prefix', type=str, help='Prefix for outfile.')
parser.add_argument('-o', '--output_dir', default='.', help='Window size for merging intra and inter replicate peaks.')
parser.add_argument('--debug', action='store_true', help='')
args = parser.parse_args()

def main():
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    count_matrix_df = pd.read_csv(args.count_matrix, sep='\t', index_col=0)

    print("Normalizing...")
    norm_matrix_df = normalize_quantiles(count_matrix_df)
    print(norm_matrix_df.head())


    print("Averaging replicate counts...")

    pat = r'\w+(?=[0-9]+$)' # regex for groups from sample ids 
    groups = OrderedDict()
    for i,sid in enumerate(list(norm_matrix_df.columns)):
        groups.setdefault(re.search(pat, sid).group(0), []).append(i)

    # get mean across group indices
    mu_norm_matrix_df = pd.DataFrame(0, index=norm_matrix_df.index, columns=groups.keys())

    for grp, idc in groups.items():
        group_cols = norm_matrix_df.columns[idc]
        mu_norm_matrix_df[grp] = norm_matrix_df[group_cols].apply(np.mean, axis=1)

    outfile = os.path.join(args.output_dir, args.prefix+'.accessibility.txt')
    mu_norm_matrix_df.to_csv(outfile, sep='\t')

    outfile = os.path.join(args.output_dir, args.prefix+'.quant_norm.pct')
    norm_matrix_df.to_csv(outfile, sep='\t')

    print("wrote to *.accessibility (averaged) and *.quant_norm.pct (w/ replicates) to: {}".format(args.output_dir)) 
    
if __name__ == '__main__':
    main()

