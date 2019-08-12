import argparse
import subprocess
import pandas as pd
import numpy as np
from rnaseqnorm import edgeR_cpm 
from collections import OrderedDict
import re
import os

parser = argparse.ArgumentParser(prog='Normalize gene count matrix using TMM procedure and take mean across replicates.')
parser.add_argument('count_matrix', type=str, help='Path to gene count matrix.')
parser.add_argument('groups', type=str, help='List of group levels for each sample column.')
parser.add_argument('prefix', type=str, help='Prefix for outfile.')
parser.add_argument('-o', '--output_dir', default='.', help='')
parser.add_argument('--debug', action='store_true', help='')
args = parser.parse_args()

def main():
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    counts_df = pd.read_csv(args.count_matrix, sep='\t', skiprows=1, index_col=0)

    print("Normalizing...")

    norm_counts_df = edgeR_cpm(counts_df)

    print("Averaging replicate counts...")

    pat = r'\w+(?=[0-9]+$)' # regex for groups from sample ids 
    groups = OrderedDict()
    for i,sid in enumerate(list(norm_counts_df.columns)):
        groups.setdefault(re.search(pat, sid).group(0), []).append(i)

    # get mean across group indices
    mu_norm_counts_df = pd.DataFrame(0, index=norm_counts_df.index, columns=groups.keys())

    for grp, idc in groups.items():
        group_cols = norm_counts_df.columns[idc]
        mu_norm_counts_df[grp] = norm_counts_df[group_cols].apply(np.mean, axis=1)

    outfile = os.path.join(args.output_dir, args.prefix+'mean_expression.tmm.txt')
    mu_norm_counts_df.to_csv(outfile, sep='\t')

    outfile = os.path.join(args.output_dir, args.prefix+'.tmm.gct')
    norm_counts_df.to_csv(outfile, sep='\t')

    print("wrote to *.expression.txt (averaged) and *.tmm.gct (w/ replicates) to: {}".format(args.output_dir)) 
    
if __name__ == '__main__':
    main()

