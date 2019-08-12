import argparse
import seaborn as sn
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import OrderedDict
import json
import re
import os

parser = argparse.ArgumentParser(prog='')
parser.add_argument('json', type=str, help='Path to motif coverage files.')
parser.add_argument('-o', '--output_dir', default='.', help='')
args = parser.parse_args()


# return dictionary of mean read coverage between groups
# {(peak_id, pos) : {'chr':'', 'start':'', 'end':'', 'group1':'' group2:''}}
def coverage_dictionary(paths):

    def sample_id(fp):
        return(os.path.split(fp)[1].split('_')[0])

    sample_files = [(sample_id(fp), fp) for group, l in paths.items() for fp in l] 
    sample_files.sort(key=lambda x: x[0])
    
    df = pd.read_csv(sample_files[0][1], sep='\t', header=None, names=['chr', 'start', 'end', 'motif_id', 'pos', sample_files[0][0]], index_col=['motif_id','chr','start','end','pos'])

    counts_df = pd.DataFrame(0, index=df.index, columns=[s[0] for s in sample_files], dtype=np.float32) 
    counts_df[sample_files[0][0]] = df

    for i,p in sample_files[1:]:

        df = pd.read_csv(p, sep='\t', header=None, names=['chr', 'start', 'end', 'motif_id', 'pos', i], index_col=['motif_id','chr','start','end','pos'])

        counts_df[i] = df[i] 

    pat = r'\w+(?=[0-9]+$)' # regex for groups from sample ids

    groups = OrderedDict()
    for i,sid in enumerate(list(counts_df.columns)):
        groups.setdefault(re.search(pat, sid).group(0), []).append(i)

    # get mean across group indices
    mu_counts_df = pd.DataFrame(0, index=counts_df.index, columns=groups.keys())

    for grp, idc in groups.items():
        group_cols = counts_df.columns[idc]
        mu_counts_df[grp] = counts_df[group_cols].apply(np.mean, axis=1)

    mu_counts_df.reset_index(inplace=True)
    mu_counts_df[['start', 'end']] = mu_counts_df[['start', 'end']].astype(np.int32)

    return(mu_counts_df)

def main():
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.json) as fp:
        paths = json.load(fp)

    motif_coverage_df = coverage_dictionary(paths)
    motif_coverage_df.to_csv('motif_coverage.tmp', sep='\t')
    mu_coverage = motif_coverage_df.groupby('pos').mean()
    mu_coverage.reset_index(inplace=True)

    sns_plot = sn.barplot(x='pos', y='Control', data=mu_coverage, color="skyblue", alpha=0.5, label="Control")
    sns_plot = sn.barplot(x='pos', y='Kdm6bKO', data=mu_coverage, color="red", alpha=0.5, label="Kdm6bKO")

    plt.setp(sns_plot.get_xticklabels(), fontsize=3)

    plt.xlabel('')
    plt.ylabel('mean insertions')


    sns_plot = sns_plot.get_figure()
    sns_plot.savefig(os.path.join(args.output_dir,'HIC_coverage.png'), dpi=300)

    # get peak sequences from 1-based peak coordinates.

if __name__ == '__main__':
    main()

