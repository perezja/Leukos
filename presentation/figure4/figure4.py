import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from collections import defaultdict
import pandas as pd
import numpy as np
import json
import os

def parse_args():

    parser = argparse.ArgumentParser(prog='')
    parser.add_argument('json', type=str, help='Figure1 JSON.')
    parser.add_argument('-o', '--output_dir', default='.', help='')
    args = parser.parse_args()

    return(args)

def a(paths, outfile):

    de_genes = pd.read_csv(paths['figure4']['a']['de_genes'], skiprows=1, sep='\t', index_col=0)
    dars = pd.read_csv(paths['figure4']['a']['dars'], skiprows=1, sep='\t', index_col=0)

    de_genes.rename(columns={'logFC':'logFC_expression'}, inplace=True)
    dars.rename(columns={'logFC':'logFC_accessibility'}, inplace=True)

    pa = pd.read_csv(paths['figure4']['a']['peak_annot'], sep='\t')

    with open(paths['figure4']['a']['ensembl2symbol']) as fp:
        ensbl2gene = {row.split('\t')[0] :row.split('\t')[1] for row in fp.read().strip().split('\n')}

    de_genes.rename(ensbl2gene, inplace=True)

    pa = pa[pa['gene'].isin(de_genes.index)]
    pa = pa[pa['peak_id'].isin(dars.index)]

    pa = pd.merge(pa, dars['logFC_accessibility'], left_on='peak_id', right_index=True)
    pa = pd.merge(pa, de_genes['logFC_expression'], left_on='gene', right_index=True)

    pa['peak_type'] = pa['logFC_accessibility'].apply(lambda x: 'KO' if (x > 0) else 'WT')
    
    # get peak ranks of logFC accessibility per { gene { peak_type }} group
    ko_open = pa[pa['logFC_accessibility'] > 0] 
    ko_open['rank'] = ko_open.groupby(['gene', 'peak_type'])['logFC_accessibility'].rank('dense') 

    wt_open = pa[pa['logFC_accessibility'] < 0] 
    wt_open['rank'] = wt_open.groupby(['gene', 'peak_type'])['logFC_accessibility'].rank('dense', ascending=False) 
    wt_open['rank'] = wt_open['rank'].apply(lambda x: -1 * x) 

    pa = pd.concat([ko_open, wt_open], ignore_index=True)

    pa.sort_values('logFC_expression', inplace=True)
    pa.dropna(inplace=True)
    pa.to_csv('pa.txt', sep=' ')

    ## Plot ##

    sns.set(style='whitegrid')

    f, axes = plt.subplots(2, 1, num='a', figsize=(15, 15), sharex=True)
    
    ## Fig4a. DARs associated with DE genes 

    # top 

    sns.barplot(x='gene', y='logFC_expression', data=pa, color='#7F8C8D', ax=axes[0], ci=None)
    axes[0].set_ylabel('logFC expression', fontsize=20)
    axes[0].set(xlabel='')

    # bottom (i)

    sns.scatterplot(x='gene', y='rank', data=pa, style='gene_region', hue='logFC_accessibility', s=65, ax=axes[1], palette='coolwarm')

    axes[1].set_ylabel('DARs', fontsize=20)

    axes[1].set_yticks([])
    axes[1].tick_params(labelsize=20)

    plt.subplots_adjust(bottom=0.20)
    plt.xticks(rotation=90)

    f.tight_layout()
    f.suptitle('DARs associated with DE genes', fontsize=30)

    sns.despine(left=True, bottom=True)

    f.savefig(outfile, dpi=300)
   
def main():

    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.json) as fp:
        paths = json.load(fp)

    aof = os.path.join(args.output_dir, 'Figure4a.png')

    a(paths, aof) 

if __name__ == '__main__':
    main()

#    palette = {"DARs":"#E67E22", "Non-DARs":"#3498DB"}
