import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
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

def b(paths, outfile):

    dar_enrich = pd.read_csv(paths['figure6']['b']['dar_enrichment'], sep='\t')

    fp_enrich = pd.read_csv(paths['figure6']['b']['footprint_enrichment'], sep='\t')

    f, axes = plt.subplots(1,2, num='b', figsize=(12, 6))

    fp_logp = fp_enrich['pval_enrichment'].map(lambda x: -1*np.log10(x))
    fp_logp = fp_logp.rename('footprint enrichments') 

    dar_logp = dar_enrich['pval_enrichment'].map(lambda x: -1*np.log10(x))
    dar_logp.sort_values(ascending=False, inplace=True)
    dar_logp = dar_logp.rename('top DAR enrichments') 
    dar_logp = dar_logp[:10]

    sns.set_style("whitegrid")

    sns.kdeplot(dar_logp, shade=True, color="#E74C3C", ax=axes[0])
    sns.kdeplot(fp_logp, shade=True, color="#3498DB", ax=axes[0])

    axes[0].set_xlabel('-log10 pval', fontsize=15)

    def label_point(x, y, val, ax):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            ax.text(point['x']+.02, point['y'], str(point['val']), fontsize=10)

    def rand_jitter(arr, c):
        stdev = c*(max(arr)-min(arr))
        return arr + stdev

    fp_enrich['pval_enrichment'] = -1*np.log10(fp_enrich['pval_enrichment'])
    fp_enrich.sort_values('pval_enrichment', ascending=False, inplace=True)
    fp_enrich.reset_index(drop=True, inplace=True)

    sns.scatterplot(x=fp_enrich.index.tolist(), y='pval_enrichment', data=fp_enrich, ax=axes[1]) 

#    label_point(pd.Series(fp_enrich.index.tolist()[:10]), fp_enrich['pval_enrichment'][:10], fp_enrich['name'][:10], axes[1])
    axes[1].set_xticks=''

    f.savefig(outfile, dpi=300)

def c(paths, outfile):

    fp_enrich = pd.read_csv(paths['figure6']['c'], sep='\t')
    hic_hit = fp_enrich[fp_enrich['name']=='ZNF416-Zf']
    hic_df = pd.melt(hic_hit, id_vars=None, value_vars=['target_freq', 'bg_freq'], var_name='enrichment group', value_name='% total footprints')
    hic_df.sort_values('enrichment group', inplace=True)

    sns.set_style("whitegrid")
    f, axes = plt.subplots(1,1, num='c', figsize=(12, 12))

    palette = ['#ABB2B9','#A569BD']
    sns.barplot(x='enrichment group', y='% total footprints', data=hic_df, palette=palette, ax=axes)

    axes.set_xlabel('', fontsize=15)
    axes.set_xticks = ''
    axes.set_xticklabels([]) 

    axes.set_ylabel('')

    f.savefig(outfile, dpi=300)

def main():


    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.json) as fp:
        paths = json.load(fp)

    bof = os.path.join(args.output_dir, 'Figure6b.png')
    cof = os.path.join(args.output_dir, 'Figure6c.png')

    b(paths, bof) 
    c(paths, cof) 

if __name__ == '__main__':
    main()

