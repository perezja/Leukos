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

    counts = pd.read_csv(paths['figure2']['a']['counts'],  index_col=0, sep='\t')
    de =  pd.read_csv(paths['figure2']['a']['de_genes'], skiprows=1, sep='\t')
    de = de['gene_id']

    with open(paths['figure2']['a']['ensembl2symbol']) as fp:
        ensbl2sym = {row.split('\t')[0]: row.split('\t')[1] for row in fp.read().strip().split('\n')}

    counts = counts.applymap(lambda x:np.log2(x+1))
    counts = counts[counts.index.isin(de)]
    counts = counts.rename(ensbl2sym)


    g = sns.clustermap(counts, cbar_kws={'label': 'log2(tmm)'}, cmap="coolwarm", vmin=np.percentile(counts.values, 2.5), vmax=np.percentile(counts.values, 97.5), figsize=(5,12))

    ax = g.ax_heatmap
    ax.set_ylabel("")
    
    g.savefig(outfile, dpi=300)

def b(paths, outfile):


    de_genes = pd.read_csv(paths['figure2']['b']['de_genes'], skiprows=1, sep='\t')

    upreg_genes = de_genes[de_genes['logFC'] > 0]['gene_id'].tolist()
    downreg_genes = de_genes[de_genes['logFC'] < 0]['gene_id'].tolist()

    go_results = pd.read_csv(paths['figure2']['b']['go_results'], skiprows=1, index_col=0, sep='\t')

    go_results['Fisher.classic'] = -1*np.log10(go_results['Fisher.classic'])
    go_results.sort_values(by=['Fisher.classic'], axis=0, ascending=False, inplace=True) 

    go2gene = defaultdict(list)
    with open(paths['figure2']['b']['ensembl2go']) as fp:
        fp.readline() # skip header
        for row in fp.read().strip().split('\n'):
            go2gene[row.split('\t')[1]].append(row.split('\t')[0])
#    print(dict(list(go2gene.items())[0:1]))

    go_uptally = defaultdict() 
    go_downtally = defaultdict() 

    for go_id, row in go_results.iterrows():

        go_uptally[go_id] = sum([gene in upreg_genes for gene in go2gene[go_id]])
        go_downtally[go_id] = sum([gene in downreg_genes for gene in go2gene[go_id]])

    go_results['Upregulated'] = go_results.index.map(go_uptally)
    go_results['Downregulated'] = go_results.index.map(go_downtally)

    go_tally = pd.melt(go_results, id_vars=['Term', 'Fisher.classic'], value_vars=['Downregulated', 'Upregulated'], var_name='Type', value_name='#Genes')
    go_results.to_csv('go_results.txt', sep='\t')

    sns.set(style='whitegrid')

    f, axes = plt.subplots(1, 2, num='b', figsize=(14, 16), sharey=True)
    f.subplots_adjust(left=0.35)
    
    # GO Analysis (pvalue) 
    sns.barplot(x='Fisher.classic', y='Term', data=go_results, color='grey', ax=axes[0], ci=None)
    axes[0].set(ylabel='', xlabel='-log10 pval (Fisher classic)')
    axes[0].tick_params(axis='y', labelsize=15, labelcolor='#2E4053')

    # Go Analysis (DE genes)
    sns.set_color_codes('pastel')
    palette = {"Upregulated":"#E74C3C", "Downregulated":"#3498DB"}
    sns.barplot(x='#Genes', y='Term', data=go_tally, hue='Type', ax=axes[1], ci=None, palette=palette)
    axes[0].set(ylabel='')
    axes[1].legend(ncol=1, loc='lower right', frameon=True)
    axes[1].legend(ncol=1, loc='lower right', frameon=True)

    f.suptitle('Gene set enrichment (GSEA)')
    sns.despine(left=True, bottom=True)

    f.savefig(outfile, dpi=300)
   
def main():

    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.json) as fp:
        paths = json.load(fp)

    aof = os.path.join(args.output_dir, 'Figure2a.png')
    bof = os.path.join(args.output_dir, 'Figure2b.png')

    a(paths, aof) 
    b(paths,bof)

if __name__ == '__main__':
    main()

