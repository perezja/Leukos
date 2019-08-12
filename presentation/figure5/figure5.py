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
    parser.add_argument('json', type=str, help='Figure5 JSON.')
    parser.add_argument('-o', '--output_dir', default='.', help='')
    args = parser.parse_args()

    return(args)

def a(paths, outfile):

    def label_point(x, y, val, ax):

        x = rand_jitter(x)
        y = rand_jitter(y)

        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            ax.text(point['x']+.01, point['y']+.01, str(point['val']), fontsize=7)

    def rand_jitter(arr):
        stdev = .02*(max(arr)-min(arr))
        return arr + np.random.randn(len(arr)) * stdev

    df = pd.read_csv(paths['figure5']['a'], sep='\t')

    intergenic = df[df['experiment']=='intergenic_dar_nondar']
    genelocal = df[df['experiment']=='genelocal_dar_nondar']
    de = df[df['experiment']=='de_dar_nondar']

    # DAR v. non-DAR motif enrichment categories
    data = [('DE-genes', de), ('gene-local', genelocal), ('intergenic', intergenic)]

    # set up figure grid
    f, axes = plt.subplots(3,1, num='a', figsize=(6, 14))

    for i,(label, df) in enumerate(data):

        ax = axes[i]

        df['pval_enrichment'] = -1*np.log10(df['pval_enrichment'])

        df.sort_values('pval_enrichment', ascending=False, inplace=True)
        df.reset_index(drop=True, inplace=True)
    
        sns.scatterplot(x=df.index.tolist(), y='pval_enrichment', data=df, ax=ax)
#        label_point(pd.Series(df.index.tolist()[:5]), df['pval_enrichment'][:5], df['name'][:5], ax)

        if i != 1:
            ax.set(ylabel='')

        ax.set_xticks=''
        ax.get_xaxis().set_visible(False)

        ax.set_title(label, fontsize=15)
        

    f.savefig(outfile, dpi=300)

def b(paths, outfile):

    with open(paths['figure5']['b']['gene2symbol']) as fp:
        ensbl2sym = {row.split('\t')[0]:row.split('\t')[1] for row in fp.read().strip().split('\n')}

    de_genes = pd.read_csv(paths['figure5']['b']['de_genes'], sep='\t', skiprows=1, index_col=0)
    de_genes.rename(ensbl2sym, inplace=True)

    df = pd.read_csv(paths['figure5']['b']['motif_density'], sep='\t', index_col=0)

    df = df[df.index.isin(de_genes.index)]

#    df = df[df.apply(np.sum, axis=1) > 0]
    df = df[df.columns[df.apply(np.sum, axis=0) > 0]]


    g = sns.clustermap(df, cbar_kws={'label': 'motif density'}, cmap="coolwarm", vmin=np.percentile(df.values, 2.5), vmax=np.percentile(df.values, 97.5), figsize=(11, 14))
    ax = g.ax_heatmap
#    ax.set_ylabel("")
    
    g.savefig(outfile, dpi=300)

def c(paths, outfile):


    tf_motif_paths = defaultdict(lambda: dict)
    with open(paths['figure5']['c']['tf_annot']) as fp:
        for row in fp.read().strip().split('\n'):
            tf_motif_paths[row.split('\t')[0]] = {'name':row.split('\t')[1], 'file':row.split('\t')[2]}
    with open(paths['figure5']['c']['dars']) as fp:
        dars = fp.read().strip().split('\n')
    with open(paths['figure5']['c']['ensembl2symbol']) as fp:
        ensbl2sym = {row.split('\t')[0]: row.split('\t')[1] for row in fp.read().strip().split('\n')}

    mu_atac_mi = pd.DataFrame(0, columns=['Control', 'Kdm6bKO', 'N', 't-stat', 'pval'], index=[gene_id for gene_id,attr in tf_motif_paths.items()])

    mu_exp_tf = pd.DataFrame(0, columns=['Control', 'Kdm6bKO'], index=[gene_id for gene_id,attr in tf_motif_paths.items()])


    mu_atac_peaks = pd.read_csv(paths['figure5']['c']['mean_peaks'], sep='\t', index_col=0)
    mu_expression = pd.read_csv(paths['figure5']['c']['mean_exp'], sep='\t', index_col=0)

    # calculate mean accessibility for each TF across DARs with motif-instance 

    # peak accessibility variance are assumed equal:
    # variances:
    # Control    3054.597110
    # Kdm6bKO    3048.227368

    for gene_id, attr in tf_motif_paths.items():
        df = pd.read_csv(attr['file'], sep='\t')
        binding_peaks = df['PositionID'].unique()
         
        # get motif binding peaks
        pk_cts = mu_atac_peaks[mu_atac_peaks.index.isin(binding_peaks)]
        
        # subset further for DARs
        pk_cts = pk_cts[pk_cts.index.isin(dars)]

        t, prob = stats.ttest_ind(pk_cts['Control'], pk_cts['Kdm6bKO'])

        mu_pk_cts = pk_cts.apply(np.mean, axis=0)
        
        mu_atac_mi.at[gene_id,['Control', 'Kdm6bKO']] = mu_pk_cts[['Control','Kdm6bKO']]
        mu_atac_mi.at[gene_id, 'N'] = pk_cts.shape[0]
        mu_atac_mi.at[gene_id, ['t-stat', 'pval']] = t, prob 

        # get expression of TF
        mu_exp_tf.at[gene_id, ['Control','Kdm6bKO']] = mu_expression.loc[gene_id,['Control','Kdm6bKO']] 
        
    # set up figure grid
    f, axes = plt.subplots(1,2, num='a', figsize=(9, 15), sharey=True)
    sns.set(style='whitegrid')

    f.subplots_adjust(wspace=1.5)
    
    mu_exp_tf.rename(ensbl2sym, inplace=True)
    mu_atac_mi.rename(ensbl2sym, inplace=True)

    mu_atac_mi.sort_values('Kdm6bKO', ascending=False, inplace=True)
    mu_exp_tf = mu_exp_tf.reindex(mu_atac_mi.index)

    mu_atac_mi = mu_atac_mi[['Control', 'Kdm6bKO']]

    mu_exp_tf = mu_exp_tf.applymap(np.log2)
    mu_atac_mi = mu_atac_mi.applymap(np.log2)

    sns.heatmap(mu_exp_tf, cmap='Reds', ax=axes[1]) 
    sns.heatmap(mu_atac_mi, cmap='Blues', ax=axes[0]) 

    mu_atac_mi.to_csv('atac_means.txt', sep='\t')

    axes[0].tick_params(axis='y', labelsize=15, labelcolor='#2E4053', rotation=0)
    axes[0].tick_params(axis='x', labelsize=15, labelcolor='#2E4053', rotation=90)
    axes[1].tick_params(axis='x', labelsize=15, labelcolor='#2E4053', rotation=90)

    f.savefig(outfile, dpi=300)

def d(paths,outfile): 

    g, ax = plt.subplots(1,1, figsize=(14, 16))
    g.subplots_adjust(left=0.35)

    go_results = pd.read_csv(paths['figure5']['c']['tf_go'], skiprows=1, index_col=0, sep='\t')

    go_results['Fisher.classic'] = -1*np.log10(go_results['Fisher.classic'])
    go_results.sort_values(by=['Fisher.classic'], axis=0, ascending=False, inplace=True) 
    go_results = go_results.head(n=30)


    sns.barplot(x='Fisher.classic', y='Term', data=go_results, color='grey', ci=None,ax=ax)
    ax.set(ylabel='', xlabel='-log10 pval (Fisher classic)')
    ax.tick_params(axis='y', labelsize=15, labelcolor='#2E4053')

    g.savefig(outfile, dpi=300)

def main():

    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.json) as fp:
        paths = json.load(fp)

    aof = os.path.join(args.output_dir, 'Figure5a.png')
    bof = os.path.join(args.output_dir, 'Figure5b.png')
    cof = os.path.join(args.output_dir, 'Figure5c.png')
    dof = os.path.join(args.output_dir, 'Figure5d.png')

#    a(paths, aof) 
#    b(paths, bof) 
    c(paths, cof) 
    print('Figure 5c complete.')
    plt.clf()
    d(paths, dof) 
    print('Figure 5d complete.')

if __name__ == '__main__':

    pd.set_option('display.max_columns', None)
    main()

