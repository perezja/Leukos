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

    cmd = 'computeMatrix -S {} -R {} '
    

def b(paths, outfile):


    dars = pd.read_csv(paths['figure3']['b'], skiprows=1, sep='\t')
    print(dars['logFC'].describe())
    print(dars['ave_logCPM'].describe())

    left, width = 0.1, 0.60
    bottom, height = 0.1, 0.60
    spacing = 0.005

    rect_scatter = [left, bottom, width, height]
    rect_dens_y = [left+width, bottom, 0.10, height]
    rect_dens_x = [left, bottom+height, width, 0.10]

    f = plt.figure(figsize=(12, 12))

    ax_scatter = plt.axes(rect_scatter)
    ax_dens_x = plt.axes(rect_dens_x)
    ax_dens_x.tick_params(direction='in', labelleft=True)

    ax_dens_y = plt.axes(rect_dens_y)
    ax_dens_y.tick_params(direction='in', labelleft=True)
 
    
    sns.scatterplot(x='ave_logCPM', y='logFC', data=dars, alpha=0.5, ax=ax_scatter)
    ax_scatter.set_ylabel('log2 FC (KO / WT)', fontsize=15)
    ax_scatter.set_xlabel('Average ATAC (log2 CPM)', fontsize=15)

    sns.kdeplot(dars['logFC'], vertical=True, ax=ax_dens_y)
    sns.kdeplot(dars['ave_logCPM'],  ax=ax_dens_x)

    f.savefig(outfile, dpi=300)



#    g = sns.JointGrid(x='ave_logCPM', y='logFC', data=dars)
#
#    g = g.plot_joint(plt.scatter)
#    g = g.plot_marginals(sns.distplot, kde=True)
#
#    g.ax_joint.set_ylabel('log2 FC (KO / WT)', fontsize=15)
#    g.ax_joint.set_xlabel('Average ATAC (log2 CPM)', fontsize=15)
#
#    g.savefig(outfile, dpi=300)
#


def c(paths, outfile):

    peak_annot = pd.read_csv(paths['figure3']['c']['peak_annot'], sep='\t')

    dars = pd.read_csv(paths['figure3']['c']['dars'], skiprows=1, index_col='peak_id', sep='\t')

    peak_annot['type'] = peak_annot['peak_id'].map(lambda x: 'DAR' if x in dars.index else 'non-DAR') 

    palette = {'DAR':'#E74C3C', 'non-DAR':'#3498DB'}

    group_counts = peak_annot.groupby(['type', 'gene_region']).size().unstack(fill_value=0).reset_index()
    group_counts = pd.melt(group_counts, id_vars='type', value_name='count')
    group_counts['proportion'] = group_counts['count'] / np.sum(group_counts['count'])

    peak_annot = peak_annot.merge(group_counts, on=['gene_region', 'type'])
    peak_annot.sort_values(by='proportion', inplace=True)

    # plot 
    sns.set(style="whitegrid")

    f = sns.catplot(y='proportion', x='gene_region', hue='type', data=peak_annot, height=6, kind='bar', palette=palette) 
    f.set(ylabel='')
    f.set(xlabel='')

    f.savefig(outfile, dpi=300)

#    peak_counts = pd.read_csv(paths['figure3']['c']['accessibility'], skiprows=1, names=['peak_id', 'WT', 'KO'], sep='\t')
#    peak_counts['FC_accessibility'] = (peak_counts['KO']+1)/(peak_counts['WT']+1)
#    peak_counts['log2FC_accessibility'] = peak_counts['FC_accessibility'].apply(lambda x: np.log2(x))
#    peak_counts = peak_counts.drop(['WT', 'KO'], axis=1)
#
#
#    gene_counts = pd.read_csv(paths['figure3']['c']['expression'], skiprows=1, names=['gene', 'WT', 'KO'], sep='\t')
#    gene_counts['FC_expression'] = (gene_counts['KO']+1)/(gene_counts['WT']+1)
#    gene_counts['log2FC_expression'] = gene_counts['FC_expression'].apply(lambda x: np.log2(x))
#    gene_counts = gene_counts.drop(['WT', 'KO'], axis=1)
#
#
#    with open(paths['figure3']['c']['ensembl2symbol']) as fp:
#        ensbl2sym = {row.split('\t')[0] : row.split('\t')[1] for row in fp.read().strip().split('\n')} 
#
#    gene_counts.set_index('gene', inplace=True)
#    gene_counts.rename(ensbl2sym, inplace=True)
#
#    peak_counts = peak_counts.merge(peak_annot, on='peak_id')
#
#    reg_df = peak_counts.merge(gene_counts, on=['gene'])
#
#    reg_df = reg_df[reg_df['gene_region'] == 'promoter-proximal']
#    reg_df.to_csv('reg_df.txt', sep=' ', index=False)
#
#    g = sns.jointplot(x='log2FC_accessibility', y='log2FC_expression', data=reg_df)
##    g.set_xlim(-8,8)
##    g.set_ylim(-8,8)

   
def main():


    args = parse_args()

    pd.set_option('display.max_columns', None)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.json) as fp:
        paths = json.load(fp)

    bof = os.path.join(args.output_dir, 'Figure3b.png')
    cof = os.path.join(args.output_dir, 'Figure3c.png')

    b(paths,bof)
    c(paths,cof)

if __name__ == '__main__':
    main()

