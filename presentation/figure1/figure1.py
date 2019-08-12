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

def a(paths, outfile):

    df = pd.read_csv(paths['figure1']['a'], index_col=0, sep='\t')

    wt_pct = df[df.columns[:3]] 
    wt_pct_log2 = wt_pct.applymap(lambda x: np.log2(x+1)) 
    ko_pct = df[df.columns[3:]] 
    ko_pct_log2 = ko_pct.applymap(lambda x: np.log2(x+1)) 

    # set up figure grid
    f, axes = plt.subplots(2,3, num='a', figsize=(14, 7))
    
    an=0
    for i,x in enumerate(wt_pct.columns[:-1]):
        for j,y in enumerate(wt_pct.columns[i+1:]): 

            r,p = stats.pearsonr(wt_pct[x], wt_pct[y])

            cur_ax = axes[0,an]
            sns.scatterplot(x=x, y=y, data=wt_pct_log2, ax=cur_ax, alpha=0.1)

            cur_ax.text(12.5,14, 'p={}'.format(round(r,2)), style='italic',horizontalalignment='center', verticalalignment='center')

            an+=1

    an=0
    for i,x in enumerate(ko_pct.columns[:-1]):
        for j,y in enumerate(ko_pct.columns[i+1:]): 

            r,p = stats.pearsonr(ko_pct[x], ko_pct[y])

            cur_ax = axes[1,an]
            sns.scatterplot(x=x, y=y, data=ko_pct_log2, ax=cur_ax, alpha=0.1)

            cur_ax.text(12.5,14, 'p={}'.format(round(r,2)), style='italic', horizontalalignment='center', verticalalignment='center')

            an+=1

    f.savefig(outfile, dpi=300)

def b(paths, outfile):

    def label_point(x, y, val, ax):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            ax.text(point['x']+.02, point['y'], str(point['val']))

    def rand_jitter(arr, c):
        stdev = c*(max(arr)-min(arr))
        return arr + stdev

    pca_rnaseq = pd.read_csv(paths['figure1']['b']['rnaseq'], sep='\t')
    pca_atac = pd.read_csv(paths['figure1']['b']['atac'], sep='\t')
 
    f, axes = plt.subplots(1,2, num='b', figsize=(14, 7))
    for i, (df, title, c) in enumerate(zip([pca_atac, pca_rnaseq], ["ATAC SEQ", "RNA SEQ"], [10,0])):

        ax = axes[i]
        ax.set_title(title)
        ax.set_xlim(-0.8, 0.8)
        ax.set_ylim(-0.8, 0.8)

        sns.scatterplot(x="PC1", y="PC2", data=df, ax=ax)
#        label_point(df['PC1']+0.1, rand_jitter(df['PC2'], c), df['sample'], ax)
        label_point(df['PC1']+0.1, df['PC2']+0.05, df['sample'], ax)

    f.savefig(outfile, dpi=300)
   
def main():

    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.json) as fp:
        paths = json.load(fp)

    sns.palplot(sns.color_palette("Paired"))

    aof = os.path.join(args.output_dir, 'Figure1a.png')
    bof = os.path.join(args.output_dir, 'Figure1b.png')

    a(paths, aof) 
    b(paths, bof) 

if __name__ == '__main__':
    main()

