import argparse
from collections import defaultdict
import pandas as pd
import numpy as np
import pickle
import re
import os

parser = argparse.ArgumentParser(prog='')
parser.add_argument('gene_to_peaks', type=str, help='Pickled gene to peaks dictionary.')
parser.add_argument('peak_motifs_dir', type=str, help='Output file from HOMER findMotifsGenome.pl ... "-find $MOTIF_FILE"')
parser.add_argument('prefix', type=str, help='')
parser.add_argument('-o', '--output_dir', default='.', help='')
args = parser.parse_args()

def main():
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.gene_to_peaks, 'rb') as fp:
        gene2peaks = pickle.load(fp)
    
    def motif_name(fn):
        return(fn.split('.peakAnnotations.txt')[0])

    motif_peak_annot_files = [(motif_name(i), os.path.join(args.peak_motifs_dir, i)) for i in os.listdir(args.peak_motifs_dir) if i.endswith('.txt')]

    gene_motif_density = defaultdict(lambda: defaultdict(float))

    for motif, mf in motif_peak_annot_files: 

        mdf = pd.read_csv(mf, sep='\t', index_col=0)
        for gene, peaks in gene2peaks.items():

            motif_hits = mdf[mdf.index.isin(peaks)] 

            # motif denisty = #(motif instances across gene peaks) / #(total gene peaks) 
            motif_density = motif_hits.shape[0] / len(peaks)

            gene_motif_density[motif][gene] = motif_density

    mdens_df = pd.DataFrame.from_dict(gene_motif_density)
        
    outfile = os.path.join(args.output_dir, args.prefix+'.txt')
    mdens_df.to_csv(outfile, sep='\t')

if __name__ == '__main__':
    main()

