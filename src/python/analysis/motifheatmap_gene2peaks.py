import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import pickle
import os

parser = argparse.ArgumentParser(prog='Generate gene to peak dictionary.')
parser.add_argument('peak_annot', type=str, help='Peak gene-annotation file')
parser.add_argument('dars', type=str, help='List of DAR peak ids')
parser.add_argument('prefix', type=str, help='')
parser.add_argument('-o', '--output_dir', default='.', help='')
args = parser.parse_args()

args.output_dir = os.path.abspath(args.output_dir)

def main():
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.dars) as fp:
        dars = fp.read().strip().split('\n')

    gene2peak = defaultdict(list) 

    df = pd.read_csv(args.peak_annot, sep='\t', index_col=[0,5])
    for (gene, peak_id), row in df.iterrows():

        if peak_id not in dars:
            continue

        gene2peak[gene].append(peak_id)

    outfile = os.path.join(args.output_dir, args.prefix+'.pkl')
    with open(outfile, 'wb') as fp:
        pickle.dump(gene2peak, fp)

if __name__ == '__main__':
    main()

