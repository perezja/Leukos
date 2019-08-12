#!/usr/bin/env python3
import numpy as np
import pandas as pd
import subprocess
import argparse
from datetime import datetime
import json
import tempfile
import shutil
import glob
import gzip
import os

class CombineExpression():
    def __init__(self, tpm_counts_json, subset=None):
        with open(tpm_counts_json, 'r') as f:
            GE = json.load(f)

        self.CTS = dict()
        self.TPM = dict()

        for exp, dir_name in GE.items():

            if subset:
                if exp not in subset:
                    continue

            print("* loading paths for '{}'.".format(exp))

            self.CTS[exp] = self.__get_counts_files(dir_name)
            self.TPM[exp] = self.__get_tpm_files(dir_name)

        self.__tmp = tempfile.mkdtemp(dir=args.temp_dir)

        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)

    def __get_counts_files(self, dir_name):

        cts = list() 

        cts = cts + glob.glob(os.path.join(dir_name,'*','*.ReadsPerGene.out.tab'))

        cts = pd.Series(cts, index=[i.split('/')[-2] for i in cts])
        cts.sort_index(inplace=True)

        return(cts)

    def __get_tpm_files(self, dir_name):

        tpm = list() 

        tpm = tpm + glob.glob(os.path.join(dir_name,'*','quant.sf'))

        tpm = pd.Series(tpm, index=[i.split('/')[-2] for i in tpm])
        tpm.sort_index(inplace=True)

        return(tpm)

    def __merge_region_counts(self, exp, paths):

        skip=lambda x: x in range(4) # header from STAR gct
        usecols=[0,1] # cols={gene_id, all_counts}
        sample_ids = paths.index.values

        df = pd.read_csv(paths[0], sep='\t', skiprows=skip, header=None, usecols=usecols, names=['gene_id', sample_ids[0]], index_col=0)
        df.index = df.index.str.replace(r'\.[0-9]+','')

        if df[sample_ids[0]].dtype == np.float64:
            dtype = np.float32
        elif df[sample_ids[0]].dtype == np.int64:
            dtype = np.int32
        else:
            dtype = df[sample_ids[0]].dtype

        gct_df = pd.DataFrame(0, index=df.index, columns=list(sample_ids), dtype=dtype)
        gct_df[sample_ids[0]] = df[sample_ids[0]].astype(dtype)
        for k, (i,p) in enumerate(zip(sample_ids[1:], paths[1:])):

            print("\rProcessing '{}': {}/{}".format(exp, (k+2), len(paths)), end='', flush=True)
            df = pd.read_csv(p, sep='\t', skiprows=skip, header=None, usecols=usecols, names=['gene_id', i], index_col=0)
            df.index = df.index.str.replace(r'\.[0-9]+','')
            gct_df[i] = df[i]
            
        with gzip.open(os.path.join(args.output_dir, exp+'.gct.gz'), 'wt', compresslevel=6) as f:
            f.write('{0}\t{1}\n'.format(gct_df.shape[0], gct_df.shape[1]))
            gct_df.to_csv(f, sep='\t', float_format='%.6g')

    def __merge_region_tpm(self, exp, paths):

        tpms = os.path.join(self.__tmp, exp+'_paths.txt')
        sample_ids = os.path.join(self.__tmp, exp+'_ids.txt')

        with open(tpms, 'w+') as f:
            f.write('\n'.join([i for i in paths]))
        with open(sample_ids, 'w+') as f:
            f.write('\n'.join([i for i in paths.index.values]))

        exc = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'merge_tpm.R')
        outfile = os.path.join(args.output_dir, exp+'.gtt.gz')
        cmd = ' '.join([exc, tpms, sample_ids, args.gtf, '-o', outfile]) 

        subprocess.check_call(cmd, shell=True)

    def merge_region_counts(self):

        print("[ {} ] Merging Counts.".format(datetime.now().strftime("%b %d %H:%M:%S")))

        for exp, paths in self.CTS.items():
            self.__merge_region_counts(exp, paths)

    def merge_region_tpm(self):

        print("[ {} ] Merging TPM.".format(datetime.now().strftime("%b %d %H:%M:%S")))

        for k, (i, p) in enumerate(self.TPM.items()):

            print("\rProcessing '{}': {}/{}".format(i, k+1, len(self.TPM.keys())),end='',flush=True)
            self.__merge_region_tpm(i, p)

    def __del__(self):
        if os.path.exists(self.__tmp):
            shutil.rmtree(self.__tmp)

parser = argparse.ArgumentParser(description='Run pipeline from RNA-Seq JSON file')
parser.add_argument('json', type=str, help='Path to the JSON file')
parser.add_argument('gtf', type=str, help='Path to the gtf file')
parser.add_argument('--subset', type=str, help='Subset of experiments.')
parser.add_argument('--mode', required=True, choices=['cts','tpm'], type=str, help='Merge counts or tpm expression.')
parser.add_argument('--temp_dir',type=str, default='./',help='Temporary directory')

parser.add_argument('-o','--output_dir',type=str, default='.',help='File with sample subset to process')

args = parser.parse_args()


def main():

    if args.subset:
        with open(args.subset) as fp:
            exp_subset = fp.read().strip().split('\n')
    else:
        exp_subset = None

    ce = CombineExpression(args.json, subset=exp_subset)

    if args.mode=='cts':
        ce.merge_region_counts()
    else:
        ce.merge_region_tpm()

if __name__ == "__main__":
    main()
