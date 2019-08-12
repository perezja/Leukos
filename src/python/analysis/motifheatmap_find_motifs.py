import argparse
import pandas as pd
import numpy as np
from common import bsub, Motif 
from collections import defaultdict
import concurrent.futures
import shutil
import tempfile
import subprocess
import time
import glob
import re
import json
import os

def concat_motifs(motifs):

    dfs = list() 
    for m in motifs:
        dfs.append(m.row())

    return(pd.concat(dfs, ignore_index=True))

def run_annotate_peaks(peak_bed, motif_outdir, motif_file, build, outfile):

    cmd = 'findMotifsGenome.pl ' + peak_bed + ' '\
        + build + ' '\
        + motif_outdir \
        + ' -find ' + motif_file \
        + ' -preparsedDir ' + motif_outdir \
        + ' > ' + outfile 

    cmd = bsub(cmd, 8, 4, docker_image='apollodorus/homer:mm10', job_name='homer_annot')
    subprocess.check_call(cmd, shell=True)

parser = argparse.ArgumentParser(prog='')
parser.add_argument('json', type=str, help='Path to motif.json')
parser.add_argument('-p', '--peak_bed', type=str, help='')
parser.add_argument('-n', '--num_top', type=int, default=10, help='Top N motifs to extract')
parser.add_argument('--annot_peaks', action='store_true', help='')
parser.add_argument('-o', '--output_dir', default='.', help='')
args = parser.parse_args()

args.json = os.path.abspath(args.json)
if args.annot_peaks:
    args.peak_bed = os.path.abspath(args.peak_bed)

args.output_dir = os.path.abspath(args.output_dir)

if args.annot_peaks:
    if not args.peak_bed:
        raise ValueError('"-p" peak_bed file not defined.')

def main():
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.json) as fp:
        paths = json.load(fp)

    pattern = r'motif[0-9]+(?!RV)\.motif' 

    homer_sets = defaultdict(list) 
    tmpdir = tempfile.mkdtemp(dir=args.output_dir)

    # processes list if computing annotations
    if args.annot_peaks:
        pcs = list()

    if args.annot_peaks:
        annot_outdir = os.path.join(args.output_dir, 'motif_peak_annotations')
        if not os.path.exists(annot_outdir):
            os.mkdir(annot_outdir)

    for key, dir_path in paths.items():

        prefix = key

        # filter files for standard .motif files
        motif_paths = np.array([mp for mp in glob.glob(os.path.join(dir_path,'homerResults', 'motif*.motif')) if re.search(pattern, os.path.split(mp)[1])]) 
        motif_num = np.array(list(map(int, [os.path.split(mf)[1].split('.')[0].split('motif')[1] for mf in motif_paths])))
        motifs = motif_paths[motif_num <= args.num_top]
       
        executor = concurrent.futures.ProcessPoolExecutor(max_workers=60)
        futures = list()

        # initialize motif objects
        for motif in motifs:

            m = Motif(motif, prefix)

            homer_sets[prefix].append(m)

            if args.annot_peaks:

                motif_tmpdir = tempfile.mkdtemp(dir=tmpdir)

                motif_fn = m.name+'.motif'
                motif_tmpfp = os.path.join(motif_tmpdir, motif_fn) 
                subprocess.check_call('cp {} {}'.format(motif, motif_tmpfp), shell=True) 

                outfile = os.path.join(annot_outdir, m.name+'.peakAnnotations.txt') 
                futures.append(executor.submit(run_annotate_peaks, args.peak_bed, motif_tmpdir, motif_tmpfp, 'mm10', outfile))

    concurrent.futures.wait(futures)

    if args.annot_peaks: 
        for po in futures:
            if po.result():
                print(po.result()) 

    shutil.rmtree(tmpdir)

    print('Finished motif annotations.')

    # create dataframe of motif information

    dfs = []
    motif_objs = []
    for prefix, motifs in homer_sets.items():
        df = concat_motifs(motifs)
        dfs.append(df)

        motif_objs = motif_objs + motifs

    mdf = pd.concat(dfs, ignore_index=True)

    outfile = os.path.join(args.output_dir, 'homer_motifs.txt')
    with open(outfile, 'w+') as fp:
        mdf.to_csv(fp, sep='\t', index=False)
    print('wrote to: {}'.format(outfile))

if __name__ == '__main__':
    pd.options.display.max_columns = 15
    main()

