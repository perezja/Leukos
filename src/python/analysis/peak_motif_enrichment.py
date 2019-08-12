import argparse
import subprocess
import pandas as pd
import numpy as np
from datetime import datetime
from common import bsub 
import re
import tempfile
import shutil
import json
import os

parser = argparse.ArgumentParser(prog='HOMER motif enrichment analysis.')
parser.add_argument('peak_features', type=str, help='Path to peak feature file (.saf, .homer)')
parser.add_argument('peak_ids', type=str, help='Path to peak ID list to test for motif enrichment.')
parser.add_argument('build', type=str, help='Genome build (e.g., mm10)')
parser.add_argument('-b', '--background_ids', type=str, help='Path to peak ID list to use as background features for enrichment analysis (e.g., shared peaks.)')
parser.add_argument('-r', '--region_size', type=str, default='given', help='The size of the region used for motif finding.')
parser.add_argument('-p', '--cpus', type=int, default=1, help='Number of CPUs to use.')
parser.add_argument('-s', '--motifs', type=int, default=25, help='Number of motifs to find')
parser.add_argument('-m', '--mismatches', type=int, default=2, help='Number of mismatches allowed in global optimization phase.')
parser.add_argument('-o', '--output_dir', default='.', help='Window size for merging intra and inter replicate peaks.')
parser.add_argument('--preparsed_dir', default=None, help='')
parser.add_argument('--debug', action='store_true', help='')
args = parser.parse_args()

args.peak_features = os.path.abspath(args.peak_features)
args.peak_ids = os.path.abspath(args.peak_ids)

if not args.preparsed_dir:
    args.preparsed_dir = args.output_dir

if args.background_ids:
    args.background_ids = os.path.abspath(args.background_ids)
args.output_dir = os.path.abspath(args.output_dir)

def prepare_features(peak_features, peak_ids, outdir, background_ids=None):

    pf_tmp = tempfile.NamedTemporaryFile(dir=outdir, delete=False)
    with open(peak_ids) as fp:
        peak_ids = fp.read().strip().split('\n')

    df = pd.read_csv(peak_features, sep='\t', header=None, names=['peak_id', 'chr', 'start', 'end', 'strand'])
    pdf = df[df['peak_id'].isin(peak_ids)]

    pdf.to_csv(pf_tmp.name, sep='\t', header=False, index=False)

    if background_ids:
        with open(background_ids) as fp:
            background_ids = fp.read().strip().split('\n')
        bg_tmp = tempfile.NamedTemporaryFile(dir=outdir, delete=False)
        bdf = df[df['peak_id'].isin(background_ids)]

        bdf.to_csv(bg_tmp.name, sep='\t', header=False, index=False)
        return(pf_tmp.name, bg_tmp.name)

    return(pf_tmp.name, None)

def run_homer(target_set, bg_set=None, mem=25, gtmp=8, docker_image="apollodorus/homer:mm10"):

# peak_features format:
#  1 - Peak ID, 2- chromosome, 3- start, 4-end, 5-strand

    cmd = 'findMotifsGenome.pl '+target_set + ' '\
      + args.build + ' ' \
      + args.output_dir \
      + ' -size '+args.region_size \
      + ' -mis ' + str(args.mismatches) \
      + ' -S ' + str(args.motifs) \
      + ' -p ' + str(args.cpus) \
      + ' -preparsedDir ' + args.preparsed_dir

    if bg_set:
        cmd = cmd + ' -bg ' + bg_set

    cmd = bsub(cmd, mem, gtmp, docker_image, 'annovar', args.debug)
    subprocess.check_call(cmd, shell=True)

    return(args.output_dir)

def main():
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)

    print("[ {} ] Preparing input.".format(datetime.now().strftime("%b %d %H:%M:%S")))   

    if args.background_ids:
        enrichment_set = prepare_features(args.peak_features, args.peak_ids, tmpdir, args.background_ids) 
    else:
        enrichment_set = prepare_features(args.peak_features, args.peak_ids, tmpdir) 

    print("[ {} ] Running HOMER.".format(datetime.now().strftime("%b %d %H:%M:%S")))

    target_set = enrichment_set[0]
    bg_set = enrichment_set[1]
    homer_outdir = run_homer(target_set, bg_set)

    print("[ {} ] Done.".format(datetime.now().strftime("%b %d %H:%M:%S")))
    print("wrote to: {}".format(homer_outdir)) 
    
    shutil.rmtree(tmpdir)
            
if __name__ == '__main__':
    main()

