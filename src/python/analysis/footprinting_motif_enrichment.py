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
parser.add_argument('footprint_features', type=str, help='Path to peak feature file (.saf, .homer)')
parser.add_argument('-b', '--background_footprints', type=str, help='Path to peak ID list to use as background features for enrichment analysis (e.g., shared peaks.)')
parser.add_argument('-r', '--region_size', type=str, default='given', help='The size of the region used for motif finding.')
parser.add_argument('-p', '--cpus', type=int, default=1, help='Number of CPUs to use.')
parser.add_argument('-s', '--motifs', type=int, default=25, help='Number of motifs to find')
parser.add_argument('-m', '--mismatches', type=int, default=2, help='Number of mismatches allowed in global optimization phase.')
parser.add_argument('-o', '--output_dir', default='.', help='Window size for merging intra and inter replicate peaks.')
parser.add_argument('--build', type=str, default='mm10', help='Genome build')
parser.add_argument('--preparsed_dir', default=None, help='')
parser.add_argument('--debug', action='store_true', help='')
args = parser.parse_args()

args.peak_features = os.path.abspath(args.footprint_features)

if not args.preparsed_dir:
    args.preparsed_dir = args.output_dir

if args.background_footprints:
    args.background_footprints = os.path.abspath(args.background_footprints)
args.output_dir = os.path.abspath(args.output_dir)

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

    if args.background_footprints:
        bg_set = args.background_footprints 
    else:
        bg_set = None 

    print("[ {} ] Running HOMER.".format(datetime.now().strftime("%b %d %H:%M:%S")))

    target_set = args.footprint_features 
    homer_outdir = run_homer(target_set, bg_set)

    print("[ {} ] Done.".format(datetime.now().strftime("%b %d %H:%M:%S")))
    print("wrote to: {}".format(homer_outdir)) 
    
            
if __name__ == '__main__':
    main()

