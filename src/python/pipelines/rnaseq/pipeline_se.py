#!/usr/bin/env python3
import argparse
import subprocess
import contextlib
import shlex
import tempfile
import shutil
import os
import sys
import logging
from datetime import datetime

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

parser = argparse.ArgumentParser(description='RNA-seq pipeline')
parser.add_argument('prefix', type=str, help='Prefix for output file (i.e. sample name).')
parser.add_argument('genome_dir', type=str, help='Path to directory storing star index files')
parser.add_argument('transcripts', type=str, help='Path to transcripts fasta file')
parser.add_argument('gtf', type=str, help='Path to gff annotation file used for building STAR and Salmon indices.')
parser.add_argument('-f', '--fqs', nargs='+', type=str, help='')
parser.add_argument('-o', '--output_dir', type=str, default='./', help='Path to output directory for all processes')
parser.add_argument('--tmp_dir', type=str, default=None, help='Path to directory for storing tmp files')

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

if args.tmp_dir is None:
    args.tmp_dir = os.path.abspath(tempfile.mkdtemp(dir=args.output_dir)) 

print("Sample: " + args.prefix)

src_dir = os.path.dirname(os.path.realpath(__file__)) 

print('Starting STAR') 

cmd = os.path.join(src_dir,'run_STAR2.py ') \
    + args.genome_dir + ' ' \
    + args.prefix  \
    + ' --fqs ' + ' '.join(args.fqs) \
    + ' -o '+args.output_dir 

print(' cmd * ' + cmd)
subprocess.check_call(cmd, shell=True)

# 4. Mark Duplicates

cmd = os.path.join(src_dir,'run_MarkDuplicates.py ') \
    + os.path.join(args.output_dir, args.prefix \
        + '.Aligned.sortedByCoord.out.bam ') \
    + args.prefix \
    + ' --tmp_dir ' + args.tmp_dir \
    + ' -o ' + args.output_dir

print(' cmd * '+cmd)
subprocess.check_call(cmd, shell=True)

# 5. Quantify

cmd = os.path.join(src_dir,'run_Salmon.py ') \
    + os.path.join(args.output_dir, args.prefix \
        + '.Aligned.toTranscriptome.out.bam ') \
    + args.transcripts + ' ' \
    + args.gtf \
    + ' -o '+args.output_dir

print(' cmd * ' + cmd)
subprocess.check_call(cmd, shell=True)

print("Sample '" + args.prefix + "'" + " processed.")

# 6. Post-processing: Delete unnecessary bams 

if os.path.exists(args.tmp_dir):
    shutil.rmtree(args.tmp_dir)
