import argparse
from read_coverage import run_coverage
from common import bsub
import concurrent.futures
import subprocess
import tempfile
import shutil
import time
import re
import os

parser = argparse.ArgumentParser(prog='')
parser.add_argument('peak_bed', type=str, help='BED file with peaks of interest.')
parser.add_argument('peak_motif_list', type=str, help='Directory with atac bams.')
parser.add_argument('fasta', type=str, help='Pickled FASTA in binary format.')
parser.add_argument('-p', '--peak_ids', type=str, help='List of peak IDs to subset BED.')
parser.add_argument('-o', '--output_dir', default='.', help='')
args = parser.parse_args()

args.peak_bed = os.path.abspath(args.peak_bed)
args.output_dir = os.path.abspath(args.output_dir)

def run_find_motif_bed(find_motif_file):

    script = os.path.dirname(os.path.abspath('footprinting_motif_bed.py'))
    cmd = 'python3 ' + script \
        + ' ' + args.peak_bed \
        + ' ' + find_motif_file \
        + ' ' + args.fasta 
        + ' -o ' + args.output_dir
         
    cmd = bsub(cmd, mem=8, gtmp=4, docker_image='apollodorus/bioinf:pr', 'motif_scan')
    print(cmd)
#    subprocess.check_call(cmd, shell=True)

def main():


    with open(args.peak_motif_list) as fp:
        motif_files = fp.read().strip().split('\n')

    futures = list()
    executor = concurrent.futures.ProcessPoolExecutor(max_workers=60)
    
    for mf in motif_files:
        futures.append(executor.submit(run_find_motif_bed, mf)) 

    for po in futures:
        if po.result():
            print(po.result())
    
