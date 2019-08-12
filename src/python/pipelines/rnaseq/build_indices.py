#!/usr/bin/env python3

import argparse
import os
import subprocess
from datetime import datetime
import contextlib

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

parser = argparse.ArgumentParser(description='Build STAR indices for RNAseq pipeline.')
parser.add_argument('fasta', type=str, help='reference genome FASTA file.')
parser.add_argument('annot', type=str, help='gtf annotation file.')
parser.add_argument('transcripts', type=str, help='FASTA of cdna transcripts in reference genome.')
parser.add_argument('-l', '--read_length', type=int, default=100, help='Output directory')
parser.add_argument('-o', '--output_dir',default='genomeDir', help='Output directory')
parser.add_argument('--lsf', action='store_true', help='Output directory')
parser.add_argument('--threads', type=int, default=4, help="")

args = parser.parse_args()

def bsub(cmd, job_name, queue="research-hpc", mem=35, gtmp=8, docker_image="apollodorus/brain-eqtl:rnaseq2"):
    """ Creates a bsub command for LSF job submission.

    Args:
        cmd: command to be run.
        queue: queue to submit job
        job_name: name of job
        mem: RAM space to request from scheduler (MB)
        gtmp: temporary space to request from scheduler (GB)
        docker_image: [username]/[repos]:[tag]

    Returns:
        string: bsub command string.
    """

    mem = str(mem)
    gtmp = str(gtmp)
    bsub = "bsub -Is -q "+queue \
        + " -J "+job_name \
        + " -M "+mem+"000000" \
        + " -R 'select[mem>"+mem+"000 && gtmp > "+gtmp+"]" \
        + " rusage[mem="+mem+"000, gtmp="+gtmp+"]'" \
        + " -a 'docker("+docker_image+")'" \
        + " /bin/bash -c '"+cmd+"'" 

    return(bsub)


if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
    
exc = '/opt/STAR-2.6.0a/bin/Linux_x86_64/STAR'

cmd = exc + ' --runMode genomeGenerate ' \
     + ' --runThreadN ' + str(args.threads) \
     + ' --genomeDir ' + args.output_dir \
     + ' --genomeFastaFiles ' + args.fasta \
     + ' --sjdbGTFfile ' + args.annot \
     + ' --sjdbOverhang ' + str(args.read_length - 1)

if args.lsf:
    cmd = bsub(cmd, job_name='build_star_index')

print('Building STAR index')

subprocess.check_call(cmd, shell=True)
