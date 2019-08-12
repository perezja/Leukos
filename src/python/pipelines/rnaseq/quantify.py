import argparse
import subprocess
import os
from datetime import datetime
import concurrent.futures
import json

def bsub(cmd, job_name, queue="research-hpc", mem=50, gtmp=15, docker_image="apollodorus/brain-eqtl:rnaseq2"):

    mem = str(mem)
    gtmp = str(gtmp)
    bsub = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false" \
        + " bsub -Is -q "+queue \
        + " -J "+job_name \
        + " -M "+mem+"000000" \
        + " -R 'select[mem>"+mem+"000 && gtmp > "+gtmp+"]" \
        + " rusage[mem="+mem+"000, gtmp="+gtmp+"]'" \
        + " -a 'docker("+docker_image+")'" \
        + " /bin/bash -c '"+cmd+"'" 

    return(bsub)

def cmd(params, lsf):
    """Sets pipeline command for rnaseq sample."""

    exc = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'pipeline2.py')

    cmd = ' '.join([exc] + params)

    if lsf:
        cmd = bsub(cmd, 'rna.seq')
    
    return(cmd) 

def run(params, lsf=True, debug=False):

    if debug: return(cmd(params, lsf=True))

    subprocess.check_call(cmd(params, lsf), shell=True)	

def arg_array(paths):

    with open(paths.get('rna.fqs')) as fp:
        fqs_dict = json.load(fp) 

    for sample_id, pair_mates in fqs_dict.items():

        fq1 = [pair[0] for pair in pair_mates]
        fq2 = [pair[1] for pair in pair_mates]

        fq1 = ' '.join(fq1)
        fq2 = ' '.join(fq2)

        yield [sample_id, paths.get('rna.genomeDir'), paths.get('rna.cdna'), paths.get('rna.gtf'), '--fq1', fq1, '--fq2', fq2, '-o', os.path.join(args.output_dir, sample_id)] 
    
def debug(paths):
    for params in arg_array(paths):
        print(run(params, debug=True))
        quit() 

parser = argparse.ArgumentParser(description='Run RNAseq pipeline.')
parser.add_argument('json', help='JSON file containing paths to expression, covariate and genotype directories.')
parser.add_argument('--max_processes', type=int, default=100, help='Suffix for sequence files.')
parser.add_argument('--lsf', action='store_true', help='Flag specifying whether to run in parallel on LSF.')
parser.add_argument('--debug', action='store_true', help='Flag specifying whether to print commands.')
parser.add_argument('-o', '--output_dir', default='./', help='Directory for output files.')

args = parser.parse_args()

def main():

    with open(args.json) as fp:
        paths = json.load(fp)

    if args.debug:
        debug(paths)

    executor = concurrent.futures.ProcessPoolExecutor(max_workers=args.max_processes)
    futures = [executor.submit(run, args, lsf=True) for args in arg_array(paths)] 
    concurrent.futures.wait(futures)

    for f in futures:
        if f.result():
            print(f.result())

if __name__ == "__main__":
    main()
