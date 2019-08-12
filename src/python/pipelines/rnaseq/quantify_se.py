import argparse
import subprocess
import concurrent.futures
from datetime import datetime
from collections import defaultdict
import json
import os
import re

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

    exc = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'pipeline_se.py')

    cmd = ' '.join([exc] + params)

    if lsf:
        cmd = bsub(cmd, 'rnaseq')
    
    return(cmd) 

def run(params, lsf=True, debug=False):

    if debug: return(cmd(params, lsf))

    subprocess.check_call(cmd(params, lsf), shell=True)	

def arg_array(paths):

    with open(paths.get('idx_map')) as fp:
        idx_map = {row.split('\t')[0]:row.split('\t')[1] for row in fp.read().strip().split('\n')} 

    fq_unordered_array = [i for i in os.listdir(paths.get('seq_dir')) if i.endswith('.fastq.gz')]
    
    fq_dict = defaultdict(list)

    def sid(x):
        return(re.match(r'[A-Za-z]+[0-9](?=_[0-9])', x).group(0))
    for fq_file in fq_unordered_array:
        srr = fq_file.split('.')[0]    
        fq_dict[sid(idx_map[srr])].append(fq_file)

    for key in fq_dict.keys():
        fq_dict[key].sort()
        fq_dict[key] = [os.path.join(paths.get('seq_dir'), i) for i in fq_dict[key]]
        
    for sample_id, fqs in fq_dict.items():
        yield [sample_id, paths.get('genome_dir'), paths.get('cdna'), paths.get('gtf'), '-o', os.path.join(args.output_dir, sample_id), '--fqs ' + ' '.join(fqs)] 
    
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
