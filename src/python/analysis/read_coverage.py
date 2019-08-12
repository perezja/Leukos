import argparse
import subprocess
from datetime import datetime
from common import bsub
import time
import tempfile
import shutil
import json
import os

def run_coverage(bam, ref_set, outfile, tmpdir, sort=False, mem=25, gtmp=4, docker_image="apollodorus/bioinf:v1"):  

    tmp = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)

    sort = '#!/bin/bash\nsamtools sort {} > {}'.format(bam, tmp.name)
    coverage = 'samtools view -b {} | bedtools coverage -a {} -b - -d > {}'.format(tmp.name, ref_set, outfile) 

    exc = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)

    cmd = [coverage]
    if sort:
        cmd = [sort, coverage]

    with open(exc.name, 'w+') as fp:
        fp.write('\n'.join(cmd))

    os.chmod(exc.name, 0o755)

    cmd = bsub(exc.name, mem, gtmp, docker_image, 'genomeCoverage')
    po = subprocess.Popen(cmd, shell=True)

    return(po)

def parse_args():

    parser = argparse.ArgumentParser(prog='Make genome-wide accessibility coverage files.')
    parser.add_argument('json', type=str, help='Path to input json listing .bam and .narrowPeak input files.')
    parser.add_argument('ref_set', type=str, help='BED file with peak set to find coverage against.')
    parser.add_argument('-o', '--output_dir', default='.', help='')
    parser.add_argument('--debug', action='store_true', help='')

    return(parser.parse_args())

def main():

    args = parse_args()

    with open(args.json) as fp:
        peak_data = json.load(fp)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)

    processes = list()
    for group, paths in peak_data.items():

        bam_files = paths.get('bams')

        for bam in bam_files:
        
            outfile = os.path.split(bam)[1].split('_')[0] + '_coverage.bed'
            outfile = os.path.join(args.output_dir, outfile)
            processes.append(run_coverage(bam, args.ref_set, outfile, tmpdir))

    print("[ {} ] Running base-pair coverage over BED reference.".format(datetime.now().strftime("%b %d %H:%M:%S")))

    while processes:
        po = processes.pop()
   
        while po.poll() is None:
            time.sleep(0.5)

    print("[ {} ] Done.".format(datetime.now().strftime("%b %d %H:%M:%S")))
    print("wrote to: {}".format(args.output_dir))

    shutil.rmtree(tmpdir)
            
if __name__ == '__main__':
    main()

