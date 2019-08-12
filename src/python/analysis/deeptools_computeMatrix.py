import argparse
import subprocess
from datetime import datetime
from common import bsub
import time
import tempfile
import json
import os

def run_bamCovereage(bam, outfile, mem=25, gtmp=4, docker_image="quay.io/bgruening/galaxy-deeptools"):  

    tmp = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)

    cmd = 'deeptools bamCoverage -b {} -o {}'.format(bam, outfile) 

    cmd = bsub(cmd, mem, gtmp, docker_image, 'bamCoverage')
    po = subprocess.Popen(cmd, shell=True)

    return(po)

def parse_args():

    parser = argparse.ArgumentParser(prog='Make genome-wide accessibility coverage files.')
    parser.add_argument('json', type=str, help='Path to input json listing .bam and .narrowPeak input files.')
    parser.add_argument('-o', '--output_dir', default='.', help='')
    parser.add_argument('--debug', action='store_true', help='')

    return(parser.parse_args())

def main():

    args = parse_args()

    with open(args.json) as fp:
        peak_data = json.load(fp)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    processes = list()
    def sid(fp):
        return(os.path.split(fp)[1].split('_')[0]

    for group, paths in peak_data.items():
        
        bam_files = paths.get('bams')

        for bam in bam_files:
            outfile = os.path.join(args.output_dir, sid(bam)+'_coverage.bigwig') 
            processes.append(run_bamCoverage(bam, outfile))

    print("[ {} ] Creating bigwigs from input BAMs.".format(datetime.now().strftime("%b %d %H:%M:%S")))

    while processes:
        po = processes.pop()
   
        while po.poll() is None:
            time.sleep(0.5)

    print("[ {} ] Done.".format(datetime.now().strftime("%b %d %H:%M:%S")))
    print("wrote to: {}".format(args.output_dir))

if __name__ == '__main__':
    main()

