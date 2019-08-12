import argparse
import subprocess
from datetime import datetime
from common import bsub
import time
import json
import os

def run_bamCoverage(bam, outfile, mem=20, gtmp=4, docker_image="apollodorus/bioinf:pr"):  

    cmd = 'bamCoverage --normalizeUsing CPM --binSize 10 -b {} -o {}'.format(bam, outfile) 

    cmd = bsub(cmd, mem, gtmp, docker_image, 'genomeCoverage')

    po = subprocess.Popen(cmd, shell=True)

    return(po)

def parse_args():

    parser = argparse.ArgumentParser(prog='Make genome-wide accessibility bigwig coverage files.')
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
    for group, paths in peak_data.items():

        bam_files = paths.get('bams')

        for bam in bam_files:
        
            outfile = os.path.join(args.output_dir, os.path.split(bam)[1].split('_')[0] + "_coverage.bw")
            processes.append(run_bamCoverage(bam, outfile))

    print("[ {} ] Running deeptools bamCoverage.".format(datetime.now().strftime("%b %d %H:%M:%S")))

    while processes:
        po = processes.pop()
   
        while po.poll() is None:
            time.sleep(0.5)

    print("[ {} ] Done.".format(datetime.now().strftime("%b %d %H:%M:%S")))
    print("wrote to: {}".format(args.output_dir))

if __name__ == '__main__':
    main()

