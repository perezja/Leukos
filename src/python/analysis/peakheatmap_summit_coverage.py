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
parser.add_argument('atac_bam_dir', type=str, help='Directory with atac bams.')
parser.add_argument('-p', '--peak_ids', type=str, help='List of peak IDs to subset BED.')
parser.add_argument('-o', '--output_dir', default='.', help='')
parser.add_argument('-d', '--debug', action='store_true', help='')
args = parser.parse_args()

args.peak_bed = os.path.abspath(args.peak_bed)
args.atac_bam_dir = os.path.abspath(args.atac_bam_dir)
args.output_dir = os.path.abspath(args.output_dir)

def get_peak_coverage(bam, refset, outfile, tmpdir):
    
    po = run_coverage(bam, refset, outfile, tmpdir=tmpdir)
    while po.poll() is None:
        time.sleep(0.5)

def peak_summit_routine(bam, refset, sample_id, tmpdir): 

    peak_coverage_outfile = os.path.join(tmpdir, sample_id + '_peakCoverage.bed')

    # returns after bedtools coverage complete 
    get_peak_coverage(bam, refset, peak_coverage_outfile, tmpdir=tmpdir) 
    
    exc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'peakheatmap_summit_bed.py')

    peak_summit_outfile = os.path.join(args.output_dir, sample_id+".summitCoverage.bed")

    cmd = 'python3 {} {} -o {}'.format(exc, peak_coverage_outfile, peak_summit_outfile) 
    cmd = bsub(cmd, mem=30, gtmp=4, docker_image='apollodorus/bioinf:pr', job_name='summitBed')

    subprocess.check_call(cmd, shell=True)

#    peak_summit_outfile = os.path.join(args.output_dir, sample_id+".summitCoverage.bed")

#    get_peak_coverage(bam, tmp.name, peak_summit_outfile, tmpdir) 

    print('Finished peak summit BED for "{}"'.format(sample_id))

def main():

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    bams = [fn for fn in os.listdir(args.atac_bam_dir) if fn.endswith('.bam')]

    def sid(fn):
        return(fn.split('.')[0].split('_')[0])

    read_dict = {sid(bf): os.path.join(os.path.abspath(args.atac_bam_dir), bf) for bf in bams}

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)

    # parallelization

    # step 1. get coverage over peaks
    launch = concurrent.futures.ProcessPoolExecutor(max_workers=20)
    processes = [launch.submit(peak_summit_routine, bam, args.peak_bed, sample_id, tmpdir) for sample_id, bam in read_dict.items()]
    concurrent.futures.wait(processes)

    for po in processes:
        if po.result():
            print(po.result())

    shutil.rmtree(tmpdir)

    print('Done.')

if __name__ == "__main__":
    main() 
