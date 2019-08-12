import argparse
import subprocess
from datetime import datetime
from collections import defaultdict
import re
import glob
import time
import tempfile
import shutil
import json
import os

def parse_args():

    parser = argparse.ArgumentParser(prog='Make genome-wide accessibility coverage files.')
    parser.add_argument('json', type=str, help='Path to input json listing .bam and .narrowPeak input files.')
    parser.add_argument('--paired', action='store_true', help='')
    parser.add_argument('-o', '--output_dir', default='.', help='')

    return(parser.parse_args())

def main():

    args = parse_args()

    with open(args.json) as fp:
        peak_data = json.load(fp)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    group_chr_files = {} 

    def chr_name(x):

        return(x.split('.nodup_')[1].replace('.txt',''))

    for group, paths in peak_data.items():

        bam_files = paths.get('bams')

        chr_files = defaultdict(list) 

        sample_bamsplit_dirs = list()
        for bam in bam_files:

            bam_fn = os.path.split(bam)[1]
            if args.paired:
                bamsplit_outdir = os.path.join(os.path.split(bam)[0], bam_fn.replace('.bam',''))
            else:
                bamsplit_outdir = os.path.join(os.path.split(bam)[0])

            sample_bamsplit_dirs.append(bamsplit_outdir)
            split_compact_prefix = os.path.join(bamsplit_outdir, bam_fn.replace('.bam', '_chr*'))

            if not glob.glob(split_compact_prefix):
                print('splitting bam for {}'.format(bam_fn.split('_')[0])) 

                if paired: 
                    subprocess.check_call('bam2compact.sh {}'.format(bam), shell=True)
                else:
                    subprocess.check_call('split.pl {}'.format(bam.replace('.bam','_compact.txt')), shell=True)
                    subprocess.check_call('paired_end_bam2split.r {}'.format(bam), shell=True)

        for sample_dir in sample_bamsplit_dirs:
            for fp in os.listdir(sample_dir): 
                if (re.search(r'_chr.*\.txt',fp) and re.match(group,fp)):
                    chr_files[chr_name(fp)].append(os.path.join(sample_dir,fp))

        group_chr_files[group] = chr_files

    # collect split _chr* files into main groups

    outfile = os.path.join(args.output_dir, 'groupings.json')
    with open(outfile,'w+') as fp:
        json.dump(group_chr_files, fp, indent=4)

    print('concatenating replicates...')
    for group, chr_dict in group_chr_files.items():

        outdir = os.path.join(args.output_dir, group)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        for chrom_label, replicates in chr_dict.items():
            outfile = os.path.join(args.output_dir, group, group+'_'+chrom_label+'.txt')
            subprocess.check_call('cat '+' '.join(replicates)+' > '+outfile, shell=True)
    
    print('done.')
        
if __name__ == '__main__':
    main()

