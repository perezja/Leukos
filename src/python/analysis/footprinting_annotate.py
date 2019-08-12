import argparse
import subprocess
import pandas as pd
import numpy as np
from datetime import datetime
from common import bsub 
import tempfile
import shutil
import json
import re
import os

parser = argparse.ArgumentParser(prog='Annotate dnase2tf output with peaks.')
parser.add_argument('json', type=str, help='JSON with group to footprint paths.')
parser.add_argument('peak_bed', type=str, help='Peak BED file.')
parser.add_argument('-w', '--window_size', type=int, default=150, help='Window size for extending center of footprints.')
parser.add_argument('-o', '--output_dir', default='.', help='Window size for merging intra and inter replicate peaks.')
parser.add_argument('--debug', action='store_true', help='')
args = parser.parse_args()

def footprint_peak_intersect(footprint_bed, peak_bed):

    # keep peak name and mean_signal (cols=4,7) 
    footprint_bed = os.path.abspath(footprint_bed)
    peak_bed = os.path.abspath(peak_bed)
    
    prefix = os.path.split(footprint_bed)[1].split('_')[0]
    outfile = os.path.join(args.output_dir, prefix + '.peakFootprint.bed')

    print('intersecting footprints with peaks...')
    cmd='bedtools intersect -a {} -b {} > {}'.format(peak_bed, footprint_bed, outfile)

    subprocess.check_call(cmd, shell=True)

    return(outfile)

def footprint_window_bed(footprint_fp, window_size, prefix):

    print('constructing BED window...')
    footprint_bed = pd.read_csv(footprint_fp, sep='\t', header=None, names=['chr', 'start', 'end', 'peak_id'])
    homer_bed = pd.read_csv(footprint_fp, sep='\t', header=None, names=['chr', 'start', 'end', 'peak_id'])

    print(footprint_bed.head())

    def construct_window(row):

        # 0-start, half-open coordinates for bedtools
        # motif center is taken as the ((motif_length/2) - 1) index position
        fp_len = row['end'] - row['start'] + 1 
        fp_center_coord = row['start'] + (fp_len/2) 

        window_start = int(fp_center_coord - (window_size + 1)) 
        window_end = int(fp_center_coord + window_size + 1) 

        return((window_start, window_end))

    footprint_bed[['start', 'end']] = footprint_bed.apply(construct_window, axis=1, result_type='expand')
    footprint_bed.sort_values(['chr', 'start'], inplace=True)

    outfile = os.path.join(args.output_dir, prefix+'.peakFootprint_window'+str(window_size)+'bp.txt')
    footprint_bed.to_csv(outfile, sep='\t', index=False, header=False)

    # specialized window for HOMER
    window_size = 50

    homer_bed[['start', 'end']] = homer_bed.apply(construct_window, axis=1, result_type='expand') 
    homer_bed['strand'] = '.'
    homer_bed = homer_bed[['peak_id', 'chr', 'start', 'end', 'strand']]

    outfile = os.path.join(args.output_dir,prefix+'.peakFootprint.homer.txt')

    homer_bed.to_csv(outfile, sep='\t', index=False, header=False)
    
    return()

def main():
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.json) as fp:
        paths = json.load(fp)

    for group, bed in paths.items():

        inter_peak_fp = footprint_peak_intersect(bed, args.peak_bed)
        footprint_window_bed(inter_peak_fp, args.window_size, group)
        
    print("wrote to: {}".format(args.output_dir)) 
    
if __name__ == '__main__':
    main()

