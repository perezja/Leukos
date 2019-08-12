#!/usr/bin/env python3
# Author: Francois Aguet

"""
changed default path to
.jar executable
"""
import argparse
import os
import struct
import subprocess
from datetime import datetime
import contextlib
import shutil

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


parser = argparse.ArgumentParser(description='run MarkDuplicates from Picard.')
parser.add_argument('input_bam', type=str, help='BAM file')
parser.add_argument('prefix', type=str, help='Prefix for output files; usually <sample_id>')
parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Output directory')
parser.add_argument('-m', '--memory', default='8', type=str, help='Memory, in GB')
parser.add_argument('--max_records_in_ram',
                    type=str,
                    default='150000',
                    help='Specifies the number of records stored in RAM before spilling to disk.')
parser.add_argument('--tmp_dir',
                    type=str,
                    default=None,
                    help='')
parser.add_argument('--optical_duplicate_pixel_distance', default=100, help='Maximum offset between two duplicate clusters. 100 (default) is appropriate for unpatterned, 2500 recommended for patterned flowcells.')
parser.add_argument('--jar', default='/opt/picard-tools/picard.jar', help='Path to Picard jar')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Starting MarkDuplicates', flush=True)

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
if not os.path.exists(args.tmp_dir):
    os.makedirs(args.tmp_dir)

outfn = os.path.split(args.input_bam)[1].replace('.bam', '.md.bam')
outfile = os.path.join(args.output_dir, outfn) 

subprocess.check_call('java -jar -Xmx'+args.memory+'g '+args.jar\
    +' MarkDuplicates I='+args.input_bam\
    +' O='+outfile\
    +' PROGRAM_RECORD_ID=null' \
    +' M='+os.path.join(args.output_dir, args.prefix+'.marked_dup_metrics.txt') \
    +' TMP_DIR='+args.tmp_dir \
    +' ASSUME_SORT_ORDER=coordinate' \
    +' MAX_RECORDS_IN_RAM='+args.max_records_in_ram \
    +' OPTICAL_DUPLICATE_PIXEL_DISTANCE='+str(args.optical_duplicate_pixel_distance), 
shell=True)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished MarkDuplicates', flush=True)
