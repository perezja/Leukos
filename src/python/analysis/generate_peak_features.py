import argparse
import subprocess
import pandas as pd
import numpy as np
from datetime import datetime
import shutil
import tempfile
import json
import os

parser = argparse.ArgumentParser(prog='Merge peaks across replicates and generate peak features.')
parser.add_argument('json', type=str, help='Path to input json listing .bam and .narrowPeak input files.')
parser.add_argument('-w', '--window_size', type=int, default=10, help='Window size for merging intra and inter replicate peaks.')
parser.add_argument('-o', '--output_dir', default='.', help='Window size for merging intra and inter replicate peaks.')
parser.add_argument('--debug', action='store_true', help='')
args = parser.parse_args()

def intra_rep_peak_merge(narrow_peak, tmpdir, debug=args.debug):

    # keep peak name and mean_signal (cols=4,7) 

    prefix = os.path.split(narrow_peak)[1].split('_')[0]
    outfile = os.path.join(tmpdir, prefix + 'intra_rep_merge.narrowPeak')

    cmd='tail -n+2 ' + narrow_peak \
      + ' | sort -k1,1 -k2,2n' \
      + ' | bedtools merge -i -' \
      + ' -d ' + str(args.window_size) \
      + ' -c 4,5,7' \
      + ' -o collapse,mean,mean' \
      + ' -delim ' + '","' \
      + ' > ' + outfile 

    if debug:
        print('* '+cmd)

    subprocess.check_call(cmd, shell=True)

    df = pd.read_csv(outfile, sep='\t', header=None, names=['chr', 'start', 'end', 'peak_names', 'score', 'signalValue']) 

    df['peak_names'] = df['peak_names'].apply(lambda x: prefix + '_' + '_'.join([peak_ids.split('_')[1] for peak_ids in x.split(',')])) 

    df.to_csv(outfile, sep='\t', header=False, index=False)
    return(outfile)

def inter_rep_peak_merge(beds, group, tmpdir, debug=args.debug): 

    tmp1 = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False) 
    subprocess.check_call('cat ' + ' '.join(beds) + ' > ' + tmp1.name, shell=True)

    tmp2 = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False) 

    cmd='sort -k1,1 -k2,2n ' + tmp1.name \
      + ' | bedtools merge -i -' \
      + ' -d ' + str(args.window_size) \
      + ' -c 4,1,5,6' \
      + ' -o collapse,count,mean,mean' \
      + ' -delim '+ '","' \
      + ' > ' + tmp2.name 

    if debug:
        print('* '+cmd)

    subprocess.check_call(cmd, shell=True)

    df = pd.read_csv(tmp2.name, sep='\t', header=None, names=['chr', 'start', 'end', 'peak_id', 'count', 'score','signalValue'])

    # filter for peaks with at least one replicate support

    mask = df['count'].map(lambda x: x >= 2)
    df = df[mask]

    # use narrowPeak format for visualization

    df['strand'] = '.'

    df = df[['chr', 'start', 'end', 'peak_id', 'score', 'strand', 'signalValue']]

    outfile = os.path.join(args.output_dir, group + '_consensus_peaks.narrowPeak')   

    desc = 'track type=narrowPeak visibility=3 db=mm10 name={}\n'.format(group)

    with open(outfile, 'w+') as fp:
        fp.write(desc)
        df.to_csv(fp, sep='\t', header=False, index=False)

    return(outfile)

def generate_peak_features(consensus_peak_files, tmpdir, debug=args.debug):

    tmp1 = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False) 
    subprocess.check_call('cat ' + ' '.join(consensus_peak_files) + ' > ' + tmp1.name, shell=True)

    tmp2 = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False) 

    cmd='sort -k1,1 -k2,2n ' + tmp1.name \
      + ' | bedtools merge -i -' \
      + ' -d ' + str(args.window_size) \
      + ' -c 4' \
      + ' -o collapse' \
      + ' -delim '+ '","' \
      + ' > ' + tmp2.name 

    if debug:
        print('* '+cmd)

    subprocess.check_call(cmd, shell=True)

    df = pd.read_csv(tmp2.name, sep='\t', header=None, names=['chr', 'start', 'end', 'peak_id'])

    # format as saf annotation file
    df['strand'] = '.'
    df = df[['peak_id', 'chr', 'start', 'end', 'strand']]

    saf = os.path.join(os.path.abspath(args.output_dir), 'peak_features.saf') 
    df.to_csv(saf, sep='\t', header=False, index=False)

    return(saf)

def main():
    
    with open(args.json) as fp:
        peak_data = json.load(fp)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    tmp_dir = tempfile.mkdtemp(dir=args.output_dir)
        
    rep_consensus_peak_files = list()

    for group, paths in peak_data.items():

        np_files = paths.get('narrowPeaks')

        intra_rep_merged_peak_files = list()

        print("[ {} ] Merging +/- {} bp peaks for {}.".format(datetime.now().strftime("%b %d %H:%M:%S"), args.window_size, group))

        for npf in np_files:

            rep_peak_file = intra_rep_peak_merge(npf, tmp_dir)
            intra_rep_merged_peak_files.append(rep_peak_file)

        print("[ {} ] Making consensus peak set for {}.".format(datetime.now().strftime("%b %d %H:%M:%S"), group))

        rep_consensus_peaks = inter_rep_peak_merge(intra_rep_merged_peak_files, group, tmp_dir)
        rep_consensus_peak_files.append(rep_consensus_peaks)

    print("[ {} ] Generating peak features (.saf) file.".format(datetime.now().strftime("%b %d %H:%M:%S")))

    saf = generate_peak_features(rep_consensus_peak_files, tmp_dir)

    print('wrote to: '+ saf)

if __name__ == '__main__':
    main()

