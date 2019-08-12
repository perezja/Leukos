import argparse
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser(prog='Callable from peakheatmap_summit_coverage.py')
parser.add_argument('peak_covereage_file', type=str, help='BED file with peaks of interest.')
parser.add_argument('-w','--window_size', type=int , default=1000, help='Length which to extend peak summit.') 
parser.add_argument('-o', '--outfile', default='out.bed', help='')
args = parser.parse_args()

def get_peak_summit_df(peak_coverage_file, outfile):

    df = pd.read_csv(peak_coverage_file, sep='\t', header=None, names=['chr', 'start', 'end', 'peak_id', '1-pos', 'reads'])
    idx = df.groupby('peak_id')['reads'].transform(max) == df['reads']
    summit_pos_groups = df[idx].groupby('peak_id')

    summit_pos_dict = defaultdict(lambda: defaultdict())
    for peak_id, group in summit_pos_groups:

        if len(group) == 1:
            summit_row = group.iloc[0] 
        elif (len(group) % 2) == 0:
            summit_row = group.iloc[int(len(group)/2) - 1]
        else: 
            summit_row = group.iloc[int(len(group)/2)]

        # 0-start 1-end bed
        summit_pos_dict[peak_id] = {'summit_pos': (summit_row['start'] - 1) + summit_row['1-pos'], 
                                    'chr' : summit_row['chr']} 

        summit_pos_dict[peak_id]['start'] = summit_pos_dict[peak_id]['summit_pos'] - args.window_size 
        summit_pos_dict[peak_id]['end'] = summit_pos_dict[peak_id]['summit_pos'] + args.window_size 

    peak_summit_df = pd.DataFrame.from_dict(summit_pos_dict, orient='index')
    peak_summit_df.index.name = 'peak_id'
    peak_summit_df.reset_index(inplace=True)
    peak_summit_df = peak_summit_df[['chr', 'start', 'end', 'peak_id']]

    peak_summit_df.to_csv(outfile, sep='\t', index=False, header=False)
 
def main():

    get_peak_summit_df(args.peak_covereage_file, args.outfile)

if __name__ == "__main__":
    main()
