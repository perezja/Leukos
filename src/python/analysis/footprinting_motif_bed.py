import argparse
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from collections import defaultdict, OrderedDict
from common import get_motif_name
import pandas as pd
import numpy as np
import pickle
import warnings
import re
import json
import os

parser = argparse.ArgumentParser(prog='')
parser.add_argument('peak_bed', type=str, help='BED file with peaks of interest.')
parser.add_argument('peak_motifs', type=str, help='Output file from HOMER findMotifsGenome.pl ... "-find $MOTIF_FILE"')
parser.add_argument('fasta', type=str, help='Pickled FASTA file')
parser.add_argument('-p', '--peak_ids', type=str, help='List of peak IDs to subset BED.')
parser.add_argument('-f','--flank_size', type=int , default=100, help='Length to extend query peak regions')
parser.add_argument('-w','--window_size', type=int , default=40, help='Length to extend motif center by in BED file')
parser.add_argument('-o', '--output_dir', default='.', help='')
args = parser.parse_args()

def peak_dictionary(fp, peak_ids = None):

    if peak_ids:
        with open(peak_ids) as peak_list:
            peak_list = peak_list.read().strip().split('\n') 

    if os.path.splitext(fp)[1]=='.saf':
        header = ['peak_id', 'chr', 'start', 'end'] 
    if os.path.splitext(fp)[1]=='.bed':
        header = ['chr', 'start', 'end', 'peak_id']

    peak_df = pd.read_csv(fp, sep='\t', usecols=list(range(0,4)) ,names=header, index_col='peak_id')

    if peak_ids:
        peak_df = peak_df[peak_df.index.isin(peak_list)]

    return(peak_df.to_dict(orient='index'))

# {(peak_id, motif_id): 'seq', 'offset', 'strand', 'score'} 
def motif_peak_dictionary(fp, peak_ids):

    motif_df = pd.read_csv(fp, sep='\t', skiprows=1, names=['peak_id', 'offset', 'seq', 'homer_id', 'strand', 'score'])
    motif_df = motif_df[motif_df['peak_id'].isin(peak_ids)] 

    motif_df['motif_id'] = motif_df['homer_id'].apply(lambda x: get_motif_name(x)) + ',' + [str(i+1) for i in range(motif_df.shape[0])] 

    motif_df = motif_df.drop('homer_id', axis=1)

    # sorting by offset is critical for motif finding
    # allows tracking along the peak sequence to avoid
    # double counting of motif hits 

    motif_df.sort_values(['peak_id', 'offset'], ascending=True, inplace=True)

    motif_df.set_index('motif_id', inplace=True)

    return(motif_df.to_dict(orient='index', into=OrderedDict()))

def motif_coordinate_dictionary(motifs, genome_coordinates, peak_location, flank_size):

    motif_loc = defaultdict(dict) 

    last_index = 0
    last_pid = None 
    last_offset = None

    for motif_id, motif_instance in motifs.items():
    
            peak_id = motif_instance['peak_id']
            loc = peak_location[peak_id]
    
            loc_start = loc['start'] - flank_size
            loc_end = loc['end'] + flank_size
    
            ps = genome_coordinates[loc['chr']][loc_start:loc_end].seq 
            ms = motif_instance['seq'] 
    
            if last_pid != peak_id:
                last_index = 0
                last_offset = None
    
            mh = re.search(ms, ps.__str__()[last_index:])
    
            # repeat region motif hit
            masked_region = False
    
            if not mh:

                if re.search(ms, ps.__str__()[last_index:], re.IGNORECASE):
                    mh = re.search(ms, ps.__str__()[last_index:], re.IGNORECASE)
                    masked_region = True
                elif (last_pid == peak_id) and (motif_instance['offset'] == last_offset): 
                    warnings.warn('"{}" in "{}" has identical position match of another motif. Skipping...'.format(motif_id, peak_id)) 
                    continue
                else:
                    print('sequence "{}" for motif "{}" is not within peak region\ntry increasing --flank_size...'.format(ms, motif_id))
                    continue
    
            # Motif coordinates in 1-based index positions
            chrom = peak_location[peak_id]['chr']
            
            if last_pid == peak_id:
                start = (loc_start+1) + last_index + mh.span()[0]
                end = loc_start + last_index + mh.span()[1]
    
            else:
                start = (loc_start+1) + mh.span()[0]
                end = loc_start + mh.span()[1]
            
            strand = motif_instance['strand']
    
            # record offset from peak center
            pk_length = loc['end'] - loc['start']
            if pk_length % 2 == 0:  # even 
                pk_center = loc['start'] + (pk_length/2) 
            else:
                pk_center = loc['start'] + (pk_length/2)+1  # odd
    
            offset = start - pk_center

            # make sure all motifs in a peak unique    
            #print('{}, {}, {}'.format(peak_id, motif_id, last_pid==peak_id))
            try:
                assert start not in [motif_loc[peak_id][mid]['start'] for mid in list(motif_loc[peak_id].keys()) if peak_id==last_pid]
            except:
                with open('motif_error.tmp', 'w+') as fp:
                    json.dump(motif_loc, fp, indent=4)
                raise Exception('sequence "{}" already matched at "{}-{}"'.format(mh.group(0), start, end))
    
            motif_loc[peak_id][motif_id] = {'chr':chrom, 'start':start, 'end':end, 'strand':strand, 'seq':ms, 'offset':offset, 'masked_region':masked_region}

            # keep track of last called motif in a peak, and
            # index for non-redundant query search space 
            last_pid = peak_id
            last_index += (mh.span()[0] + 1) 
            last_offset = motif_instance['offset']

    return(motif_loc)
     
def motif_bed_file(motif_coord, window_size, prefix):

    motif_coord_df = pd.DataFrame.from_dict({(peak_id, motif_id):motif_coord[peak_id][motif_id] for peak_id in motif_coord.keys() \
                                                for motif_id in motif_coord[peak_id].keys()}, orient='index') 
    motif_coord_df.index.names = ['peak_id', 'motif_id']
    motif_coord_df.reset_index(inplace=True)

    motif_coord_df = motif_coord_df[['chr', 'start', 'end', 'peak_id', 'motif_id', 'strand', 'offset', 'seq', 'masked_region']]
    motif_coord_df.to_csv(prefix+'.coord.txt', sep='\t', index=False)
    
    motif_len = len(motif_coord_df.at[0, 'seq'])

    if (motif_len % 2) != 0:
        raise ValueError('motif length ({}) is not even'.format(motif_len))

    def construct_window(row):

        # 0-start, half-open coordinates for bedtools
        # motif center is taken as the ((motif_length/2) - 1) index position
        motif_center_coord = row['start'] + ((motif_len/2)) 

        window_start = int(motif_center_coord - (window_size + 1)) 
        window_end = int(motif_center_coord + window_size + 1) 

        return((window_start, window_end))

    motif_coord_df[['start', 'end']] = motif_coord_df.apply(construct_window, axis=1, result_type='expand')
    motif_coord_df = motif_coord_df[['chr', 'start', 'end', 'motif_id']]
    motif_coord_df.to_csv(prefix + '.coord.bed', sep='\t', index=False, header=False)

def main():
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    peak_loc = peak_dictionary(args.peak_bed, args.peak_ids)
    motifs = motif_peak_dictionary(args.peak_motifs, list(peak_loc.keys()))

    print("reading genome...")
    with open(args.fasta, 'rb') as fp:
        genome_coord = pickle.load(fp) # SeqIO::SeqRecord() object 
    print('done reading.')

    motif_loc = motif_coordinate_dictionary(motifs, genome_coordinates=genome_coord, peak_location=peak_loc, flank_size=args.flank_size) 

    prefix = os.path.join(args.output_dir, list(motifs.keys())[0].split(',')[0])
    motif_bed_file(motif_loc, window_size=args.window_size, prefix=prefix)
        
    print('wrote {}*'.format(prefix))

if __name__ == '__main__':
    main()

