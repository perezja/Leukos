import argparse
import os
import re
import json
import subprocess

def format_input_json(replicates):

    input_json = dict()

    input_json["atac.genome_tsv"] = os.path.abspath(args.genome_tsv) 
    input_json["atac.paired_end"] = True 

    with open(args.sample_adapter_lookup) as fp:
        adapter_dict = {sa.split('\t')[0] : sa.split('\t')[1] for sa in fp.read().strip().split('\n')}
   
    fastq_key_base = "atac.fastqs_rep" 
    adapter_key_base = "atac.adapters_rep"

    for n, rep in enumerate(replicates):
        for k in range(2):
            fq_replicate_key = fastq_key_base + str(n+1) + '_R' + str(k+1) 
            adapter_key = adapter_key_base + str(n+1) + '_R' + str(k+1) 

            fq_dir = os.path.abspath(args.fastq_dir)
            input_json[fq_replicate_key] = [os.path.join(fq_dir, rep[k])] 
            input_json[adapter_key] = [adapter_dict[os.path.splitext(rep[k])[0]]]

    sample_id = prefix(rep[k])
   
    input_json["atac.title"] = "{} (paired-end)".format(sample_id)

    input_json["atac.description"] = "{}".format(sample_id)
    if args.description:
        input_json["atac.description"] = "{}".format(args.description)

    input_json["atac.pipeline_type"] = "atac"

    input_json["atac.auto_detect_adapter"] = False

    input_json["atac.align_only"] = False
    input_json["atac.multimapping"] = 4
    input_json["atac.enable_xcor"] = False
    input_json["atac.enable_idr"] = False
    input_json["atac.idr_thresh"] = 0.05
    input_json["atac.disable_ataqc"] = False 
    
    return(sample_id, input_json)

def prefix(item):
    match = re.match(r'^[a-zA-Z0-9]*', item)
    # last character in prefix is replicate number
    return(match.group(0)[:-1])


def get_sample_set(idx, s):

    for k in range(idx+1, len(s)):
        if prefix(s[k][0]) == prefix(s[idx][0]):
            continue
        else:
            return((k, s[idx:k]))  

    return((len(s), s[idx:len(s)])) 


parser = argparse.ArgumentParser(description='Run ENCODE ATAC-seq pipeline using Cromwell.')
parser.add_argument('fastq_dir', help='List of FASTQ replicate prefix names (e.g., <rep1 rep2 rep3> will read <rep1_R1 rep1_R2 rep2_R1 rep2_R2 rep3_R1 rep3_R2> as inputs (case-insensitive).')
parser.add_argument('genome_tsv', help='Genome data TSV file (https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/genome/download_genome_data.sh.)')
parser.add_argument('sample_adapter_lookup', default=None, help='Lookup table linking sample to adaptor sequence.')
parser.add_argument('-d', '--description', default=None, help='')
parser.add_argument('-o', '--output_dir', default='.', help='Directory for output files.')
args = parser.parse_args()

def main():

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    fqs = [i for i in os.listdir(args.fastq_dir) if (i.endswith('.fq') or i.endswith('.fastq') or i.endswith('.gz'))] 
    # sort by fn prefix
    fqs.sort(key=lambda x: os.path.splitext(os.path.split(x)[1])[0])

    mates = [(mate1, mate2) for n, mate2 in enumerate(fqs) for mate1 in fqs[:n] if mate1.split('_')[0] == mate2.split('_')[0]] 
    idx = 0
    sample_set = list()

    while (idx < (len(mates)-1)):
        idx, samples = get_sample_set(idx, mates)
        sample_set.append(samples)
         
    # sample = list of mate pair replicates
    for sample in sample_set:
        sample_id, input_json = format_input_json(sample) 

        outfile = os.path.join(args.output_dir, sample_id+'_atac.json')
        with open(outfile, 'w+') as fp:
            json.dump(input_json, fp, sort_keys=True, indent=4, separators=(',', ': '))
        print('wrote to: {}'.format(outfile))
 
if __name__ == "__main__":
    main()
