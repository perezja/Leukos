import argparse
import subprocess
import pandas as pd
import numpy as np
from datetime import datetime
from common import explode, bsub
import re
import time
import tempfile
import shutil
import json
import os

# note, for conservation annotations, the <build>_phastConsElements60way.txt.gz 
# file has to be manually donwloaded from ucsc. Additionally, the gene-based 
# annotation has to be downloaded but can we done with ANNOVAR. Both files must
# be saved in the same database directory. 

parser = argparse.ArgumentParser(prog='Gene annotation of peak features using ANNOVAR.')
parser.add_argument('features', type=str, help='Path to peak feature file (.saf)')
parser.add_argument('db', type=str, help='Path to database directory (from "-downb" command in ANNOVAR).')
parser.add_argument('build', type=str, help='Genome build (e.g., mm10)')
parser.add_argument('prefix', type=str, help='Prefix for outfile.')
parser.add_argument('-d', '--distance', type=int, default=2000, help='Distance threshold to define upstream/downstream of a gene.')
parser.add_argument('-o', '--output_dir', default='.', help='Window size for merging intra and inter replicate peaks.')
parser.add_argument('--debug', action='store_true', help='')
args = parser.parse_args()

args.features = os.path.abspath(args.features)
args.db = os.path.abspath(args.db)
args.output_dir = os.path.abspath(args.output_dir)

def bsub(cmd, mem, gtmp, docker_image, job_name, debug=args.debug):

    mem = str(mem)
    gtmp=str(gtmp)
    bsub = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false" \
    + " bsub -Is -q research-hpc" \
    + " -J " + job_name \
    + " -M "+mem+"000000" \
    + " -R 'select[mem>"+mem+"000 && gtmp > "+gtmp+"]" \
    + " rusage[mem="+mem+"000, gtmp="+gtmp+"]'" \
    + " -a 'docker("+docker_image+")'" \
    + " /bin/bash -c '"+cmd+"'"

    if debug:
        print('* '+bsub)

    return(bsub)

def prepare_avinput(saf, outdir):

    df = pd.read_csv(saf, sep='\t', header=None, usecols=[0,1,2,3], names=['peak_id', 'chr', 'start', 'end'])

    df['null1'] = 0
    df['null2'] = 0

    df = df[['chr', 'start', 'end', 'null1', 'null2', 'peak_id']]
    
    tmp = tempfile.NamedTemporaryFile(dir=outdir, delete=False)
    df.to_csv(tmp.name, sep='\t', header=False, index=False)

    return(tmp.name)

def run_gene_annotation(avinput, dist, build, db, outdir, mem=25, gtmp=8, docker_image="apollodorus/annovar:latest"):

    prefix = os.path.split(avinput)[1].split('_')[0] 
    outfile = os.path.join(outdir, prefix + '_annot')
    
    cmd = 'annotate_variation.pl -out '+outfile \
      + ' -build ' + build + ' ' \
      + avinput \
      + ' --geneanno' \
      + ' -neargene '+str(dist) \
      + ' -dbtype refGene ' \
      + db 

    cmd = bsub(cmd, mem, gtmp, docker_image, 'annovar')
    po = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

    return((po, outfile + '.variant_function'))

def format_gene_annotation(variant_function_file):

    # ANNOVAR outputs ref and obs alleles (null cols)
    # which for peak annotation is irrelevant 

    # ['exonic', 'splicing', 'intronic', ncRNA'] gives gene name.
    # ['intergenic', 'upstream', 'downstream'] gives neighboring genes 
    #     or gene(s) hits with distances.
    # ['UTR5', 'UTR3'] gives gene name with refGene tag support in parantheses
    #     e.g., (NM_021431:c.*1592_*1724delins0).

    colnames = ['gene_region', 'gene_desc', 'chr', 'start', 'end', 'null1', 'null2', 'peak_id']

    df = pd.read_csv(variant_function_file, sep='\t', header=None, names=colnames)

    df = df.drop(['null1', 'null2'], axis=1)

    # two types of annotation cases with 
    # multiple gene hits. Format each separately.

    # first case: multiple up/downstream gene hits 
    spec_df = df[df['gene_region'] == 'upstream;downstream']
    df = df.drop(spec_df.index.tolist())

    spec_df['gene_region'] = spec_df['gene_region'].map(lambda x: x.split(';'))
    spec_df['gene_desc'] = spec_df['gene_desc'].map(lambda x: x.split(';'))
    spec_df = explode(spec_df, ['gene_desc', 'gene_region'])

    df = pd.concat([df, spec_df], ignore_index=True, sort=True)

    spec_cols = ['upstream', 'downstream']

    # second case: equidistant up/downstream gene hits
    equid_df = df[df['gene_desc'].str.contains(',') & df['gene_region'].isin(spec_cols)] 
    df = df.drop(equid_df.index.tolist()) 

    # add (dist=<int>) suffix to each ',' separated gene_desc value
    equid_df['gene_desc'] = equid_df['gene_desc'].map(lambda x: [x.split(',')[i] + re.search(r'(\(dist\=[0-9]+\))', x.split(',')[-1]).group(0) for i in range(len(x.split(',')) - 1)] + [x.split(',')[-1]]) 

    equid_df = explode(equid_df, ['gene_desc'])

    df = pd.concat([df, equid_df], ignore_index=True, sort=True)

    # make distance into separate columns 
    dist_cols = ['upstream', 'downstream', 'intergenic']
    dist_df = df[df['gene_region'].isin(dist_cols)]   
    df = df.drop(dist_df.index.tolist()) 

    dist_df['dist'] = dist_df['gene_desc'].map(lambda x: re.search(r'(?<=\(dist=)[0-9A-Z]+(?=\))', x).group(0))
    df['dist'] = 0

    df = pd.concat([df, dist_df], ignore_index=True, sort=True)

    # make gene as separate column 
    df['gene'] = df['gene_desc'].map(lambda x: re.search(r'(^[A-Za-z0-9]+)', x).group(0))

    # drop 'NONE' intergenic flanking rows
    mask = df['gene'].map(lambda x: x == 'NONE')
    df = df[[not i for i in mask]]

    # final formatting
    df = df[['gene', 'gene_region', 'chr', 'start', 'end', 'peak_id', 'dist']]
    df['dist'] = df['dist'].astype(np.int32)

    return(df) 

def regroup_gene_regions(df):

    peak_annot = df 

    peak_annot['gene_region'] = peak_annot['gene_region'].map(lambda x: 'promoter-proximal' if (x=='UTR5' or x=='upstream') else x)
    peak_annot['gene_region'] = peak_annot['gene_region'].map(lambda x: 'UTR3' if (x=='UTR3' or x=='downstream') else x)
    peak_annot['gene_region'] = peak_annot['gene_region'].map(lambda x: 'exonic' if (x=='exonic' or x=='ncRNA_exonic') else x)
    peak_annot['gene_region'] = peak_annot['gene_region'].map(lambda x: 'intronic' if (x=='intronic' or x=='splicing' or x=='ncRNA_intronic' or x=='ncRNA_splicing') else x)

    return(peak_annot)

def main():
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)

    print("[ {} ] Preparing input.".format(datetime.now().strftime("%b %d %H:%M:%S")))   
    avinput = prepare_avinput(args.features, tmpdir) 

    print("[ {} ] Running ANNOVAR.".format(datetime.now().strftime("%b %d %H:%M:%S")))
    po, fp = run_gene_annotation(avinput, args.distance, args.build, args.db, tmpdir)

    # wait for processes to finish
    while po.poll() is None:
        time.sleep(0.5) 

    print("[ {} ] Formatting Annotation.".format(datetime.now().strftime("%b %d %H:%M:%S")))    

    annot_df = format_gene_annotation(fp)
    annot_df = regroup_gene_regions(annot_df)

    fn = args.prefix + '.gene_annot.txt' 
    outfile = os.path.join(args.output_dir, fn)

    annot_df.to_csv(outfile, sep='\t', index=False)

    print("[ {} ] Done.".format(datetime.now().strftime("%b %d %H:%M:%S")))
    print("wrote to: {}".format(outfile)) 
    
    shutil.rmtree(tmpdir)
            
if __name__ == '__main__':
    main()

