import numpy as np
import pandas as pd
import re

def get_motif_name(homer_motif_id):

    name = re.search(r'(?<=BestGuess:)[^/]*(?=/)', homer_motif_id).group(0)
    name = re.sub(r'[,(]', '-', name)

    return(re.sub(r'[)]', '', name))

class Motif():
    def __init__(self, homer_motif_file, experiment):

        self.experiment = experiment 

        with open(homer_motif_file) as fp:

            header = fp.readline().strip()

            elem, desc, log_odds_threshold, logP_enrichment, null, occurence_info = header.split('\t')

            # make bash-friendly name for file prefixes
           
            self.name = self.get_motif_name(desc)
            self.desc = desc 

            self.lod_threshold = log_odds_threshold
            self.lpval = -1 * float(logP_enrichment)

            t, b, p =  occurence_info.split(',')

            self.target_hits = float(re.search(r'(?<=T:)[0-9.]+(?!\))', t).group(0))
            self.target_freq = float(re.search(r'(?<=\()[0-9]+\.[0-9]+', t).group(0))/100
            self.bg_hits = float(re.search(r'(?<=B:)[0-9.]+(?!\))', b).group(0))
            self.bg_freq = float(re.search(r'(?<=\()[0-9]+\.[0-9]+', b).group(0))/100
            self.pval = float(re.search(r'(?<=P:)-?\d+(?:\.\d*)?(?:[eE][+\-]?\d+)?', p).group(0))

            self.motif = re.search(r'(?<=>)[ACGTRYKMSWBDHVN]+', elem).group(0)
            self.dim = (len(self.motif), 4)

            self.pwm = np.array([line.split('\t') for line in fp.read().strip().split('\n')], dtype=np.float32) 

    def get_motif_name(self, homer_motif_id):

        name = re.search(r'(?<=BestGuess:)[^/]*(?=/)', homer_motif_id).group(0)
        name = re.sub(r'[,(]', '-', name)

        return(re.sub(r'[)]', '', name))

    def row(self):

        col = pd.Series([self.name, self.desc, self.experiment, self.motif, self.dim[0], self.target_hits, self.target_freq, self.bg_hits, self.bg_freq, self.lod_threshold, self.pval], index=['name', 'desc', 'experiment', 'motif', 'length', 'target_hits', 'target_freq', 'background_hits', 'bg_freq', 'log_threshold', 'pval_enrichment']) 
        return(col.to_frame().transpose())

    def __repr__(self):
        return('motif: {}\ntarget hits: {}\nprop(targets): {}\nbackground hits: {}\nprop(bg): {}\npval(enrichment): {}'.format(self.motif, self.target_hits, self.target_freq, self.bg_hits, self.bg_freq, self.pval))
    def __str__(self):
        return('name: {} ; dim: {}'.format(self.name, self.dim))


def bsub(cmd, mem, gtmp, docker_image, job_name, debug=False):

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
        quit()

    return(bsub)

def explode(df, lst_cols, fill_value='', preserve_index=False):
    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
        and len(lst_cols) > 0
        and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values    
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
                col:np.repeat(df[col].values, lens)
                for col in idx_cols},
                index=idx)
             .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
                            for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (res.append(df.loc[lens==0, idx_cols], sort=False)
                  .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:        
        res = res.reset_index(drop=True)
    return res


