from subprocess import PIPE, Popen
import random as rd
import tqdm
import pickle
import numpy as np
import multiprocessing
import modules
from Bio import SeqIO
from tqdm.contrib.concurrent import process_map
import pandas as pd
from datetime import datetime
import itertools
import RNA

bases = ['A','G','C','T']
safe = ['A','G','C','T', '-', '.']
blanks = ['-','.']
species = ['Cyanobacteria']

def process_val(raw_input):
    val = []
    s=0
    for char in raw_input:
        if char=='-':
            s=1
        if s:
            if char!=')':
                val.append(char)
            else:
                s=0
    val = float(''.join(val))
    return val

def mutate(pxn, full_seq):
    new_s = full_seq.copy()

    new_l = bases.copy()
    new_l.remove(new_s[pxn])

    rand = rd.randint(0,2)
    new_s[pxn] = new_l[rand]

    return new_s

def MFE_func_all(s, pxns):
    ref_input = ''.join(s)
    (nonsense, MFE_val) = RNA.fold(ref_input)
    ref_val = MFE_val

    MFE_list = []
    for x in pxns:
        if x=='BLANK':
            MFE_list.append(np.nan)
        else:
            inp = ''.join(mutate(x,s))
            (nonsense2, MFE_val2) = RNA.fold(inp)
            value_for_dict = MFE_val2
            MFE_list.append(value_for_dict-ref_val)

    return MFE_list, ref_val

def new_MFE_func(s, pairs):
    ref_input = ''.join(s)
    ref_val = RNA.fold(ref_input)[1]

    MFE_list = []
    for x in pairs:
        if x=='BLANK':
            MFE_list.append(np.nan)
        else:
            inp1 = ''.join(mutate(x[0], s))
            inp2 = ''.join(mutate(x[1], s))
            val1 = RNA.fold(inp)[1]
            val2 = RNA.fold(inp)[1]
            value_for_dict = (abs(val1-ref_val) + abs(val2-ref_val)) / 2
            MFE_list.append(value_for_dict)

    return MFE_list, ref_val

def run_sequence(seq, sites):
    site_nums = []
    pure_seq_l = []

    for i in range(len(seq)):
        if seq[i] not in blanks:
            site_nums.append(i)
            pure_seq_l.append(seq[i])

    sites_2 = []
    for site in sites:
        try:
            sites_2.append(site_nums.index(site))
        except:
            sites_2.append('BLANK')

    MFE_vals, rv = MFE_func_all(pure_seq_l, sites_2)

    return MFE_vals, rv

def new_run_sequence(seq, sites):
    pairs = list(itertools.combinations(sites,2))

    site_nums = []
    pure_seq_l = []

    for i in range(len(seq)):
        if seq[i] not in blanks:
            site_nums.append(i)
            pure_seq_l.append(seq[i])

    pairs_2 = []
    for pair in pairs:
        try:
            pairs_2.append((site_nums.index(pair[0]), site_nums.index(pair[1])))
        except:
            pairs_2.append('BLANK')

    MFE_vals, rv = new_MFE_func(pure_seq_l, pairs_2)
    return MFE_vals, rv

if __name__=='__main__':
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print("STARTING:", date_time)

    b = []
    for spec in species:
        with open(f'{spec}Results_E_0.232.pickle', 'rb') as handle:
            b.append(pickle.load(handle)[0.95])

        sites = list(modules.flatten(b))

        sites=sites
        pairs = list(itertools.combinations(sites,2))

        with open(f'data/SILVA_138.1_{spec}.fasta') as handle:
            data = list(SeqIO.parse(handle, "fasta"))

        bad_seqs = []
        for sequence_ind in range(100):
            s=0
            for char in str(data[sequence_ind].seq):
                if char not in safe:
                    s=1
            if s:
                bad_seqs.append(sequence_ind)

        results_l = []
        index_col = []

        pool_input = []
        c=0
        j=0
        while c<20:
            if j not in bad_seqs:
                c += 1
                index_col.append(j)
                pool_input.append((str(data[j].seq), sites))
            j += 1

        with multiprocessing.Pool() as pool:
            pool_output = pool.starmap(new_run_sequence, pool_input)

        formatted_output = [[y]+x for x,y in pool_output]
        columns_labels = ['REFERENCE']+[str(_) for _ in pairs]
        results_df = pd.DataFrame(formatted_output, index=index_col, columns=columns_labels)

        results_df.to_excel(f'New_{spec}_MFE50.xlsx')
        print('____________________________________________________________________________')
        print(f'{spec} IS FINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        now = datetime.now()
        date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print("End Time: ", date_time)
        print('____________________________________________________________________________')
