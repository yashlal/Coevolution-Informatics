from subprocess import PIPE, Popen
import random as rd
import tqdm
import p_tqdm
import pickle
import numpy as np
import multiprocessing
import modules
from Bio import SeqIO
from tqdm.contrib.concurrent import process_map
import pandas as pd
from datetime import datetime
import RNA

bases = ['A','G','C','T']
safe = ['A','G','C','T', '-', '.']
blanks = ['-','.']
species = ['Bacteroidota']

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
		(nonsense2, MFE_val2) = RNA.fold(ref_input)
		value_for_dict = MFE_val
		MFE_list.append(value_for_dict-ref_val)

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

if __name__=='__main__':
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print("STARTING:", date_time)

    b = []
    for spec in species:
        with open(f'Results\{spec}Results_E_0.116.pickle', 'rb') as handle:
            b.append(pickle.load(handle)[0.95])

        sites = list(modules.flatten(b))

        sites=sites

        with open(f'data\SILVA_138.1_{spec}.fasta') as handle:
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

        for j in range(50):
            if j not in bad_seqs:
                index_col.append(j)
                pool_input.append((str(data[j].seq), sites))

        with multiprocessing.Pool() as pool:
            pool_output = pool.starmap(run_sequence, pool_input)
        formatted_output = [[y]+x for x,y in pool_output]
        columns_labels = ['REFERENCE']+[str(_) for _ in sites]
        results_df = pd.DataFrame(formatted_output, index=index_col, columns=columns_labels)

        results_df.to_excel(f'MFE/{spec}_MFE50.xlsx')
        print('____________________________________________________________________________')
        print(f'{spec} IS FINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        now = datetime.now()
        date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print("End Time: ",date_time)
        print('____________________________________________________________________________')
