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
import math

bases = ['A','G','C','T']
safe = ['A','G','C','T', '-', '.']
pair_combs = ['AA','AG','AC','AT','GA','GG','GC','GT','CA','CG','CC','CT','TA','TG','TC','TT']
blanks = ['-','.']
species = ['Cyanobacteria']

def setup(raw_seq, pair):
    editable_seq = list(raw_seq)

    new_coords = []
    clean_seq = []
    for i in range(len(editable_seq)):
        bp = editable_seq[i]
        if bp in bases:
            new_coords.append(i)
            clean_seq.append(bp)

    converted_pair = (new_coords.index(pair[0])), new_coords.index(pair[1])

    return clean_seq, converted_pair

def all_mutations_pair(full_seq, pxns):
    i,j = pxns
    all_seqs = []

    for comb in pair_combs:
        new_seq = full_seq.copy()
        new_seq[i] = comb[0]
        new_seq[j] = comb[1]

        all_seqs.append(new_seq)
    all_seqs.append(full_seq)
    return all_seqs

def run_pair(all_seqs_list):
    vals = []
    for seq in all_seqs_list:
        val = RNA.fold(''.join(seq))[1]
        vals.append(val)
    return vals

def mp_func(pair, data_list):
    flag = False
    c = 0
    while (not flag) and (c<2000):
        if (data_list[c][pair[0]] in bases) and (data_list[c][pair[1]] in bases):
            flag = True
        else:
            c += 1
    if c==2000:
        return (pair, 'SKIP')
    else:
        final_seq = data_list[c]
        clean_seq, converted_pair = setup(final_seq, pair)
        all_seqs = all_mutations_pair(clean_seq, converted_pair)
        vals = run_pair(all_seqs)

        return (pair, vals)


if __name__=='__main__':
    #printing start time
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print("STARTING:", date_time)

    for spec in species:

        data_list = []
        data = SeqIO.parse(f'data/SILVA_138.1_{spec}.fasta',"fasta")
        for sample in data:
            data_list.append(str(sample.seq))
        data_list = data_list[0:2000]
        with open(f'Results/AlgE0.232/{spec}Results_E_0.232.pickle', 'rb') as handle:
            b = pickle.load(handle)[0.95]

        sites = list(modules.flatten(list(filter(lambda x: type(x)==list, b))))
        b = list(modules.flatten(b))
        nonsites = []
        for x in b:
            if x not in sites:
                nonsites.append(x)

        pairs_real_real = list(itertools.combinations(sites,2))
        pairs_real_fake = list(itertools.product(sites, nonsites))
        pairs_fake_fake = list(itertools.combinations(nonsites,2))

        cap = 10000

        pairs_real_real = rd.sample(pairs_real_real, min(cap, len(pairs_real_real)))
        pairs_real_fake = rd.sample(pairs_real_fake, min(cap, len(pairs_real_fake)))
        pairs_fake_fake = rd.sample(pairs_fake_fake, min(cap, len(pairs_fake_fake)))

        all_pairs = pairs_real_real + pairs_real_fake + pairs_fake_fake

        pool_input = [(input_pair, data_list) for input_pair in all_pairs]
        with multiprocessing.Pool() as pool:
            pool_output = pool.starmap(mp_func, pool_input)

        with open(f'Results/Dump/{spec}_PW_MFE.pickle', 'wb') as handle2:
            pickle.dump(pool_output, handle2)

        print('____________________________________________________________________________')
        print(f'{spec} IS FINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        now = datetime.now()
        date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print("End Time: ", date_time)
        print('____________________________________________________________________________')
