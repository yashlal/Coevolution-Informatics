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

def setup(raw_seq, all_pairs):
    editable_seq = list(raw_seq)

    new_coords = []
    clean_seq = []
    for i in range(len(editable_seq)):
        bp = editable_seq[i]
        if bp in bases:
            new_coords.append(i)
            clean_seq.append(bp)

    converted_pairs = []
    for pair in all_pairs:
        if (pair[0] in new_coords) and (pair[1] in new_coords):
            converted_pairs.append((new_coords.index(pair[0])), new_coords.index(pair[1]))

    return clean_seq, converted_pairs

def all_mutations_pair(pxns, full_seq):
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
    return vals

def mp_func(pair, seq):
    all_seqs = all_mutations_pair(pair, seq)
    vals = run_pair(all_seqs)

    return (pair, vals)


if __name__=='__main__':
    #printing start time
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print("STARTING:", date_time)

    b = []
    for spec in species:

        data_list = []
        data = SeqIO.parse(f'data/SILVA_138.1_{spec}.fasta',"fasta")
        for sample in data:
            data_list.append(str(sample.seq))

        flag=False
        ind=0
        while flag==False:
            ind += 1
            if all([p in safe for p in data_list[ind]]):
                flag=True

        chosen_seq = seq[ind]

        with open(f'{spec}Results_E_0.232.pickle', 'rb') as handle:
            b.append(pickle.load(handle)[0.95])

        sites = list(modules.flatten(list(filter(lambda x: type(x)==list, b))))
        b = list(modules.flatten(b))
        nonsites = []
        for x in b:
            if x not in sites:
                nonsites.append(x)

        sites = list(filter(lambda x: x in bases, sites))
        nonsites = list(filter(lambda x: x in bases, nonsites))

        l = 10

        pairs_real_real = rd.sample(list(itertools.combinations(sites,2)), l)
        pairs_real_fake = rd.sample(list(itertools.product(sites, nonsites)), l)
        pairs_fake_fake = rd.sample(list(itertools.combinations(nonsites,2)), l)

        all_pairs = pairs_real_real + pairs_real_fake + pairs_fake_fake

        clean_seq, converted_pairs = setup(chosen_seq, all_pairs)


        pool_input = [(input_pair, clean_seq) for input_pair in converted_pairsirs]
        with multiprocessing.Pool() as pool:
            pool_output = pool.starmap(mp_func, pool_input)

        formatted_output = [(tup[0], tup[1], *vals) for tup in pool_output]

        with open(f'Results/Dump/{spec}_PW_MFE.pickle', 'wb') as handle2:
            pickle.dump(formatted_output, handle2)

        print('____________________________________________________________________________')
        print(f'{spec} IS FINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        now = datetime.now()
        date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print("End Time: ", date_time)
        print('____________________________________________________________________________')
