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
import modules
import algorithm
import math

bases = ['A','G','C','T']
pair_combs = ['AA','AG','AC','AT','GA','GG','GC','GT','CA','CG','CC','CT','TA','TG','TC','TT']

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

def ecoli_mp_func(ecoli_seq, ecoli_pair):
    all_seqs = all_mutations_pair(ecoli_seq, ecoli_pair)
    vals = run_pair(all_seqs)

    return (ecoli_pair, vals)

def run_specs(species=['Fusobacteriota', 'Cyanobacteria', 'Bacteroidota'], cap=3000):
    for spec in species:
        with open(f'data/SILVA_138.1_{spec}.pickle', 'rb') as handle1:
            data_list = pickle.load(handle1)
        with open(f'Results/AlgE0.232/{spec}Results_E_0.232.pickle', 'rb') as handle2:
            b = pickle.load(handle2)[0.95]

        sites = list(modules.flatten(list(filter(lambda x: type(x)==list, b))))

        pairs_real_real = list(itertools.combinations(sites,2))

        pairs_real_real = rd.sample(pairs_real_real, min(cap, len(pairs_real_real)))

        pool_input = [(input_pair, data_list) for input_pair in pairs_real_real]
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

def run_ecoli(path='data/4ybb.fasta', max_iter=6, total=30000):
    ecoli_seq_raw = str(list(SeqIO.parse(path, 'fasta'))[0].seq)
    ecoli_psxns = []
    ecoli_seq = []

    for ind in range(len(ecoli_seq_raw)):
        if ecoli_seq_raw[ind] in bases:
            ecoli_psxns.append(ind)
            ecoli_seq.append(ecoli_seq_raw[ind])

    all_ecoli_pairs = list(itertools.combinations(ecoli_psxns, 2))
    all_ecoli_pairs = [(ecoli_psxns.index(a), ecoli_psxns.index(b)) for a,b in all_ecoli_pairs]
    rd.shuffle(all_ecoli_pairs)
    rd.shuffle(all_ecoli_pairs)
    rd.shuffle(all_ecoli_pairs)
    chunskize = total / max_iter
    ecoli_pairs = rd.sample(all_ecoli_pairs, total)

    iter_num=0
    while iter_num<max_iter:
        itervar1 = int(iter_num*chunskize)
        itervar2 = int((iter_num+1)*chunskize)
        iter_pairs = ecoli_pairs[itervar1:itervar2]

        pool_input = [(ecoli_seq, ecoli_pair) for ecoli_pair in iter_pairs]
        with multiprocessing.Pool() as pool:
            pool_output = pool.starmap(ecoli_mp_func, pool_input)

        with open(f'Results/Dump/ecoli_{iter_num}.pickle', 'wb') as ecoli_handle:
            pickle.dump(pool_output, ecoli_handle)

        print('____________________________________________________________________________')
        print(f'E.Coli {iter_num+1}/{max_iter} IS FINISHED')
        now = datetime.now()
        date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print("End Time: ", date_time)
        print('____________________________________________________________________________')
        iter_num += 1

def spec_pairwise_MI_MFE_mp_func(col1, col2, seq, pxns):
    #1. get Ha, Hb, Hab
    entr_result = modules.mutual_inf_MP_new(col1, col2)

    #2. Do mfe run on the seq
    clean_seq, converted_pair = setup(seq, pxns)
    all_seqs = all_mutations_pair(clean_seq, converted_pair)
    vals = run_pair(all_seqs)

    return list(modules.flatten([pxns, entr_result, vals]))

def spec_pairwise_MI_MFE(ecoli_path='data/4ybb.fasta', spec='Fusobacteriota', n=10000, save=5000):
    ecoli_seq_raw = str(list(SeqIO.parse(ecoli_path, 'fasta'))[0].seq)
    ecoli_psxns = []

    for ind in range(len(ecoli_seq_raw)):
        if ecoli_seq_raw[ind] in bases:
            ecoli_psxns.append(ind)

    all_pairs = list(itertools.combinations(ecoli_psxns, 2))

    rd.shuffle(all_pairs)
    rd.shuffle(all_pairs)
    rd.shuffle(all_pairs)

    full_sample = rd.sample(all_pairs, n)
    iter_number = int(n/save)

    spec_path = f'data/SILVA_138.1_{spec}.fasta'
    MSA = [str(sample.seq) for sample in list(SeqIO.parse(spec_path, 'fasta'))]
    MSA = algorithm.preprocess(MSA)

    for loop in range(iter_number):
        lower_bound = int(loop*save)
        upper_bound = int((loop+1)*save)

        loop_pairs = full_sample[lower_bound:upper_bound]

        pool_input = []
        for pair in loop_pairs:
            c=0
            flag=False
            while (not flag) and (c<len(MSA)):
                if (MSA[c][pair[0]] in bases) and (MSA[c][pair[1]] in bases):
                    col1 = ''.join([seq[pair[0]] for seq in MSA])
                    col2 = ''.join([seq[pair[1]] for seq in MSA])
                    pool_input.append((col1, col2, MSA[c], pair))
                    flag=True
                else:
                    c += 1

        with multiprocessing.Pool() as pool:
            pool_output = pool.starmap(spec_pairwise_MI_MFE_mp_func, pool_input)

        with open(f'Results/Dump/2{spec}_{loop}.pickle', 'wb') as handle:
            pickle.dump(pool_output, handle)

        now = datetime.now()
        date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print(f'Finished {spec}: {loop+1}/{iter_number}')
        print(date_time)
        print('______________________________________________')

if __name__=='__main__':
    #printing start time
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")

    run_specs()
