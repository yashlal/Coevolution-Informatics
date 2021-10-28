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

def new_MFE_func(data, pair):
    c = 0
    j = 0
    sum = 0

    while c<20 and j<len(data):
        s = data[j]

        if (s[pair[0]] in bases) and (s[pair[1]] in bases):
            c += 1

            site_nums = []
            pure_seq_l = []
            for i in range(len(s)):
                if s[i] in bases:
                    site_nums.append(i)
                    pure_seq_l.append(s[i])

            pair_2 = (site_nums.index(pair[0]), site_nums.index(pair[1]))

            ref_input = ''.join(pure_seq_l)
            ref_val = RNA.fold(ref_input)[1]

            inp1 = ''.join(mutate(pair[0], pure_seq_l))
            inp2 = ''.join(mutate(pair[1], pure_seq_l))
            val1 = RNA.fold(inp1)[1]
            val2 = RNA.fold(inp2)[1]
            sum += ((abs(val1-ref_val) + abs(val2-ref_val)) / 2)
        j += 1

    return sum/(c+1)

if __name__=='__main__':
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print("STARTING:", date_time)

    b = []
    for spec in species:
        with open(f'{spec}Results_E_0.232.pickle', 'rb') as handle:
            b.append(pickle.load(handle)[0.95])

        sites = list(filter(lambda x: type(x)==list, b))
        nonsites = []
        for x in sites:
            if x not in sites:
                nonsites.append(x)

        l = math.floor(len(nonsites)/12)
        nonsites = rd.sample(nonsites, l)
        all_sites=sites+nonsites
        pairs = list(itertools.combinations(all_sites,2))

        data_list = []
        data = SeqIO.parse(f'data/SILVA_138.1_{spec}.fasta',"fasta")
        for sample in data:
            data_list.append(str(sample.seq))

        pool_input = [(data_list, pair), for pair in pairs]
        with multiprocessing.Pool() as pool:
            pool_output = pool.starmap(new_MFE_func, pool_input)

        formatted_output = [(pairs[i], pool_output[i]) for i in range(len(pairs))]

        with open(f'{spec}_PW_MFE.pickle', 'wb') as handle2:
            pickle.dump(formatted_output, handle2)

        print('____________________________________________________________________________')
        print(f'{spec} IS FINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        now = datetime.now()
        date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
        print("End Time: ", date_time)
        print('____________________________________________________________________________')
