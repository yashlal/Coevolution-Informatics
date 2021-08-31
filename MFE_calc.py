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


bases = ['A','G','C','T']
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
    p = Popen('C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe', stdin=PIPE, stdout=PIPE)
    reference = list(p.communicate(ref_input.encode())[0].decode())
    ref_val = process_val(reference)

    MFE_list = []
    for x in pxns:
        if x=='BLANK':
            MFE_list.append(np.nan)
        else:
            inp = ''.join(mutate(x,s))
            p = Popen('C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe', stdin=PIPE, stdout=PIPE)
            ans = list(p.communicate(inp.encode())[0].decode())
            value_for_dict = process_val(ans)
            MFE_list.append(value_for_dict-ref_val)

    return MFE_list, ref_val

def run_sequence(seq):
    site_nums = []
    pure_seq_l = []

    for i in range(len(seq)):
        if seq[i] not in blanks:
            site_nums.append(i)
            pure_seq_l.append(seq[i])

    sites_2 = []
    nonsites_2 = []
    for site in sites:
        try:
            sites_2.append(site_nums.index(site))
        except:
            sites_2.append('BLANK')
    for nonsite in nonsites:
        try:
            nonsites_2.append(site_nums.index(nonsite))
        except:
            nonsites_2.append('BLANK')

    MFE_pos, rv = MFE_func_all(pure_seq_l, sites_2)
    MFE_neg, rv = MFE_func_all(pure_seq_l, nonsites_2)

    return MFE_pos, MFE_neg, rv

blanks = ['-','.']
species = ['Bacteroidota']
b = []
for spec in species:
    with open(f'Results\{spec}Results_E_0.116.pickle', 'rb') as handle:
        b.append(pickle.load(handle)[0.8])

sites = []
nonsites = []

for i in b[0]:
    if type(i)==list:
        sites.append(i)
    else:
        nonsites.append(i)

sites=list(modules.flatten(sites))

sites=sites
nonsites=nonsites

with open('data\SILVA_138.1_Bacteroidota.fasta') as handle:
    data = list(SeqIO.parse(handle, "fasta"))[0:50]

bad_seqs = [0, 1, 5, 17, 18, 19, 24, 26, 29, 31, 36, 39, 40, 48, 56, 61, 72, 84, 92, 93, 98]

if __name__=='__main__':

    results_l = []
    index_col = []

    pool_input = []

    for j in range(35):
        if j not in bad_seqs:
            index_col.append(j)
            print(f'Sequence {j}/50')
            pool_input.append(str(data[j].seq))

    pool_output = process_map(run_sequence, pool_input)
    formatted_output = [[z]+x+y for x,y,z in pool_output]
    columns_labels = ['REFERENCE']+['S'+str(_) for _ in sites]+['NS'+str(n_) for n_ in nonsites]
    results_df = pd.DataFrame(formatted_output, index=index_col, columns=columns_labels)

    results_df.to_excel('BD_MFE25.xlsx')
