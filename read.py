from Bio import SeqIO
import pandas as pd
import numpy as np
from tqdm import tqdm
from numba import njit

data_list = []
data = SeqIO.parse("data/SILVA_138.1_Fusobacteriota.fasta","fasta")
for sample in data:
    data_list.append(sample.seq)



edgecases = ['.', '-', 'N', 'Y', 'M', 'S', 'K', 'R', 'W', 'V']
dict = {'A': 0, 'G':1, 'C':2, 'T':3}
l = ['A', 'G', 'C', 'T']

# returns 2d list representing probabilities: the ith index gives the probabilities for that position
# each probability list has 5 elements: ordered for the probailities of A G C T BLANK
def calc_probs(data_list):
    prob = []
    pbar = tqdm(range(len(data_list[0])))
    pbar.set_description('Calculating Probabilities')
    for c in pbar:
        col = [x[c] for x in data_list]
        prob.append([0,0,0,0,0])
        for i in col:
            if i in edgecases:
                prob[-1][-1] += 1
            else:
                prob[-1][dict[i]] += 1

        prob[-1] = [x/sum(prob[-1]) for x in prob[-1]]

    return prob

# returns list of all entropy values
# entropy is shannon entropy with logbase e
def get_entr(prob_list):
    pbar2 = tqdm(range(len(prob_list)))
    pbar2.set_description('Calculating Entropies')
    entr_list = []
    for i in pbar2:
        s = 0
        for el in prob_list[i]:
            if el!=0:
                s -= (el*np.log(el))
        entr_list.append(s)
    return entr_list

prob_list = calc_probs(data_list=data_list[12000:20000])
entr_list = get_entr(prob_list=prob_list)
