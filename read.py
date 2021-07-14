from Bio import SeqIO
import pandas as pd
import numpy as np
from tqdm import tqdm
from modules import *

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

# returns dictionary where key is index and value is the probability
# only consider sites where H(X) > epsilon
def filter_stable_sites(prob_list, epsilon):
    pbar2 = tqdm(range(len(prob_list)))
    pbar2.set_description('Generating Filtered Data')
    indices = []
    new_prob_list = []
    for i in pbar2:
        se = shannon_entr(prob_list[i])
        if se>epsilon:
            new_prob_list.append(prob_list[i])
            indices.append(i)
    return indices, new_prob_list

#creates a list of 3-element lists where the first two elements are the corresponding column indices and the third is the I/min value
#for the future we can merge these lists and keep the last element to be the sorting value
def gen_mut_inf_mat(indices, prob_list):
    ln = len(prob_list)
    pbar3 = tqdm(range(ln))
    pbar3.set_description('Generating MI Matrix')
    mut_inf_list = []
    for i in pbar3:
        for j in range(i):
            minimum = min(shannon_entr(prob_list[i]), shannon_entr(prob_list[j]))
            mut_inf_list.append([indices[i], indices[j], (mutual_inf(prob_list[i], prob_list[j]) / minimum)])
    #sort by last value
    return sorted(mut_inf_list, key=lambda x:x[-1])

if __name__=='__main__':
    epsilon = 0.1
    prob_list = calc_probs(data_list=data_list[1000:1200])
    stable_indices, filtered_prob_list = filter_stable_sites(prob_list, epsilon)
    sorted_inf = gen_mut_inf_mat(stable_indices, filtered_prob_list)
    print(sorted_inf)
