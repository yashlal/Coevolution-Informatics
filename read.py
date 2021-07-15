from Bio import SeqIO
import pandas as pd
import time
import numpy as np
from tqdm import tqdm
from modules import *

data_list = []
data = SeqIO.parse("data/SILVA_138.1_Fusobacteriota.fasta","fasta")
for sample in data:
    data_list.append(str(sample.seq))

# returns dictionary where key is index and value is the probability
# only consider sites where H(X) > epsilon
def filter_stable_sites(data, epsilon):
    pbar = tqdm(range(len(data)))
    pbar.set_description('Generating Filtered Data')
    indices = []
    new_data = []
    t1 = 0
    t2 = 0
    for i in pbar:
        start = time.time()
        col = [x[i] for x in data]
        se = shannon_entr(col)
        if se>epsilon:
            indices.append(i)
            new_data.append(col)
    return indices, new_data

#creates a list of 3-element lists where the first two elements are the corresponding column indices and the third is the I/min value
#for the future we can merge these lists and keep the last element (index= -1) to be the sorting value
def gen_mut_inf_mat(indices, cols):
    init_u = []
    pbar2 = tqdm(range(len(indices)))
    pbar2.set_description('Doing Mutual Info Calcs')
    for i in pbar2:
        for j in range(i):
            mut_inf_ij = mutual_inf([cols[i]], [cols[j]])
            init_u.append([indices[i], indices[j], mut_inf_ij])
    return sorted(init_u, key=lambda x:x[-1])

if __name__=='__main__':
    epsilon = 0.0232
    indices, proper_data = filter_stable_sites(data=data_list[10000:12000], epsilon=epsilon)
    print(len(indices))
    init_sorted = gen_mut_inf_mat(indices[-10:], proper_data[-10:])
    print(init_sorted)
