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

if __name__=='__main__':
    epsilon = 0.0232
    indices, proper_data = filter_stable_sites(data=data_list, epsilon=epsilon)
    print(len(indices))
