from Bio import SeqIO
import pandas as pd
import numpy as np
from tqdm import tqdm

data_list = []
data = SeqIO.parse("data/SILVA_138.1_Fusobacteriota.fasta","fasta")
for sample in data:
    data_list.append(sample.seq)

#2d list representing probabilities: the ith index gives the probabilities for that position
#each probability list has 5 elements: ordered for the probailities of A G C T BLANK

edgecases = ['.', '-', 'N', 'Y', 'M', 'S', 'K', 'R', 'W', 'V']
dict = {'A': 0, 'G':1, 'C':2, 'T':3}
l = ['A', 'G', 'C', 'T']
prob = []

pbar = tqdm(range(len(data_list[0])))
for c in pbar:
    col = [x[c] for x in data_list]
    prob.append([0,0,0,0,0])
    for i in col:
        if i in edgecases:
            prob[-1][-1] += 1
        else:
            prob[-1][dict[i]] += 1

    prob[-1] = [x/sum(prob[-1]) for x in prob[-1]]
