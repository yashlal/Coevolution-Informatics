from tqdm import tqdm
from Bio import SeqIO
import time
import numpy as np
from tqdm import tqdm
from modules import *
import csv
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

problems = ['.', '-', 'N', 'Y', 'M', 'S', 'K', 'R', 'W', 'V']
bases_dict = {'A': 0, 'G':1, 'C':2, 'T':3, 'B':4, '.':4, '-':4, 'N':4, 'Y':4, 'M':4, 'S':4, 'K':4, 'R':4, 'W':4, 'V':4}

data_list = []
data = SeqIO.parse("data/SILVA_138.1_Fusobacteriota.fasta","fasta")
for sample in data:
  data_list.append(str(sample.seq))

def preprocess(data):
    new_data = []
    pbar_prcs = tqdm(range(len(data)))
    pbar_prcs.set_description('Preprocessing Data Blanks')
    for ind in pbar_prcs:
        new_s = data[ind]
        for char in problems:
            new_s = new_s.replace(char, 'B')
        new_data.append(new_s)
    return new_data

# returns dictionary where key is index and value is the probability
# only consider sites where H(X) > epsilon
def filter_stable_sites(data, epsilon):
    pbar = tqdm(range(len(data)))
    pbar.set_description('Generating Filtered Data')
    indices = []
    new_data = []
    for i in pbar:
        col = [x[i] for x in data]
        se = shannon_entr(col)
        if se>epsilon:
            indices.append(i)
            new_data.append(col)
    return indices, new_data

def numericalize(cols):
    new_cols = []
    for col in cols:
        new_cols.append([])
        for el in col:
            new_cols[-1].append(bases_dict[el])
    return new_cols

if __name__=='__main__':
    epsilon, gamma = 0.0232, 0.5
    path = 'KMC.csv'

    prcsd_data = preprocess(data_list)
    inds, cols = filter_stable_sites(data=prcsd_data, epsilon=epsilon)

    new_cols = numericalize(cols)
    print('Running PCA')
    transformed_data = PCA(n_components=2).fit_transform(new_cols)
    print('Running KMC')
    labels = KMeans().fit_predict(transformed_data)

    fields = ['Column Index', 'KMC Cluster']
    rows = []
    for i in range(len(labels)):
        rows.append([inds[i], labels[i]])
    with open(path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        writer.writerows(rows)
