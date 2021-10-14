from algorithm import getdata, preprocess
from tqdm import tqdm
from Bio import SeqIO
from modules import *
import time
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from itertools import combinations
import pickle

species = 'Cyanobacteria'
epsilon=0.116
def load_data():
    datafile = f'data/SILVA_138.1_{species}.fasta'
    filename = f'{species}Results_E_{epsilon}.pickle'

    data_list = getdata(datafile)
    prcsd_data = preprocess(data_list)

    return prcsd_data

def get_sites():
    gamma = 0.90
    with open(f'Results/{species}Results_E_0.116.pickle', 'rb') as handle1:
        b1 = pickle.load(handle1)

    sites = list(filter(lambda x: type(x)==list, b1[gamma]))
    return sites

def ordering(data, sites):
    t_len = len(sites)
    results = []
    for cluster in sites:
results.append([])
	for i in range(len(cluster)):
	new_cluster = cluster.copy()
	    new_cluster.remove(cluster[i])
	    full_data = [data[j] for j in cluster]
	    mod_data = [data[k] for k in new_cluster]
	    print(full_data, mod_data)
    return results
#data = load_data()
#sites = get_sites()
data = ['AGTCGTAGGT', 'TTTAGACGTG', 'TGTACAAGGA', 'AAGTGTCCCG', 'AAGTCGCATG'] 
sites = [[1,2],[2,4],[0,2,4]]
RES = ordering(data, sites)
print(RES)
