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

    print(len(prcsd_data), len(prcsd_data[0]))
    return prcsd_data

def get_sites():
    gamma = 0.90
    with open(f'Results/{species}Results_E_0.116.pickle', 'rb') as handle1:
        b1 = pickle.load(handle1)

    sites = list(filter(lambda x: type(x)==list, b1[gamma]))
    return sites

data = load_data()
sites = get_sites()
print(sites)
