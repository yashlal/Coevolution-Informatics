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
import random

species_list = ['Fusobacteriota', 'Cyanobacteria', 'Bacteroidota']
epsilon = 0.116
gammas = [0.95, 0.90, 0.85]

def load_data(species, epsilon):
    datafile = f'data/SILVA_138.1_{species}.fasta'
    filename = f'{species}Results_E_{epsilon}.pickle'

    data_list = getdata(datafile)
    prcsd_data = preprocess(data_list)

    return prcsd_data

def get_sites(species, gamma):
    with open(f'Results/{species}Results_E_0.116.pickle', 'rb') as handle1:
        b1 = pickle.load(handle1)

    sites = list(filter(lambda x: type(x)==list, b1[gamma]))
    return sites

def ordering(data, sites):
	t_len = len(sites)
	results = []
	for cluster in sites:
		print(str(int(sites.index(cluster))+1) + f'/{t_len}')
		results.append([])
		for i in range(len(cluster)):
			new_cluster = cluster.copy()
			new_cluster.remove(cluster[i])
			X = []
			for j in cluster:
				X.append([d[j] for d in data])
			Y = []
			for k in new_cluster:
				Y.append([d2[k] for d2 in data])
			if type(X[0])!=list:
				X=[X]
			if type(Y[0])!=list:
				Y=[Y]
			H_x = joint_entr(*X)
			H_y = joint_entr(*Y)
			results[-1].append(H_x-H_y)
	return results

for species in species_list:
	data = load_data(species=species, epsilon=epsilon)
	for gamma in gammas:
		sites = get_sites(species=species, gamma=gamma)
		print(sites)
		print(len(data), len(data[0]))
		RES = ordering(data, sites)

		with open(f'Results/InfRankings/{species}_E{epsilon}_G{gamma}_Ranking.pickle', 'wb') as handle:
			pickle.dump(RES, handle)
