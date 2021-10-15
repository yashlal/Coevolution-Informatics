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
# import mendotaMFE

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

def entr_ordering(data, sites):
	t_len = len(sites)
	results = []
	for cluster in sites:
		print(str(int(sites.index(cluster))+1) + f'/{t_len}')
		results.append([])
		unsorted_l = []

		X = []
		for j in cluster:
			X.append([d[j] for d in data])
		H_all = joint_entr(*X)

		for i in range(len(cluster)):
			new_cluster = cluster.copy()
			new_cluster.remove(cluster[i])

			Y = []
			for k in new_cluster:
				Y.append([d2[k] for d2 in data])
			H_y = joint_entr(*Y)

			H_ind = shannon_entr([d3[i] for d3 in data])

			statistic = (H_ind-(H_all-H_y))/(H_ind)
			unsorted_l.append(statistic)

		sorted_l = unsorted_l.copy()
		sorted_l.sort()

		for el in sorted_l:
			results[-1].append(cluster[unsorted_l.index(el)])

	return results

def mfe_ordering(data, sites):
	results = []
	for cluster in sites:
		results.append([])
		for step in range(len(cluster)):
			print(results)
			calc_list = []
			for el in cluster:
				seqs = []
				c=0
				i=0
				while c<20:
					seq = data[i]
					if not any([seq[el]=='B']):
						c += 1
						i += 1
						seqs.append(seq)
						continue
					i += 1

				s=0
				for seq in seqs:
					ref_val = RNA.fold(seq)
					exp_val = RNA.fold(mutate(el, seq))
					s += abs(ref_val-exp_val)

				calc_list.append((el,s/50))

				calc_list.sort(key=lambda x: x[1])

			cluster.remove(calc_list[-1][0])
			results.append(calc_list[-1][0])

	for x in results:
		x.reverse()

	return results

if __name__=='__main__':
	# for species in species_list:
	# 	data = load_data(species=species, epsilon=epsilon)
	# 	for gamma in gammas:
	# 		sites = get_sites(species=species, gamma=gamma)
	# 		print(sites)
	# 		print(len(data), len(data[0]))
	# 		RES = entr_ordering(data, sites)
	# 		with open(f'Results/InfRankings/{species}_E{epsilon}_G{gamma}_Ranking.pickle', 'wb') as handle:
	# 			pickle.dump(RES, handle)

	data = load_data(species=species_list[0], epsilon=epsilon)
	sites = get_sites(species=species_list[0], gamma=gammas[0])
	print(mfe_ordering(data,sites))
