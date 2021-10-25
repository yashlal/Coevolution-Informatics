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
# import RNA
import mendotaMFE

species_list = ['Cyanobacteria']
epsilon = 0.116
gammas = [0.95, 0.90, 0.85]

def load_data(species, epsilon):
    datafile = f'data/SILVA_138.1_{species}.fasta'
    filename = f'{species}Results_E_{epsilon}.pickle'

    data_list = getdata(datafile)
    prcsd_data = preprocess(data_list)

    return prcsd_data

def get_sites(species, gamma):
    with open(f'Results/AlgE{epsilon}/{species}Results_E_{epsilon}.pickle', 'rb') as handle1:
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

			denom = (H_ind+H_y)/2
			statistic = (H_all-H_y)
			unsorted_l.append(statistic)

		sorted_l = unsorted_l.copy()
		sorted_l.sort()

		for el in sorted_l:
			results[-1].append(cluster[unsorted_l.index(el)])

	return results

def MP_MFE_func(seqs, cluster):
	results = []

	for step in range(len(cluster)):
		calc_list = []
		for postn in cluster:
			s=0
			for seq in seqs:

				new_seq1 = list(seq).copy()
				for el in results:
					new_seq1 = mendotaMFE.mutate(el, new_seq1)

				new_seq2 = new_seq1.copy()
				new_seq1 = list(filter(lambda z: z!='B', new_seq1))

				new_seq2 = mendotaMFE.mutate(postn, new_seq2)
				new_seq2 = list(filter(lambda z: z!='B', new_seq2))

				ref_val = RNA.fold(''.join(new_seq1))[1]
				exp_val = RNA.fold(''.join(new_seq2))[1]

				s += abs(ref_val-exp_val)

			calc_list.append((postn, s/len(seqs)))

		calc_list.sort(key=lambda x: x[1])

		cluster.remove(calc_list[0][0])
		results.append(calc_list[0][0])

	results.reverse()

	return results

if __name__=='__main__':
	for species in species_list:
		data = load_data(species=species, epsilon=epsilon)
		for gamma in gammas:
			sites = get_sites(species=species, gamma=gamma)
			# print(sites)
			# print(len(data), len(data[0]))
			RES = entr_ordering(data, sites)
			with open(f'Results/InfRankings/R4/{species}_E{epsilon}_G{gamma}_Ranking.pickle', 'wb') as handle:
				pickle.dump(RES, handle)

	# for species in species_list:
	# 	data = load_data(species=species, epsilon=epsilon)
	# 	for gamma in gammas:
	# 		sites = get_sites(species=species, gamma=gamma)
	# 		seqs = []
	# 		for cluster in sites:
	# 			seqs.append([])
	# 			for seq in data:
	# 				if all([seq[el]!='B' for el in cluster]):
	# 					seqs[-1].append(seq)
	#
	# 		seqs = [x[0:20] for x in seqs]
	#
	# 		pool_input = [(seqs[i], sites[i]) for i in range(len(sites))]
	# 		with multiprocessing.Pool() as pool:
	# 			all_results = pool.starmap(MP_MFE_func, pool_input)
	#
	# 		with open(f'{species}_E{epsilon}_G{gamma}_Ranking.pickle', 'wb') as handle:
	# 			pickle.dump(all_results, handle)
