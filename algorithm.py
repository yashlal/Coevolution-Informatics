from tqdm import tqdm
from Bio import SeqIO
import time
import numpy as np
from tqdm import tqdm
from modules import *
import matplotlib.pyplot as plt
import multiprocessing
from itertools import combinations
import pickle

problems = ['.', '-', 'N', 'Y', 'M', 'S', 'K', 'R', 'W', 'V', 'D', 'H']
bases_dict = {'A': 0, 'G':1, 'C':2, 'T':3, 'B':4, '.':4, '-':4, 'N':4, 'Y':4, 'M':4, 'S':4, 'K':4, 'R':4, 'W':4, 'V':4, 'D':4, 'H':4}

def getdata(path):
    data_list = []
    data = SeqIO.parse(path,"fasta")
    for sample in data:
        data_list.append(str(sample.seq))
    return data_list

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
    lg = int(len(data[0]))
    pbar = tqdm(range(lg))
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

#creates a list of 3-element lists where the first two elements are the corresponding column indices and the third is the I/min value
#for the future we can merge these lists and keep the last element (index= -1) to be the sorting value
def gen_mut_inf_mat(indices, cols):
    init_u = []
    pbar2 = tqdm(range(len(indices)))
    pbar2.set_description('Doing Mutual Info Calcs')
    for i in pbar2:
        for j in range(i):
            mut_inf_ij = mutual_inf(cols[i], cols[j])
            init_u.append([indices[i], indices[j], mut_inf_ij])
    return sorted(init_u, key=lambda x:x[-1])

def mutinf_setup(indices, cols):
    lg = len(cols)
    custom_indices = combinations(list(range(lg)),2)
    myargs = [(cols[tup1], cols[tup2], indices[tup1], indices[tup2]) for tup1, tup2 in custom_indices]

    with multiprocessing.Pool() as pool:
        MI_calcs = pool.starmap(mutual_inf_MP, myargs)

    return sorted(MI_calcs, key=lambda x:x[-1])

def list_comp(s,l):
    f = list(flatten(l))[:-1]
    if s.isdisjoint(f):
        return l

def alg(MI_list_, indices_, cols_, gammas_):
    gammas_ = sorted(gammas_, reverse=True)
    gamma_ind = 0
    results_dict = {}

    MI_list = MI_list_.copy()
    indices = indices_.copy()
    cols = cols_.copy()

    t = 0
    max_values = []

    print('\n')
    print('Beginning Algorithm\n')

    while (len(MI_list)>0) and (gamma_ind<len(gammas_)):
        start = time.time()
        max_values.append(MI_list[-1][-1])

        print(f'On time t={t}')
        print(f'Max MI Term is {MI_list[-1][-1]}')

        if max_values[-1] < gammas_[gamma_ind]:
            results_dict[gammas_[gamma_ind]] = indices.copy()
            print('\n')
            print('-------------------------------------------------------')
            print('-------------------------------------------------------')
            print('-------------------------------------------------------')
            print(f'REACHED GAMMA {gammas_[gamma_ind]} AT TIME STEP {t}')
            print('-------------------------------------------------------')
            print('-------------------------------------------------------')
            print('-------------------------------------------------------')
            print('\n')
            gamma_ind += 1

        #join the sites
        ind_bound_site = []
        cols_bound_site = []
        for el in MI_list[-1][:-1]:
            i = indices.index(el)
            col = cols[i]

            if type(el)!=list:
                ind_bound_site.append(el)
                cols_bound_site.append(col)
            else:
                [ind_bound_site.append(x) for x in el]
                [cols_bound_site.append(c) for c in col]

            indices.remove(el)
            cols.remove(col)
        indices.append(ind_bound_site)
        cols.append(cols_bound_site)
        #delete MI calcs with the old columns
        f1 = set(flatten(MI_list[-1][:-1]))

        myargs1 = [(f1,j) for j in MI_list.copy()]
        myargs2 = [(cols[-1], cols[k], indices[-1], indices[k]) for k in range(len(cols)-1)]

        with multiprocessing.Pool() as pool:
            MI_list = list(filter(None, pool.starmap(list_comp, myargs1)))
            results = pool.starmap(mutual_inf_MP, myargs2)

        MI_list = sorted(MI_list+results, key=lambda x:x[-1])
        print(f'Iteration took {round((time.time()-start), 3)} seconds')
        t += 1

    return results_dict, max_values

if __name__=='__main__':
    epsilon, gammas = 0.116, [0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
    all_species = ['Cyanobacteria']

    for species in all_species:
        print(f'-------------------------------RUNNING SPECIES {species}!-------------------------------')
        datafile = f'data/SILVA_138.1_{species}.fasta'
        filename = f'{species}Results_E_{epsilon}.pickle'

        data_list = getdata(datafile)
        prcsd_data = preprocess(data_list)
        inds, cols = filter_stable_sites(data=prcsd_data, epsilon=epsilon)
        print(f'{len(inds)} Unstable Sites Found!')
        mi_init = mutinf_setup(indices=inds, cols=cols)
        rd, mv = alg(MI_list_=mi_init, indices_=inds, cols_=cols, gammas_=gammas)

        with open(filename, 'wb') as handle:
            pickle.dump(rd, handle, protocol=pickle.HIGHEST_PROTOCOL)

        plt.plot(mv)
        plt.savefig(f'{species}_MV_E{epsilon}.png')
