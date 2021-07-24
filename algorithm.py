from tqdm import tqdm
from Bio import SeqIO
import time
import numpy as np
from tqdm import tqdm
from modules import *
import csv
import matplotlib.pyplot as plt

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

fields = ['MI_max_term']
with open('output.csv', 'w') as f:
    csvwriter = csv.writer(f)
    csvwriter.writerow(fields)

def alg(MI_list, gamma, indices, cols):
    t = 0
    max_values = []
    print('\n')
    print('Beginning Algorithm\n')
    while (len(MI_list)>0) and (MI_list[-1][-1]>=gamma):
        max_values.append(MI_list[-1][-1])
        start = time.time()
        print(f'On time t={t}')
        print(f'Max MI Term is {MI_list[-1][-1]}')

        with open('output.csv', 'a') as f:
            csvwriter = csv.writer(f)
            csvwriter.writerow([MI_list[-1][-1]])

        #join the sites
        ind_bound_site = []
        cols_bound_site = []
        for el in MI_list[-1][:-1]:
            if type(el)!=list:
                i = indices.index(el)
                col = cols[i]
                ind_bound_site.append(el)
                cols_bound_site.append(col)
                indices.remove(el)
                cols.remove(col)
            else:
                i = indices.index(el)
                col = cols[i]
                [ind_bound_site.append(x) for x in el]
                [cols_bound_site.append(c) for c in col]
                indices.remove(el)
                cols.remove(col)
        indices.append(ind_bound_site)
        cols.append(cols_bound_site)

        #delete MI calcs with the old columns
        for j in MI_list.copy():
            f1 = set(flatten(MI_list[-1][:-1]))
            f2 = list(flatten(j[:-1]))
            if not f1.isdisjoint(f2):
                MI_list.remove(j)

        #generate new MI elements
        for k in range(len(cols)-1):
            mut_inf_ = mutual_inf(cols[k], cols[-1])
            MI_list.append([indices[k], indices[-1], mut_inf_])

        MI_list = sorted(MI_list, key=lambda x:x[-1])

        print(f'Iteration took {round((time.time()-start), 3)} seconds')
        print("\n")

        t += 1
    return MI_list, indices, cols, max_values

if __name__=='__main__':
    epsilon, gamma = 0.0232, 0.5
    prcsd_data = preprocess(data_list)
    inds, cols = filter_stable_sites(data=prcsd_data, epsilon=epsilon)
    print(f'{len(inds)} Unstable Sites Found!')
    mi_init = gen_mut_inf_mat(indices=inds, cols=cols)
    mi_final, inds_final, cols_final, mv = alg(MI_list=mi_init, gamma=gamma, indices=inds, cols=cols)
    print(mi_final, inds_final)
    plt.plot(mv)
    plt.savefig('MaxValues.png')
    plt.show()
