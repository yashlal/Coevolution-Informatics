import pickle
import modules
import random as rd
import itertools
from Bio import SeqIO
import time
from distanceROC import get_reference_dist
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import scipy.stats
import math

pair_combs = ['AA','AG','AC','AT','GA','GG','GC','GT','CA','CG','CC','CT','TA','TG','TC','TT']
pairs = ['A','G','C','T']

def mfe_metrics(values_list):
    base = values_list[-1]
    base_pairs = pair_combs[values_list[:-1].index(base)]
    k = pairs.index(base_pairs[0])
    l = pairs.index(base_pairs[1])
    sum = 0

    possibilities = []

    for pair1 in pairs:
        if pair1!=base_pairs[0]:
            for pair2 in pairs:
                if pair2!=base_pairs[1]:
                    i = pairs.index(pair1)
                    j = pairs.index(pair2)

                    E_ab = values_list[4*i+j]-base

                    E_a = values_list[4*i+l]-base
                    E_b = values_list[4*k+j]-base

                    metric_val = E_ab-(E_a+E_b)
                    sum += metric_val

                    choices = [E_ab-E_a, E_ab-E_b]

                    possibilities.append(min(choices))
    return abs(sum/9), min(possibilities)

def get_dist_and_mfe_metric(items, coords, ref_l):
    all_dist = []
    all_mfe = []
    for item in items:
        pair, values = item

        try:
            dist = np.linalg.norm(coords[ref_l.index(pair[0])]-coords[ref_l.index(pair[1])])
            all_dist.append(dist)

            mfe = mfe_metrics(values)
            all_mfe.append(mfe)
        except ValueError:
            continue

    return all_dist, all_mfe

if __name__=='__main__':
    spec = 'Cyanobacteria'
    with open(f'Results/Dump/{spec}_PW_MFE.pickle', 'rb') as handle:
        b = pickle.load(handle)

    coords, ref_l = get_reference_dist()
    count = 0
    n=10000

    ss = b[0:n]
    sn = b[n:2*n]
    nn = b[2*n:]

    datasets = [ss,sn,nn]

    for dataset in datasets:
        dist, mfe = get_dist_and_mfe_metric(dataset, coords, ref_l)
        x = dist
        y1 = [mfe_tuple[0] for mfe_tuple in mfe]
        y2 = [mfe_tuple[1] for mfe_tuple in mfe]

        sample = y2

        # sample = list(filter(lambda x:x>=0.001, sample))
        bins_input = np.arange(math.floor(min(sample)), math.ceil(max(sample)), 0.05)

        plt.hist(sample, bins=bins_input, density=True)
        P = scipy.stats.expon.fit(sample)

        rX = bins_input
        rP = scipy.stats.expon.pdf(rX, *P)
        plt.plot(rX,rP)
        plt.show()
        plt.clf()
