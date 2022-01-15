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

                    possibilities.append(abs(E_ab))
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

def ecoli_get_dist_and_mfe_metric(items, coords, ref_l):
    all_dist = []
    all_mfe = []
    for item in items:
        pair, values = item

        try:
            dist = np.linalg.norm(coords[pair[0]]-coords[pair[1]])
            all_dist.append(dist)

            mfe = mfe_metrics(values)
            all_mfe.append(mfe)
        except ValueError:
            continue

    return all_dist, all_mfe

def spec_figures(spec='Cyanobacteria', n=10000):

    with open(f'Results/Dump/{spec}_PW_MFE.pickle', 'rb') as handle:
        b = pickle.load(handle)

    coords, ref_l = get_reference_dist()
    count = 0

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

def ecoli_figures(type='hist'):
    all_ecoli_data = []
    for n in range(6):
        with open(f'Results/Dump/ecoli_{n}.pickle', 'rb') as handle:
            b = pickle.load(handle)
            [all_ecoli_data.append(_) for _ in b]

    coords, ref_l = get_reference_dist()
    ecoli_dist, ecoli_mfe = ecoli_get_dist_and_mfe_metric(all_ecoli_data, coords, ref_l)

    var1 = ecoli_dist
    var2 = [el[0] for el in ecoli_mfe]
    var3 = [el[1] for el in ecoli_mfe]

    vars = [var1,var2,var3]
    paths = ['Results/NewFigs/EColi/DistHist', 'Results/NewFigs/EColi/MFEAddHist', 'Results/NewFigs/EColi/MFECompHist']

    if type=='hist':
        for n in range(3):
            var = vars[n]

            if n == 0:
                d = 0.1
            else:
                d = 0.05
                var = list(filter(lambda x: x>0.001, var))
            bins_input = np.arange(math.floor(min(var)), math.ceil(max(var)), d)

            a, b, c = plt.hist(var, bins=bins_input)

            plt.savefig(f'{paths[n]}/hist_nozero.png')
            plt.savefig(f'{paths[n]}/hist_nozero.eps', format='eps')
            plt.clf()

    elif type=='loglin':
        for n in range(3):
            var = vars[n]

            if n == 0:
                d = 0.1
            else:
                d = 0.05
                var = list(filter(lambda x: x>0.001, var))
            bins_input = np.arange(math.floor(min(var)), math.ceil(max(var)), d)

            a, b, c = plt.hist(var, bins=bins_input)
            plt.clf()

            X, Y = [], []
            for i in range(len(a)):
                if a[i]!=0:
                    X.append(b[i])
                    Y.append(np.log(a[i]))

            plt.scatter(X,Y)
            plt.show()
    else:
        pass

if __name__=='__main__':
    ecoli_figures('loglin')
