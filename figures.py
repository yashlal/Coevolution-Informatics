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
import seaborn as sb
import pandas as pd

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

    return max(possibilities)

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
        except IndexError:
            print(pair)

    return all_dist, all_mfe

def ecoli_get_dist_and_mfe_metric(items, coords, ref_l, flag=True):
    all_dist = []
    all_mfe = []
    for item in items:
        pair, values = item

        try:
            dist = np.linalg.norm(coords[pair[0]]-coords[pair[1]])
            all_dist.append(dist)
            if flag:
                mfe = mfe_metrics(values)
                all_mfe.append(mfe)
            else:
                mfe = new_mfe_metric(values)
                all_mfe.append(mfe)
        except ValueError:
            continue

    return all_dist, all_mfe

def old_spec_figures(spec='Cyanobacteria', n=10000):

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
    names = ['Distance', 'MFE Metric 1', 'MFE Metric 2']

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
            plt.xlabel(f'{names[n]}')
            plt.ylabel('Frequency')
            plt.savefig(f'{paths[n]}/hist_nozero.png')
            plt.clf()


    elif type=='loglin':
        for n in range(3):

            labels = [('Frequency', 'Distance'), ('MFE Metric 1', 'log(frequency)'), ('MFE Metric 2', 'log(frequency)')]

            var = vars[n]

            if n == 0:
                d = 0.1
            else:
                d = 0.3
                var = list(filter(lambda x: x>0.001, var))
            bins_input = np.arange(math.floor(min(var)), math.ceil(max(var)), d)

            a, b, c = plt.hist(var, bins=bins_input)
            plt.clf()

            X, Y = [], []
            for i in range(len(a)):
                if a[i]!=0:
                    X.append(b[i])
                    Y.append(np.log(a[i]))

            X=np.array(X)
            Y=np.array(Y)
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X, Y)
            print(slope,intercept,r_value)
            plt.plot(X, slope*X+intercept)
            plt.scatter(X,Y)
            plt.xlabel(f'{labels[n][0]}')
            plt.ylabel(f'{labels[n][1]}')
            plt.show()

    elif type=='loglog':
        for n in range(3):
            labels = [('Frequency', 'Distance'), ('log(MFE Metric 1)', 'log(frequency)'), ('log(MFE Metric 2)', 'log(frequency)')]

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
                    X.append(np.log(b[i]))
                    Y.append(np.log(a[i]))

            X=np.array(X)
            Y=np.array(Y)
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X, Y)
            print(slope,intercept,r_value)
            plt.plot(X, slope*X+intercept)
            plt.scatter(X,Y)
            plt.xlabel(f'{labels[n][0]}')
            plt.ylabel(f'{labels[n][1]}')
            plt.show()

    elif type=='scattershade':
        data = []
        metric_labels = ['Additivity', 'Compensating Mutation']
        for n in range(2):
            var = vars[n+1]
            X = vars[0]
            Y = vars[n+1]

            width = 10
            bin_edges = np.arange(0,230,width)

            for edge_int in range(len(bin_edges)):
                for i in range(len(Y)):
                    if bin_edges[edge_int]<=X[i]<(bin_edges[edge_int]+width):
                        data.append((bin_edges[edge_int], Y[i], metric_labels[n]))

        df = pd.DataFrame(data, columns=['Distance', 'mean(MFE Metric)', 'MFE Metric'])
        sb.lineplot(x='Distance', y='mean(MFE Metric)', ci='sd', hue='MFE Metric', data=df)
        plt.show()

    elif type=='old_scatter':
        for n in range(2):
            X=vars[0]
            Y=vars[n+1]
            plt.scatter(X,Y)
            plt.show()

    else:
        pass

def spec_figures(specs=['Bacteroidota', 'Fusobacteriota', 'Cyanobacteria']):
    all_spec_data = []
    for spec in specs:
        with open(f'Results/Dump/{spec}_sites_MFE.pickle', 'rb') as handle:
            b = pickle.load(handle)
            [all_spec_data.append(_) for _ in b]

    all_spec_data = list(filter(lambda x: x[1]!='SKIP', all_spec_data))

    all_vals = []
    for el in all_spec_data:
        vals = el[1]
        all_vals.append(vals)

    all_spec_data2 = []
    for spec in specs:
        for n in range(4):
            with open(f'Results/Dump/2{spec}_{n}.pickle', 'rb') as handle:
                b = pickle.load(handle)
                [all_spec_data2.append(_) for _ in b]

    all_vals2 = []
    for el in all_spec_data2:
        mi = el[2]+el[3]-el[4]

        if mi<0.001:
            vals = el[5:]
            all_vals2.append(vals)

    all_mfe1 = [mfe_metrics(val_el) for val_el in all_vals]
    all_mfe2 = [mfe_metrics(val_el2) for val_el2 in all_vals2]

    all_max = max([max(all_mfe1), max(all_mfe2)])
    all_min = min([min(all_mfe1), min(all_mfe2)])
    custom_step = (all_max-all_min)/10

    custom_bins = np.arange(all_min, all_max, custom_step)

    plt.hist(all_mfe1, bins=custom_bins, density=True, label='Sites', fill=False, edgecolor='red')
    plt.hist(all_mfe2, bins=custom_bins, density=True, label='MI Zero', fill=False, edgecolor='blue')
    plt.legend(loc='upper right')
    plt.show()

if __name__=='__main__':
    spec_figures()
