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

                    possibilities.append(metric_val)

    return sum/9, min(possibilities)

def new_mfe_metric(values_list):
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

                    if E_ab<0.1:
                        possibilities.append('skip')
                    else:
                        E_a = values_list[4*i+l]-base
                        E_b = values_list[4*k+j]-base

                        metric_val = (E_a+E_b)-E_ab
                        sum += metric_val

                        possibilities.append(abs(metric_val))

    return possibilities

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

def spec_figures(spec='Fusobacteriota', type='scattershade'):
    all_spec_data = []
    for n in range(4):
        with open(f'Results/Dump/{spec}_{n}.pickle', 'rb') as handle:
            b = pickle.load(handle)
            [all_spec_data.append(_) for _ in b]

    all_vals = []
    all_mi = []
    for el in all_spec_data:
        mi = el[0]+el[1]-el[2]
        vals = el[3:]
        all_vals.append(vals)
        all_mi.append(mi)

    spec_mfe = [mfe_metrics(el_vals) for el_vals in all_vals]

    var1 = []
    var2 = [el[0] for el in spec_mfe]
    var3 = [el[1] for el in spec_mfe]
    var4 = all_mi

    vars = [var1,var2,var3,var4]
    paths = ['Results/NewFigs/EColi/DistHist', 'Results/NewFigs/EColi/MFEAddHist', 'Results/NewFigs/EColi/MFECompHist']
    names = ['Distance', 'MFE Metric 1', 'MFE Metric 2', 'Mutual Information']

    if type=='scattershade':
        metric_labels = ['Additivity', 'Compensating Mutation']
        data=[]
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

    elif type=='mi_scattershade':
        metric_labels = ['Additivity', 'Compensating Mutation']
        data=[]
        for n in range(2):
            var = vars[n+1]
            X = vars[3]
            Y = vars[n+1]
            iter_list=[]
            # for iterater in range(len(Y)):
            #     if Y[iterater]<0.001:
            #         iter_list.append(iterater)

            width = 0.1
            bin_edges = np.arange(0,1,width)

            for edge_int in range(len(bin_edges)):
                for i in range(len(Y)):
                    if bin_edges[edge_int]<=X[i]<(bin_edges[edge_int]+width):
                        data.append((bin_edges[edge_int], Y[i], metric_labels[n]))

        df = pd.DataFrame(data, columns=['Mutual Information', 'mean(MFE Metric)', 'MFE Metric'])
        sb.lineplot(x='Mutual Information', y='mean(MFE Metric)', ci='sd', hue='MFE Metric', data=df)
        plt.show()
    elif type=='histogram':
        var = vars[3]

        var = list(filter(lambda x: x>0.01, var))
        bins_input = np.arange(math.floor(min(var)), math.ceil(max(var)), 0.05)

        a, b, c = plt.hist(var, bins=bins_input)
        plt.xlabel(f'{names[n]}')
        plt.ylabel('Frequency')
        plt.show()
    elif type=='scatter':
        for n in range(2):
            plt.scatter(vars[3],vars[n+1])
            plt.show()
            plt.clf()

def new_scatter():
    all_ecoli_data = []
    for n in range(6):
        with open(f'Results/Dump/ecoli_{n}.pickle', 'rb') as handle:
            b = pickle.load(handle)
            [all_ecoli_data.append(_) for _ in b]

    coords, ref_l = get_reference_dist()
    ecoli_dist, ecoli_mfe = ecoli_get_dist_and_mfe_metric(all_ecoli_data, coords, ref_l, flag=False)

    var1 = ecoli_dist
    var1 = [(x,x,x,x,x,x,x,x,x) for x in var1]
    var1 = list(modules.flatten(var1))
    var2 = ecoli_mfe
    var2 = list(modules.flatten(var2))

    newvar1 = []
    newvar2 = []
    for q in range(len(var2)):
        if var2[q]!='skip':
            newvar1.append(var1[q])
            newvar2.append(var2[q])

    X=np.array(newvar1)
    Y=np.array(newvar2)
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(X, Y)
    plt.plot(X, slope*X+intercept)
    print(slope, intercept, r_value)
    plt.scatter(X, Y)
    plt.show()

if __name__=='__main__':
    new_scatter()
