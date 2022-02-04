import pickle
import matplotlib.pyplot as plt
from Bio import SeqIO
import modules
import itertools
import figures
import numpy as np
import scipy

def hist_sites(specs, epsilon, gamma):
    all_mfe = []
    for spec in specs:
        with open(f'Results/Dump/2{spec}_sites_MFE_E{epsilon}_G{gamma}.pickle', 'rb') as handle:
            b = pickle.load(handle)
            b = list(filter(lambda x: x[1]!='SKIP', b))

            [all_mfe.append(figures.mfe_metrics(el[1])) for el in b]

    return all_mfe

def hist_nonsites(specs, epsilon, gamma):
    pairs_used = []
    all_mfe = []

    for spec in specs:
        sites = []
        with open(f'Results/AlgE{epsilon}/{spec}Results_E_{epsilon}.pickle', 'rb') as handle:
            r = pickle.load(handle)[gamma]

            for el in r:
                if type(el)==list:
                    [sites.append(p) for p in list(itertools.permutations(el,2))]


        for n in range(4):
            with open(f'Results/Dump/2{spec}_{n}.pickle', 'rb') as handle:
                b = pickle.load(handle)

                for el in b:
                    pair1 = (el[0], el[1])
                    pair2 = (el[1], el[0])

                    if (pair1 not in sites) and (pair1 not in pairs_used):
                        vals = el[5:]
                        all_mfe.append(figures.mfe_metrics(vals))

                        pairs_used.append(pair1)
                        pairs_used.append(pair2)
    return all_mfe

def hist_stable_sites(specs, epsilon):
    pairs_used = []
    all_mfe = []

    for spec in specs:
        unstable_sites = []
        with open(f'Results/AlgE{epsilon}/{spec}Results_E_{epsilon}.pickle', 'rb') as handle:
            r = pickle.load(handle)[0.95]
            r = list(modules.flatten(r))
            [unstable_sites.append(e) for e in r]


        for n in range(4):
            with open(f'Results/Dump/2{spec}_{n}.pickle', 'rb') as handle:
                b = pickle.load(handle)

                for el in b:
                    pair1 = (el[0], el[1])
                    pair2 = (el[1], el[0])

                    if (pair1 not in pairs_used) and ((pair1[0] not in unstable_sites) and (pair1[1] not in unstable_sites)):
                        vals = el[5:]
                        all_mfe.append(figures.mfe_metrics(vals))

                        pairs_used.append(pair1)
                        pairs_used.append(pair2)
    return all_mfe

all_data = []

if __name__=='__main__':
    epsilons = [0.116,0.232,1.16]
    gamma = [0.95, 0.90, 0.85, 0.80, 0.75, 0.70]
    for epsilon in epsilons:
        for gamma in gammas:

            specs = ['Fusobacteriota', 'Cyanobacteria', 'Bacteroidota']

            all_mfe_sites = hist_sites(['Fusobacteriota', 'Cyanobacteria', 'Bacteroidota'], epsilon, gamma)
            all_mfe_nonsites = hist_nonsites(['Fusobacteriota', 'Cyanobacteria', 'Bacteroidota'], epsilon, gamma)
            all_mfe_unstable = hist_stable_sites(['Fusobacteriota', 'Cyanobacteria', 'Bacteroidota'], epsilon)
            all_mfe_stable = all_mfe_sites+all_mfe_nonsites
            # print(len(all_mfe_sites), len(all_mfe_nonsites))
            #
            # bins_max = max(max(all_mfe_sites), max(all_mfe_nonsites))
            # bins_min = min(min(all_mfe_sites), min(all_mfe_nonsites))
            # bins_step = (bins_max-bins_min)/20
            # bins=np.arange(bins_min, bins_max, bins_step)
            #
            # plt.hist(all_mfe_sites, bins, alpha=0.35, label='Sites', edgecolor='black', density=True)
            # plt.hist(all_mfe_nonsites, bins, alpha=0.35, label='Non-Sites', edgecolor='black', density=True)
            # plt.hist(all_mfe_unstable,bins, alpha=0.35, label='Stable Sites', density=True)
            # plt.legend(loc='upper right')
            #
            # plt.show()

            hp1 = epsilon
            hp2 = gamma

            a = len(all_mfe_sites)
            b = len(all_mfe_nonsites)
            c = len(all_mfe_unstable)

            ks_v, p1 = scipy.stats.kstest(all_mfe_sites, all_mfe_nonsites)
            ks_v, p2 = scipy.stats.kstest(all_mfe_stable, all_mfe_unstable)

            l = [hp1,hp2,a,b,c,p1,p2]
            all_data.append(l)

with open('p_values.pickle', 'rb') as handle:
    pickle.dump(all_data, handle)
