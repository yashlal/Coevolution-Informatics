import pickle
from algorithm import getdata, preprocess
from scipy.stats import wasserstein_distance
from tqdm import tqdm
from collections import Counter

problems = ['.', '-', 'N', 'Y', 'M', 'S', 'K', 'R', 'W', 'V', 'D', 'H']
bases = ['A','G','C','T']
bases_dict = {'A': 0, 'G':1, 'C':2, 'T':3, 'B':4, '.':4, '-':4, 'N':4, 'Y':4, 'M':4, 'S':4, 'K':4, 'R':4, 'W':4, 'V':4, 'D':4, 'H':4}

species = ['Cyanobacteria', "Bacteroidota"]
results = []
gamma = 0.8
datasets = []

for spec in species:
    with open(f'Results\{spec}Results_E_0.116.pickle', 'rb') as handle:
        b = pickle.load(handle)
        results.append(b[gamma])
    datasets.append(preprocess(getdata(f'data\SILVA_138.1_{spec}.fasta')))

for x in results.copy():
    for y in x.copy():
        try:
            len(y)
        except TypeError:
            results[results.index(x)].remove(y)

def get_distro(d, col_inds):
    cols = []
    for ind in col_inds:
        cols.append([])
        for s in d:
            cols[-1].append(s[ind])
    rows = [row for row in zip(*cols)]
    probs_ar = dict(Counter(rows))
    mysum = sum(list(probs_ar.values()))
    probs_ar = {k:v/mysum for k,v in probs_ar.items()}

    return dict(probs_ar)

distros = []
pbar=tqdm(range(len(results[0])))
for site_ind1 in pbar:
    site = results[0][site_ind1]
    distros.append([get_distro(datasets[0],site), get_distro(datasets[1],site), site])

pbar2=tqdm(range(len(results[1])))
for site_ind2 in pbar2:
    site = results[1][site_ind2]
    distros.append([get_distro(datasets[0],site), get_distro(datasets[1],site), site])

with open('Results/CBFB.pickle', 'wb') as handle:
    pickle.dump(distros, handle, protocol=pickle.HIGHEST_PROTOCOL)
