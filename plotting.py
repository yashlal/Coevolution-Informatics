import pickle
from Bio import SeqIO
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import algorithm
import modules
import time

with open('Results\FusobacteriotaResults_E_0.116.pickle', 'rb') as handle:
    b = pickle.load(handle)

# converts columns values to numbers to do PCA + clustering
def numericalize(cols):
    new_cols = []
    for col in cols:
        new_cols.append([])
        for el in col:
            new_cols[-1].append(algorithm.bases_dict[el])
    return new_cols

# takes a set of KMC labels and returns them as a dictionary of form site_size:n_sites (see histograms in Plots)
def labels_to_hist_dict(labels):
    dict1 = {}
    dict2 = {}

    for el in labels:
        dict1[el] = dict1.get(el,0) + 1
    for site in dict1.items():
        dict2[site[1]] = dict2.get(site[1], 0) + 1

    return dict2

def algplots(new_cols):

    for x in range(len(list(b.items()))):
        mytuple = list(b.items())[x]
        sites = [_ for _ in mytuple[1] if type(_)==list]
        gamma = round(mytuple[0], 2)
        n_sites = len(sites)

        iter_cols_variable = [new_cols[inds.index(i)] for i in list(modules.flatten(sites))]
        transformed_data = PCA(n_components=2).fit_transform(iter_cols_variable)
        labels = KMeans(n_clusters=n_sites, n_init=50).fit_predict(transformed_data)
        graphing_dict = labels_to_hist_dict(labels)

        hist_dict = {}
        for site in sites:
            hist_dict[len(site)] = hist_dict.get(len(site), 0) + 1

        print(gamma, f'Algorithm:{hist_dict}', f'KMC:{graphing_dict}')

        plt.bar(list(hist_dict.keys()), list(hist_dict.values()), width=0.8, color='maroon')
        plt.bar(list(graphing_dict.keys()), list(graphing_dict.values()), width=0.5, color='navy')
        plt.title(f'Gamma={gamma} with {n_sites} Sites')
        plt.xlabel('Size of Sites')
        plt.ylabel('Number of Sites')

        colors = {'Algorithm':'maroon', 'KMC':'navy'}
        labels = list(colors.keys())
        handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]
        plt.legend(handles, labels)
        plt.savefig(f'Plots/Fusobacteriota/FusobacteriotaGamma{gamma}.png')
        plt.clf()


if __name__=='__main__':
    epsilon = 0.0232
    filepath = 'data/SILVA_138.1_Fusobacteriota.fasta'

    data_list = algorithm.getdata(filepath)
    prcsd_data = algorithm.preprocess(data_list)
    inds, cols = algorithm.filter_stable_sites(data=prcsd_data, epsilon=epsilon)
    new_cols = numericalize(cols)

    algplots(new_cols)
