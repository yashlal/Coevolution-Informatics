from Bio import SeqIO
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import histograms
import algorithm
import time

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

if __name__=='__main__':
    itervar = histograms.algplots()
    epsilon = 0.0232
    filepath = 'data/SILVA_138.1_Fusobacteriota.fasta'

    data_list = algorithm.getdata(filepath)
    prcsd_data = algorithm.preprocess(data_list)
    inds, cols = algorithm.filter_stable_sites(data=prcsd_data, epsilon=epsilon)
    new_cols = numericalize(cols)

    for tup in zip(*itervar):
        # Too tired to make up variable names
        # a is the number of cluster, b is the list of indices for the columns to be clustered, c is the gamma value
        a,b,c = tup
        iter_cols_variable = [new_cols[inds.index(i)] for i in b]

        print('Running PCA')
        transformed_data = PCA(n_components=2).fit_transform(iter_cols_variable)

        labels = KMeans(n_clusters=a).fit_predict(transformed_data)
        graphing_dict = labels_to_hist_dict(labels)
        print(c, graphing_dict)
        plt.bar(list(graphing_dict.keys()), list(graphing_dict.values()), width=0.8, color='navy')
        plt.title(f'KMC on Gamma={c} with {a} Sites')
        plt.xlabel('Size of Sites')
        plt.ylabel('Number of Sites')
        plt.savefig(f'Plots/KMCGamma{c}.png')
        plt.clf()
