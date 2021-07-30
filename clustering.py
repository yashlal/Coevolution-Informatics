from Bio import SeqIO
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import histograms
import algorithm
import time

def numericalize(cols):
    new_cols = []
    for col in cols:
        new_cols.append([])
        for el in col:
            new_cols[-1].append(algorithm.bases_dict[el])
    return new_cols

if __name__=='__main__':
    itervar = histograms.algplots()
    epsilon = 0.0232
    filepath = 'data/SILVA_138.1_Fusobacteriota.fasta'

    data_list = algorithm.getdata(filepath)
    prcsd_data = algorithm.preprocess(data_list)
    inds, cols = algorithm.filter_stable_sites(data=prcsd_data, epsilon=epsilon)
    new_cols = numericalize(cols)

    # print('Running PCA')
    # transformed_data = PCA(n_components=2).fit_transform(new_cols)
    #
    # for n_clus in lns:
    #     print(f'Running KMC with {n_clus} Clusters')
    #     labels = KMeans(n_clusters=n_clus).fit_predict(transformed_data)

    for tup in zip(*itervar):
        iter_cols_variable = [new_cols[inds.index(i)] for i in tup[1]]
        print(iter_cols_variable)
