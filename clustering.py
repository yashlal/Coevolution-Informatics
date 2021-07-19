from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from Bio import SeqIO
from tqdm import tqdm
import sys

bases_dict = {'A': 0, 'G':1, 'C':2, 'T':3, '.':4, '-':4, 'N':4, 'Y':4, 'M':4, 'S':4, 'K':4, 'R':4, 'W':4, 'V':4}


data_list = []
data = SeqIO.parse("data/SILVA_138.1_Fusobacteriota.fasta","fasta")
for sample in data:
  data_list.append(str(sample.seq))

def preprocess(data):
    new_data = []
    pbar_prcs = tqdm(range(len(data)))
    pbar_prcs.set_description('Preprocessing Data Blanks')
    for ind in pbar_prcs:
        new_s = []
        for char in data[ind]:
            new_s.append(bases_dict[char])
        new_data.append(new_s)
    return new_data

new_data_list = preprocess(data_list)
print('Transposing Rows and Columns')
transposed_data = list(map(list, zip(*new_data_list)))
print('Running PCA')
transformed_data = PCA(n_components=2).fit_transform(transposed_data)
print('Running KMC')
labels = KMeans().fit_predict(transformed_data)
print(labels)
