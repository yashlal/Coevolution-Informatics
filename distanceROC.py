from Bio import PDB
from Bio import SeqIO
from Bio.PDB.Polypeptide import PPBuilder
import pickle
from modules import flatten
from itertools import combinations
import numpy as np
import scikitplot as skplt
import matplotlib.pyplot as plt
import random as rd

blanks = ['-','.']
species_list = ['Fusobacteriota', 'Bacteroidota', 'Cyanobacteria']
epsilons = [0.116, 0.232, 1.16]
gammas = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]

parser = PDB.MMCIFParser()
structure = parser.get_structure('4ybb', 'data/4ybb.cif')
model = structure[0]
chain = model['AA']
reds = [r for r in chain]
coords = [r["C5'"].get_coord() for r in reds[0:1534]]

path1 = 'data/4ybb.fasta'
data_list1 = []
data1 = SeqIO.parse(path1, 'fasta')
for record in data1:
    data_list1.append(str(record.seq))
ref = data_list1[0]
ref_l = []
for char_ind in range(len(ref)):
    if ref[char_ind] not in blanks:
        ref_l.append(char_ind)

for species in species_list:
    for epsilon in epsilons:
        for gamma in gammas:
            path2 = f'Results/AlgE{epsilon}/{species}Results_E_{epsilon}.pickle'
            with open(path2,'rb') as handle:
                b=pickle.load(handle)

            sites = list(flatten(list(filter(lambda x: type(x)==list, b[gamma]))))
            comb = list(combinations(sites,2))
            rd.shuffle(comb)

            all_dist= []
            labels = []
            n=0
            l = len(comb)
            while (n<len(comb)) and (n<10000):
                print(f'{n}/{l}')
                tup = comb[n]
                try:
                    all_dist.append(np.linalg.norm(coords[ref_l.index(tup[0])]-coords[ref_l.index(tup[1])]))
                    c=0
                    for el in b[gamma]:
                        if type(el)==list:
                            if (tup[0] in el) and (tup[1] in el):
                                c = 1
                    labels.append(c)
                except:
                    pass

                n += 1

            new_dist = [i/max(all_dist) for i in all_dist]
            y_probas = [[1-x,x] for x in new_dist]
            print('Graphing')
            ax_inst = skplt.metrics.plot_roc(y_true=labels, y_probas=y_probas, classes_to_plot=[], plot_macro=False, plot_micro=True)
            ax_inst.set(xlabel=None)
            ax_inst.set(ylabel=None)
            ax_inst.set(title=None)
            plt.savefig(f'dROC/E{epsilon}/S{species}_G{gamma}_E{epsilon}.eps', format='eps')
