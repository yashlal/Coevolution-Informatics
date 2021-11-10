import distanceROC
from itertools import combinations
import random as rd
from Bio import SeqIO
import pickle
import modules
import mendotaMFE
import RNA
import numpy as np

blanks = ['-','.']

gamma = 0.85
epsilon = 0.232
species_list = ['Fusobacteriota', 'Bacteroidota', 'Cyanobacteria']

if __name__=='__main__':
    coords, ref_l = distanceROC.get_reference_dist()
    stats = []
    for species in species_list:

        stats.append([])

        data_list = []
        data = SeqIO.parse(f'data/SILVA_138.1_{species}.fasta',"fasta")
        for sample in data:
            data_list.append(str(sample.seq))
        seq = list(data_list[0])

        with open(f'Results/AlgE{epsilon}/{species}Results_E_{epsilon}.pickle', 'rb') as handle:
            b = pickle.load(handle)[gamma]

        c = b.copy()
        all_sites = list(modules.flatten(c))

        comb = list(combinations(all_sites,2))
        rd.shuffle(comb)
        rd.shuffle(comb)
        rd.shuffle(comb)

        tups = []
        for tup in comb:
            i,j=tup
            if (seq[i] not in blanks) and (seq[j] not in blanks):
                if (i in ref_l) and (j in ref_l):
                    if len(tups)<400:
                        tups.append(tup)
                    else:
                        break

        for pair_ind in range(len(tups)):
            n_tot = len(tups)
            print(f"{pair_ind} / {n_tot}")
            i,j = tups[pair_ind]

            seqinp = ''.join(list(filter(lambda x: x not in blanks, seq)))
            inp1 = ''.join(list(filter(lambda x: x not in blanks, mendotaMFE.mutate(i, seq))))
            inp2 = ''.join(list(filter(lambda x: x not in blanks, mendotaMFE.mutate(j, seq))))

            mfe1 = abs(RNA.fold(inp1)[1] - RNA.fold(seqinp)[1])
            mfe2 = abs(RNA.fold(inp2)[1] - RNA.fold(seqinp)[1])
            mfe_stat = (mfe1 + mfe2)/2

            coord1 = coords[ref_l.index(i)]
            coord2 = coords[ref_l.index(j)]
            dist_stat = np.linalg.norm(coord1-coord2)

            stats[-1].append((mfe_stat, dist_stat))

            with open('Results/distance_mfe_corr.pickle', 'wb') as handle2:
                pickle.dump(stats, handle2)
