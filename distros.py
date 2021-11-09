from modules import flatten
from Bio import SeqIO
import pickle

def get_distros(data, sites, nonsites):

    i = {'A':0, 'G':1, 'C':2, 'T':3}
    c1 = [0,0,0,0,0]
    c2 = [0,0,0,0,0]
    for seq in data:
        for site in sites:
            try:
                c1[i[seq[site]]] += 1
            except:
                c1[-1] += 1

        for nonsite in nonsites:
            try:
                c2[i[seq[nonsite]]] += 1
            except:
                c2[-1] += 1

    c3 = [round(k/sum(c1), 4) for k in c1]
    c4 = [round(j/sum(c2), 4) for j in c2]

    return c3, c4

species = ['Fusobacteriota', 'Cyanobacteria', 'Bacteroidota']
epsilons = [0.116, 0.232, 1.16]
gammas = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]

for spec in species:
    print('\n')
    data = []
    path1 = SeqIO.parse(f'data/SILVA_138.1_{spec}.fasta', 'fasta')
    for record in path1:
        data.append(str(record.seq))

    for epsilon in epsilons:
        for gamma in gammas:

            if (spec=='Cyanobacteria') and ((epsilon==1.16) and (gamma==0.95)):
                continue
            path2 = f'Results/AlgE{epsilon}/{spec}Results_E_{epsilon}.pickle'
            with open(path2, 'rb') as handle:
                b = pickle.load(handle)[gamma]

            sites = list(flatten(list(filter(lambda x: type(x)==list, b))))
            nonsites = list(flatten(list(filter(lambda x: type(x)!=list, b))))

            c3,c4 = get_distros(data=data, sites=sites, nonsites=nonsites)
            print(f'Species: {spec} with Epsilon {epsilon} with Gamma {gamma}')
            print(c3,c4)
