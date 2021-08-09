from Bio import SeqIO
import pickle
import json

mypath = 'data/SILVA_138.1_Fusobacteriota_plus_Ecoli.fasta'
bases = ['A', 'G', 'C', 'T']

def ecoli_setup(path):
    data = SeqIO.parse(path,"fasta")
    data_list = []
    for sample in data:
        data_list.append(str(sample.seq))

    ecoli = []
    l = list(data_list[-1])

    for i in range(len(l)):
        if l[i] in bases:
            ecoli.append(i)
    return ecoli

def get_plot_inds(ecoli_l, sites):
    output = []

    for arg in sites:
        print(f'Site of Size {len(arg)}')

        try:
            output.append(sorted([ecoli_l.index(ind)+1 for ind in arg]))
        except ValueError:
            print('Site Error')
            continue

    return output

with open('Results\FusobacteriotaResults_E_0.116.pickle', 'rb') as handle:
    b = pickle.load(handle)

all_sites = []

for x in b.values():
    for i in x:
        if type(i)!=int:
            if len(i)>2:
                all_sites.append(i)

ecoli_l = ecoli_setup(mypath)
all_proper_sites = get_plot_inds(ecoli_l=ecoli_l, sites=all_sites)
print(all_proper_sites)

with open('Results/FB_VP.txt', 'w') as f:
    for item in all_proper_sites:
        f.write("%s\n" % item)
