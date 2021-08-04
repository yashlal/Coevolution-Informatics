from Bio import SeqIO
import pickle

path = 'data/SILVA_138.1_Fusobacteriota_plus_Ecoli.fasta'
inds = [1012, 1011, 1020, 1019, 1008, 1007, 1013, 1010]
bases = ['A', 'G', 'C', 'T']

data = SeqIO.parse(path,"fasta")
data_list = []
for sample in data:
    data_list.append(str(sample.seq))

ecoli = []
l = list(data_list[-1])

for i in range(len(l)):
    if l[i] in bases:
        ecoli.append(i)

print(f'Site of Size {len(inds)}')

new_inds = sorted([ecoli.index(ind) for ind in inds])

print(new_inds)
