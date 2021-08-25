import pickle
from Bio import SeqIO

blanks = ['-','.']
species = ['Cyanobacteria', "Bacteroidota"]
b = []
for spec in species:
    with open(f'Results\{spec}Results_E_0.116.pickle', 'rb') as handle:
        b.append(pickle.load(handle)[0.8])

sites = []
for i in b[0]:
    if i in b[1]:
        try:
            if len(i)>1:
                sites.append(i)
        except:
            pass

with open('data\SILVA_138.1_Cyanobacteria.fasta') as handle:
    data = list(SeqIO.parse(handle, "fasta"))
    sequence = str(data[0].seq)

site_nums = []
pure_seq_l = []
for char in sequence:
    if char not in blanks:
        pure_seq_l.append(char)
pure_seq = ''.join(pure_seq_l)
print(pure_seq)
print(len(pure_seq))
for i in range(len(sequence)):
    if sequence[i] not in blanks:
        site_nums.append(i)

new_sites = []
for site in sites:
    new_sites.append([])
    for location in site:
        new_sites[-1].append(site_nums.index(location))

print(new_sites)

with open('VP/CBBD.txt', 'w') as f:
    for item in new_sites:
        f.write("%s\n" % item)
