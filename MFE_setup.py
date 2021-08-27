import pickle
from Bio import SeqIO
import modules
import MFE_calc

blanks = ['-','.']
species = ['Cyanobacteria', "Bacteroidota"]
b = []
for spec in species:
    with open(f'Results\{spec}Results_E_0.116.pickle', 'rb') as handle:
        b.append(pickle.load(handle)[0.8])

sites = []
nonsites = []
for i in b[0]:
    if i in b[1]:
        try:
            if len(i)>1:
                sites.append(i)
        except:
            pass
    else:
        if type(i)!=list:
            nonsites.append(i)

sites=list(modules.flatten(sites))

with open('data\SILVA_138.1_Cyanobacteria.fasta') as handle:
    data = list(SeqIO.parse(handle, "fasta"))
    sequence = str(data[0].seq)

site_nums = []
pure_seq_l = []

for i in range(len(sequence)):
    if sequence[i] not in blanks:
        site_nums.append(i)
        pure_seq_l.append(sequence[i])

sites_2 = []
nonsites_2 = []
for site in sites:
    sites_2.append(site_nums.index(site))
for nonsite in nonsites:
    try:
        nonsites_2.append(site_nums.index(nonsite))
    except:
        pass

MFE_pos = MFE_calc.main_MFE_func(pure_seq_l, sites_2)
MFE_neg = MFE_calc.main_MFE_func(pure_seq_l, nonsites_2)

print(MFE_pos, MFE_neg)
