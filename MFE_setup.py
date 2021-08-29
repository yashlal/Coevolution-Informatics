import pickle
from Bio import SeqIO
import modules
import MFE_calc
import pandas as pd

blanks = ['-','.']
species = ['Cyanobacteria']
b = []
for spec in species:
    with open(f'Results\{spec}Results_E_0.116.pickle', 'rb') as handle:
        b.append(pickle.load(handle)[0.8])

sites = []
nonsites = []

for i in b[0]:
    if type(i)==list:
        sites.append(i)
    else:
        nonsites.append(i)

sites=list(modules.flatten(sites))

def run_sequence(seq):
    site_nums = []
    pure_seq_l = []

    for i in range(len(seq)):
        if seq[i] not in blanks:
            site_nums.append(i)
            pure_seq_l.append(seq[i])

    sites_2 = []
    nonsites_2 = []
    for site in sites:
        try:
            sites_2.append(site_nums.index(site))
        except:
            sites_2.append('BLANK')
    for nonsite in nonsites:
        try:
            nonsites_2.append(site_nums.index(nonsite))
        except:
            nonsites_2.append('BLANK')

    MFE_pos, rv = MFE_calc.MFE_func_all(pure_seq_l, sites_2)
    MFE_neg, rv = MFE_calc.MFE_func_all(pure_seq_l, nonsites_2)

    return MFE_pos, MFE_neg, rv

with open('data\SILVA_138.1_Cyanobacteria.fasta') as handle:
    data = list(SeqIO.parse(handle, "fasta"))

bad_seqs = [51, 74, 75, 76, 93, 113, 118, 120]

if __name__=='__main__':
    results_l = []
    index_col = []

    for j in range(50):
        if j not in bad_seqs:
            index_col.append(j)
            print(f'Sequence {j}/50')
            sequence = str(data[j].seq)
            p_vals, n_vals, ref = run_sequence(sequence)
            print(p_vals, n_vals, ref)
            results_l.append([ref]+p_vals+n_vals)
    columns_labels = ['REFERENCE']+[str(_) for _ in sites]+[str(n_) for n_ in nonsites]
    results_df = pd.DataFrame(results_l, index=index_col, columns=columns_labels)

    results_df.to_excel('CB_MFE50.xlsx')
    print(results_df)
