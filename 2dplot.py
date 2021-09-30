import pandas as pd
import matplotlib.pyplot as plt
import scikitplot as skplt
import pickle
import pandas as pd

species = 'Cyanobacteria'
with open(f'Results/{species}Results_E_0.116.pickle', 'rb') as handle1:
    b1 = pickle.load(handle1)[0.75]

df = pd.read_excel(f'MFE/{species}_MFE50.xlsx', index_col=0)
df = df.abs()
means = df.mean()[1:].tolist()
inds_list = df.columns[1:].tolist()
inds_list = [int(k) for k in inds_list]
sites = list(filter(lambda elm: isinstance(elm, list), b1))

mfe_vals = []
sep_vals = []


for site in sites:
    mfe_val = 0
    sep_val = 0
    l = len(site)
    l2 = int(l*(l-1)/2)
    for el in site:
        # try:
        site_index = inds_list.index(el)
        print(site_index)
        mfe_val += means[site_index]
        for el2 in site[1:]:
            sep_val += abs(el2-el)
        mfe_val = mfe_val / l
        sep_val = sep_val / l2

        if sep_val >= 25000:
            print(site)

        mfe_vals.append(mfe_val)
        sep_vals.append(sep_val)


print(max(sep_vals))

plt.scatter(mfe_vals, sep_vals)
plt.xlabel('Average MFE per Mutation in Site')
plt.ylabel('Average Pairwise Distance in Site')
plt.title('MFE vs. Distance (Average) in CB sites at Gamma 0.95')
plt.show()
