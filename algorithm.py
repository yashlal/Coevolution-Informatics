from modules import *
from tqdm import tqdm
#creates a list of 3-element lists where the first two elements are the corresponding column indices and the third is the I/min value
#for the future we can merge these lists and keep the last element (index= -1) to be the sorting value
def gen_mut_inf_mat(indices, cols):
    init_u = []
    pbar2 = tqdm(range(len(indices)))
    pbar2.set_description('Doing Mutual Info Calcs')
    for i in pbar2:
        for j in range(i):
            mut_inf_ij = mutual_inf(cols[i], cols[j])
            init_u.append([[indices[i], indices[j]], mut_inf_ij])
    return sorted(init_u, key=lambda x:x[-1])

def alg(MI_list, gamma, indices, cols):
    t = 0
    while t<1:
        #join the sites
        ind_bound_site = []
        cols_bound_site = []
        for el in MI_list[-1][0]:
            i = indices.index(el)
            col = cols[i]
            ind_bound_site.append(el)
            cols_bound_site.append(col)
            indices.remove(el)
            cols.remove(col)
        indices.append(ind_bound_site)
        cols.append(cols_bound_site)

        #delete MI calcs with the old columns
        for j in MI_list.copy():
            if not set(MI_list[-1][0]).isdisjoint(j):
                MI_list.remove(j)

        #generate new MI elements
        for k in range(len(cols)-1):
            mut_inf_ = mutual_inf(cols[k], cols[-1])
            init_u.append([[indices[k], indices[-1]], mut_inf_])

        t += 1

    return indices, cols
