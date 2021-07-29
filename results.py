#file for generating some graphs and results using working algorithm.py file
from algorithm import *
from Bio import SeqIO
import pickle

def get_data():
    data_list = []
    data = SeqIO.parse("data/SILVA_138.1_Fusobacteriota.fasta","fasta")
    for sample in data:
        data_list.append(str(sample.seq))
    return data_list


if __name__=='__main__':
    results_dict = {}

    data_list = get_data()

    filename = 'gamma_steps.pickle'
    epsilon = 0.0232

    prcsd_data = preprocess(data_list)
    inds, cols = filter_stable_sites(data=prcsd_data, epsilon=epsilon)
    print(f'{len(inds)} Unstable Sites Found!')

    mi_init = gen_mut_inf_mat(indices=inds, cols=cols)
    gamma = 0.70
    for step in range(6):
        mi_final, inds_final, cols_final, mv = alg(MI_list_=mi_init, gamma=gamma, indices_=inds, cols_=cols)
        results_dict[gamma] = inds_final
        gamma += 0.05


    with open(filename, 'wb') as handle:
        pickle.dump(results_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
