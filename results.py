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
    species = 'Fusobacteriota'

    data_list = get_data()

    epsilon = 0.00232
    filename = f'Results/{species}Results_E_{epsilon}.pickle'

    prcsd_data = preprocess(data_list)
    inds, cols = filter_stable_sites(data=prcsd_data, epsilon=epsilon)
    print(f'{len(inds)} Unstable Sites Found!')

    mi_init = gen_mut_inf_mat(indices=inds, cols=cols)
    gamma = 0.70
    for step in range(6):
        mi_final, inds_final, cols_final, mv = alg(MI_list_=mi_init, gamma=gamma, indices_=inds, cols_=cols)
        results_dict[gamma] = inds_final
        gamma += 0.05
        plt.plot(mv)
        plt.savefig(f'Plots/Dump/S_{species}_E_{epsilon}_G_{gamma}.png')


    with open(filename, 'wb') as handle:
        pickle.dump(results_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
