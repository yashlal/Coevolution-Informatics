#file for generating some graphs and results using working algorithm.py file
from algorithm import *
from Bio import SeqIO
import pickle

def get_data(species):
    data_list = []
    data = SeqIO.parse(f"data/SILVA_138.1_{species}.fasta","fasta")
    for sample in data:
        data_list.append(str(sample.seq))
    return data_list

if __name__=='__main__':
    results_dict = {}
    species = 'Fusobacteriota'

    print('Reading datafile...')
    data_list = get_data(species)

    epsilon = 0.0232
    filename = f'Results/{species}Results_E_{epsilon}.pickle'

    prcsd_data = preprocess(data_list)
    inds, cols = filter_stable_sites(data=prcsd_data, epsilon=epsilon)
    print(f'{len(inds)} Unstable Sites Found!')

    start=time.time()
    mi_init = mutinf_setup(indices=inds, cols=cols)
    elapsed = time.time()-start
    print(f'MI Setup took {elapsed} seconds')
    gamma = 0.70
    for step in range(6):
        mi_final, inds_final, cols_final, mv = alg(MI_list_=mi_init, gamma=gamma, indices_=inds, cols_=cols)
        results_dict[gamma] = inds_final
        plt.plot(mv)
        plt.savefig('Plots/Dump/FB_MV_E{epsilon}_G{gamma}.png')
        gamma += 0.05

    with open(filename, 'wb') as handle:
        pickle.dump(results_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
