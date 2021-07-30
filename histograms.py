import pickle
import matplotlib.pyplot as plt
import modules

with open('gamma_steps.pickle', 'rb') as handle:
    b = pickle.load(handle)

def algplots():
    list_n_sites = []
    all_bound_cols = {}

    for x in range(1,6):
        mytuple = list(b.items())[x]
        sites = [_ for _ in mytuple[1] if type(_)==list]
        gamma = round(mytuple[0], 2)

        n_sites = len(sites)
        list_n_sites.append(n_sites)

        all_bound_cols[gamma] = list(modules.flatten(sites))

        hist_dict = {}
        for site in sites:
            hist_dict[len(site)] = hist_dict.get(len(site), 0) + 1

        print(gamma, hist_dict)

        plt.bar(list(hist_dict.keys()), list(hist_dict.values()), width=0.8, color='g')
        plt.title(f'Gamma={gamma} yields {n_sites} Sites')
        plt.xlabel('Size of Sites')
        plt.ylabel('Number of Sites')
        plt.savefig(f'Plots/AlgGamma{gamma}.png')
        plt.clf()

    return list_n_sites, list(all_bound_cols.values())

if __name__=='__main__':
    itervar = algplots()
