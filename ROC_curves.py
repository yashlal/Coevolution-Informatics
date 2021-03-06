import pandas as pd
import matplotlib.pyplot as plt
import scikitplot as skplt
import pickle

all_species = ['Fusobacteriota', 'Bacteroidota', 'Cyanobacteria']
epsilon = 0.116

gammas = [0.70,0.75,0.80,0.85,0.90,0.95]

def main():
    for species in all_species:
        with open(f'Results/AlgE{epsilon}/{species}Results_E_{epsilon}.pickle', 'rb') as handle1:
            b1 = pickle.load(handle1)
        # plt.figure(figsize=(15,7))
        # grid = plt.GridSpec(2, 3, wspace=0.3, hspace=0.4)

        for gamma in gammas:
            # n = gammas.index(gamma)
            # col, row = n%3,n//3
            # ax = plt.subplot(grid[row,col])

            true_sites = []
            for i in b1[gamma]:
                if type(i)==list:
                    for el in i:
                        true_sites.append(el)
            df = pd.read_excel(f'MFE/{species}_MFE50.xlsx', index_col=0)
            df = df.abs()
            means = df.mean()[1:].tolist()
            labels = []
            col_labels = df.columns.tolist()[1:]
            for i in col_labels:
                testval=int(i)
                if testval in true_sites:
                    labels.append(1)
                else:
                    labels.append(0)

            new_list = [[means[k], labels[k]] for k in range(len(means))]
            col_labels_2 = ['Mean MFE', 'Classification']
            df2 = pd.DataFrame(new_list, columns=col_labels_2, index=df.columns[1:])
            df2 = df2.dropna()
            df2['Mean MFE'] = df2['Mean MFE'].abs()
            df2['Mean MFE'] = df2['Mean MFE'] / max(means)
            df2 = df2.sort_values('Mean MFE', ascending=False)

            y_true = df2['Classification'].tolist()
            y_prob = df2['Mean MFE'].tolist()
            y_probas = [[1-x,x] for x in y_prob]

            ax_inst = skplt.metrics.plot_roc(y_true=y_true, y_probas=y_probas, classes_to_plot=[], plot_macro=False, plot_micro=True)
            ax_inst.set(xlabel=None)
            ax_inst.set(ylabel=None)
            ax_inst.set(title=None)
            # plt.show()
            plt.savefig(f'ROC/E{epsilon}/S{species}_G{gamma}_E{epsilon}.eps', format='eps')

if __name__=='__main__':
    main()
