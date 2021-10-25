import pickle
import modules

def hamm_dist(x1,x2):
    dist = 0
    count = 0
    for i in range(len(x1)):
        for j in range(len(x1[i])):
            count +=1
            if x1[i][j]!=x2[i][j]:
                dist += 1
    return (dist, count, dist/count)

species_list = ['Cyanobacteria']
epsilon = 0.116
gammas = [0.95, 0.90, 0.85]

x_sets = ['Old', 'New']
y_sets = ['R1', 'R2', 'R3', 'R4']

for x_set in x_sets:
    for y_set in y_sets:
        print('\n')
        for species in species_list:
            for gamma in gammas:
                with open(f'Results/MFERankings/{x_set}/{species}_E{epsilon}_G{gamma}_Ranking.pickle', 'rb') as handle1:
                    b1 = pickle.load(handle1)
                    b1 = list(filter(lambda x: type(x)==list, b1))
                    # for i in b1:
                    #     i.reverse()

                with open(f'Results/InfRankings/{y_set}/{species}_E{epsilon}_G{gamma}_Ranking.pickle', 'rb') as handle2:
                    b2 = pickle.load(handle2)
                    # for j in b2:
                    #     j.reverse()

                tup1 = hamm_dist(b1,b2)
                print(species, gamma, tup1)
