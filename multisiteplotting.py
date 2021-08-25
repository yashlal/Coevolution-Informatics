import pickle

with open('Results/CBFB.pickle', 'rb') as handle:
    f = pickle.load(handle)

dists = []
for site in f:
    dist = 0
    c = 0
    for k,v in site[0].items():
        try:
            dist += abs(site[1][k]-v)
            c += 1
        except:
            dist += site[0][k]
            c += 1
    for k,v in site[1].items():
        if k in site[0].keys():
            pass
        else:
            dist += site[1][k]
            c += 1
    dist = dist / c
    dists.append(dist)

print(dists)
