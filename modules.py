import numpy as np

edgecases = ['.', '-', 'N', 'Y', 'M', 'S', 'K', 'R', 'W', 'V']
dict = {'A': 0, 'G':1, 'C':2, 'T':3}
l = ['A', 'G', 'C', 'T']

# entropy is shannon entropy with logbase e
def shannon_entr(l):
    sh_entr = 0
    for el in l:
        if el!=0:
            sh_entr -= el*np.log2(el)
    return sh_entr

# Calculate H(X,Y) (joint entropy) where *args is any n columns, the zip function forms a row iterator
# The idea is that again we can convert AGCTB into indices 01234 so eg if we have 3 columns then we store the probabilities in a 5,5,5 size array
# In that example, a row AAA would += the 0,0,0 element and a GCT would += 123 position
# ind_l is the list that stores the indices ie 0,0,0 or 1,2,3 to +=1 to, ar[tuple(ind_l)] += 1 is some magic that allows the indexing to work
def joint_entr(*args):
    size = tuple([5 for i in range(len(args))])
    ar = np.zeros(size)

    for row in zip(*args):
        ind_l = []
        for el in row:
            if el in edgecases:
                ind_l.append(4)
            else:
                ind_l.append(dict[el])
        ar[tuple(ind_l)] += 1
    ar = ar / (len(args[0]))
    #this line prevents log(0)= -infinity errors
    entr = np.sum(-ar*np.log2(ar, where=(ar!=0)))
    return entr
