import numpy as np

bases_dict = {'A': 0, 'G':1, 'C':2, 'T':3, '.':4, '-':4, 'N':4, 'Y':4, 'M':4, 'S':4, 'K':4, 'R':4, 'W':4, 'V':4}


# entropy is shannon entropy with logbase 2
def shannon_entr(col):
    ar = np.zeros(5)
    for el in col:
        ar[bases_dict[el]] += 1
    ar = ar/len(col)
    entr = np.sum(-ar*np.log2(ar, where=(ar!=0)))
    return entr

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
            ind_l.append([bases_dict[el]])
        ar[tuple(ind_l)] += 1
    ar = ar / (len(args[0]))
    print(ar)
    #this line prevents log(0)= -infinity errors
    entr = np.sum(-ar*np.log2(ar, out=np.zeros_like(ar), where=(ar!=0)))
    return entr

# X and Y are lists where each element in the list is a columns eg if we had AGC and no others columns then X = [['A', 'G', 'C']]
# Even when one column is present we want the double list for preserving star operator functionality and for consistency
# Function returns mutual information normalized by division by the minimum of H_x and H_y
def mutual_inf(X, Y):
    if type(X[0])!=list:
        X=[X]
    if type(Y[0])!=list:
        Y=[Y]
    H_x = joint_entr(*X)
    H_y = joint_entr(*Y)
    H_xy = joint_entr(*X, *Y)
    print(H_x,H_y,H_xy)
    minimum = min(H_x, H_y)
    return ((H_x+H_y-H_xy)/minimum)
