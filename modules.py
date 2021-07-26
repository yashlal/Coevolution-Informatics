import numpy as np

bases_dict = {'A': 0, 'G':1, 'C':2, 'T':3, 'B':4, '.':4, '-':4, 'N':4, 'Y':4, 'M':4, 'S':4, 'K':4, 'R':4, 'W':4, 'V':4}


# entropy is shannon entropy with logbase 2
def shannon_entr(col):
    ar = np.zeros(5)
    for el in col:
        ar[bases_dict[el]] += 1
    ar = ar/len(col)
    entr = np.sum(-ar*np.log2(ar, where=(ar!=0)))
    return entr

# Calculate H(X,Y) (joint entropy) where *args is any n columns, the zip function forms a row iterator
# 1/len(args[0]) is the prob of 1 divided by number of rows, added each time a rows is encountered
# the try except steps the probability of a row if it exists else it adds it
def joint_entr(*args):
    prob_dict = {}
    entr = 0
    for row in zip(*args):
        s = ''.join(row)
        try:
            prob_dict[s] += 1/len(args[0])
        except KeyError:
            prob_dict[s] = 1/len(args[0])
    for el in prob_dict.values():
        entr += -el*(np.log2(el))
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
    minimum = min(H_x, H_y)
    return ((H_x+H_y-H_xy)/minimum)

# this is a bit of wizardry that takes a list or tupl with any degree of nesting (irregular) or not and flattens it
# ie list(flatten([1,2,3,[4,5,[6,[7,[8]]]]])) returns [1,2,3,4,5,6,7,8]
# is useful for the disjoint testing in the algorithm when removing old MI terms
flatten = lambda *n: (e for a in n for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))
