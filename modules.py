import numpy as np
from collections import Counter
import time

bases_dict = {'A': 0, 'G':1, 'C':2, 'T':3, 'B':4, '.':4, '-':4, 'N':4, 'Y':4, 'M':4, 'S':4, 'K':4, 'R':4, 'W':4, 'V':4, 'D':4, 'H':4}


# entropy is shannon entropy with logbase 2
def shannon_entr(col):
    l = [0,0,0,0,0]
    itr_amnt = 1/len(col)
    for x in col:
        l[bases_dict[x]] += itr_amnt
    entr = sum([-v*np.log2(v) for v in l if v])
    return entr

# Calculate H(X,Y) (joint entropy) where *args is any n columns, the zip function forms a row iterator
# 1/len(args[0]) is the prob of 1 divided by number of rows, added each time a rows is encountered
# the try except steps the probability of a row if it exists else it adds it
def joint_entr(*args):
    prob_dict = {}
    prob_dict = {row: prob_dict.get(row,0) + (1/len(args[0])) for row in zip(*args)}

    probs_ar = np.array(list(prob_dict.values()))
    entr = sum(-probs_ar*(np.log2(probs_ar, out=np.zeros_like(probs_ar), where=(probs_ar!=0))))
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

# this is a bit of StackOverflow wizardry that takes a list or tuple with any degree of nesting and flattens it
# ie list(flatten([1,2,3,[4,5,[6,[7,[8]]]]])) returns [1,2,3,4,5,6,7,8]
# is useful for the disjoint testing in the algorithm when removing old MI terms
# Credit to samplebias from https://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists
flatten = lambda *n: (e for a in n for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))

# a form of the MI function that allows multiprocessing
def mutual_inf_MP(X, Y, ind1, ind2):
    if type(X[0])!=list:
        X=[X]
    if type(Y[0])!=list:
        Y=[Y]
    H_x = joint_entr(*X)
    H_y = joint_entr(*Y)
    H_xy = joint_entr(*X, *Y)
    minimum = min(H_x, H_y)
    calc = ((H_x+H_y-H_xy)/minimum)

    return [ind1, ind2, calc]

def inv_across_envs(*args):
    cnt = Counter(flatten(args))
    shared_items = [k for k, v in cnt.items() if v > 1]
    return shared_items
