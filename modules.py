import numpy as np

# entropy is shannon entropy with logbase e
def shannon_entr(l):
    sh_entr = 0
    for el in l:
        if el!=0:
            sh_entr -= el*np.log(el)
    return sh_entr

#Calculate H(X,Y) (join entropy) using px*py for p(x,y)
def joint_etr(l1,l2):
    entr = 0
    for p1 in l1:
        for p2 in l2:
            if p1*p2 != 0:
                entr -= (p1*p2)*np.log(p1*p2)
    return entr

#calculate H(X|Y) using chain rule: H(X|Y) = H(X,Y)-H(Y)
def cond_entr(X, Y):
    h_xy = joint_etr(X,Y)
    h_y = shannon_entr(Y)
    return h_xy-h_y

#mutual information using I(X;Y) = H(X)-H(X|Y)
def mutual_inf(X,Y):
    h_x = shannon_entr(X)
    h_xgiveny = cond_entr(X,Y)
    return h_x-h_xgiveny
