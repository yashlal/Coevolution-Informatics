from Bio import SeqIO
import pickle
import json

mypath = 'data/SILVA_138.1_Fusobacteriota_plus_Ecoli.fasta'
bases = ['A', 'G', 'C', 'T']

def ecoli_setup(path):
    data = SeqIO.parse(path,"fasta")
    data_list = []
    for sample in data:
        data_list.append(str(sample.seq))

    ecoli = []
    l = list(data_list[-1])

    for i in range(len(l)):
        if l[i] in bases:
            ecoli.append(i)
    return ecoli

def get_plot_inds(ecoli_l, sites):
    output = []

    for arg in sites:
        print(f'Site of Size {len(arg)}')

        try:
            output.append(sorted([ecoli_l.index(ind)+1 for ind in arg]))
        except ValueError:
            print('Site Error')
            continue

    return output

# with open('Results\FusobacteriotaResults_E_0.116.pickle', 'rb') as handle:
#     b = pickle.load(handle)

all_sites = [[10353, 13882], [40955, 41445], [22567, 25431], [25509, 26976], [25499, 26983], [10356, 13881], [41502, 41533], [15666, 22460], [8358, 9819], [26161, 26802], [21773, 21934], [40730, 40779], [11884, 13878], [27173, 27638], [21771, 21966], [21335, 22088], [2071, 5264], [21767, 21975], [8417, 8599], [1045, 1047], [9884, 10260], [40318, 40871], [5289, 5459], [31189, 40147], [43104, 43106], [40962, 41433], [40726, 40782], [2083, 5248], [21294, 22112], [6208, 6325], [31195, 40140], [21340, 22084], [9883, 10262], [5342, 5408], [2520, 3163], [9880, 10272]]

# for x in b.values():
#     for i in x:
#         if type(i)!=int:
#             if len(i)>2:
#                 all_sites.append(i)

ecoli_l = ecoli_setup(mypath)
all_proper_sites = get_plot_inds(ecoli_l=ecoli_l, sites=all_sites)
print(all_proper_sites)

with open('Results/FB_BT.txt', 'w') as f:
    for item in all_proper_sites:
        f.write("%s\n" % item)
