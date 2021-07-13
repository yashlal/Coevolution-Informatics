from Bio import SeqIO
import pandas as pd
import numpy as np


data_list = []
data = SeqIO.parse("SILVA_138.1_Fusobacteriota.fasta","fasta")
for sample in data:
    data_list.append(sample.seq)

for column_iter in range(len(data_list[0])):
    col = [x[column_iter] for x in data_list]
    print(col)
