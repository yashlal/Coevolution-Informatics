from Bio import PDB
from Bio import SeqIO
from Bio.PDB.Polypeptide import PPBuilder

parser = PDB.MMCIFParser()
structure = parser.get_structure('4ybb', '4ybb.cif')
model = structure[0]
chain = model['AA']

reds = [r for r in chain]


for i in range(len(reds)):
    try:
        print(reds[i]["C5'"].get_coord())
    except:
        print(i)
