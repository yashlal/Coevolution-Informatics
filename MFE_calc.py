from subprocess import PIPE, Popen
import random as rd
import tqdm
import pickle

bases = ['A','G','C','T']

def process_val(raw_input):
    val = []
    s=0
    for char in raw_input:
        if char=='-':
            s=1
        if s:
            if char!=')':
                val.append(char)
            else:
                s=0
    val = float(''.join(val))
    return val

def mutate(pxn, full_seq):
    new_s = full_seq.copy()

    new_l = bases.copy()
    new_l.remove(new_s[pxn])

    rand = rd.randint(0,2)
    new_s[pxn] = new_l[rand]

    return new_s

def main_MFE_func(s, pxns):
    ref_input = ''.join(s)
    p = Popen('C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe', stdin=PIPE, stdout=PIPE)
    reference = list(p.communicate(ref_input.encode())[0].decode())
    ref_val = process_val(reference)

    MFE_sum = 0

    pbar = tqdm.tqdm(range(len(pxns)))
    for dummy_iter in pbar:

        x = pxns[dummy_iter]
        inp = ''.join(mutate(x, s))

        p = Popen('C:\Program Files (x86)\ViennaRNA Package\RNAfold.exe', stdin=PIPE, stdout=PIPE)
        ans = list(p.communicate(inp.encode())[0].decode())

        value_for_dict = process_val(ans)
        MFE_sum += value_for_dict-ref_val

    avg_mfe_pert = MFE_sum/len(pxns)
    return avg_mfe_pert
