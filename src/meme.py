import pandas as pd
import numpy as np
import argparse
from collections import defaultdict

#parser = argparse.ArgumentParser()
#parser.add_argument("--str_input", action="store_true", help = "indicates string inputs for comparison")
#parser.add_argument("--file-input", action="store_true", help="indicates file inputs for comparison")
#parser.add_argument("--trainf", action="store_true", help="input train file name")
#parser.add_argument("--testf", action="store_true", help="input test file name")
#args = parser.parse_args()

MAXLEN = 113

def read_data(filename):
    """
    reads files into a 
    defaultdict of identifier : list[sequence]
    """
    f = open(filename,'r')
    lines = f.readlines()
    d = defaultdict()

    for i in range(0,len(lines),3):
        ident = lines[i].strip()
        seq = lines[i+1].strip().upper() + lines[i+2].strip().upper()
        d[ident] = list(seq)

    return d


#if args.file:
#    Dtrain = read_data(args.trainf)
#    Dtest = read_data(args.testf)
    

#if args.str_input:
#    None



def makeCountMatrix(data_dict):
    vec_a = np.zeros(MAXLEN)
    vec_t = np.zeros(MAXLEN)
    vec_c = np.zeros(MAXLEN)
    vec_g = np.zeros(MAXLEN)

    for i in range(MAXLEN):
        a, t, c, g = 0, 0, 0, 0
        for seq in data_dict.values():
            try:
                if seq[i] == "A":
                    a += 1
                if seq[i] == "T":
                    t += 1
                if seq[i] == "C":
                    c += 1
                if seq[i] == "G":
                    g += 1
            except:
                pass
            vec_a[i], vec_t[i], vec_c[i], vec_g[i] = a, t, c, g

    return np.stack((vec_a, vec_t, vec_c, vec_g))


def addPseudo(count_matrix, pseudo_count=1):
    pseudo_matrix = np.zeros((4, MAXLEN)) + pseudo_count
    return count_matrix + pseudo_matrix

def makeFrequencyMatrix():
    return None

def entropy():
    return None

def makeWMM():
    return None

def scanWMM():
    return None

def Estep():
    return None

def Mstep():
    return None
