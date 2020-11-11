import pandas as pd
import numpy as np
from numpy import log2
import argparse
from collections import defaultdict

#parser = argparse.ArgumentParser()
#parser.add_argument("--str_input", action="store_true", help = "indicates string inputs for comparison")
#parser.add_argument("--file-input", action="store_true", help="indicates file inputs for comparison")
#parser.add_argument("--trainf", action="store_true", help="input train file name")
#parser.add_argument("--testf", action="store_true", help="input test file name")
#args = parser.parse_args()

MAXLEN = 113 # maximum length of a sequence
k = 10 # Motif width
IND_TO_NUC = dict( zip( range(4), ['A','C','G','T'] ) )
NUC_TO_IND = dict( zip( ['A','C','G','T'], range(4) ) )
def convert_nuc_to_ind(segment):
    return [NUC_TO_IND[x] for x in segment]

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



def makeCountMatrix(seq):
    """
    given a sequence, i.e. [a,t,c,g,a,a,t ... ]
    make a count matrix
    """
    l = len(seq)
    vec_a = np.zeros(l)
    vec_c = np.zeros(l)
    vec_g = np.zeros(l)
    vec_t = np.zeros(l)

    for i in range(l):
        a, c, g, t = 0, 0, 0, 0
        try:
            if seq[i] == "A":
                a += 1
            if seq[i] == "C":
                c += 1
            if seq[i] == "G":
                g += 1
            if seq[i] == "T":
                t += 1
        except:
            pass
        vec_a[i], vec_c[i], vec_g[i], vec_t[i] = a, c, g, t

    return np.stack((vec_a, vec_c, vec_g, vec_t))


def addPseudo(count_matrix, pseudo_count=(1,1,1,1)):                                                                                                                                                                                                                                                          
    pseudo_vec_a = np.zeros((MAXLEN)) + pseudo_count[0]
    pseudo_vec_c = np.zeros((MAXLEN)) + pseudo_count[1]
    pseudo_vec_g = np.zeros((MAXLEN)) + pseudo_count[2]
    pseudo_vec_t = np.zeros((MAXLEN)) + pseudo_count[3]
    pseudo_matrix = np.stack((pseudo_vec_a, pseudo_vec_c, pseudo_vec_g, pseudo_vec_t))
    return count_matrix + pseudo_matrix

def makeFrequencyMatrix(count_matrix):
    normalization = np.sum(count_matrix, axis=0)
    return count_matrix/normalization

def entropy(signal_freq_vector, background_freq):
    return log2(signal_freq_vector/background_freq)

def makeWMM(frequency_matrix, background_vec = (0.25, 0.25, 0.25, 0.25)):
    entropies = []

    for i in range(4):
        entropies.append(entropy(frequency_matrix[i], background_vec[i]))

    return np.stack(entropies)

def scanWMM(seq, motif_wmm):
    scores = []
    for i in range(MAXLEN-k+1):
        #print(seq[i:i+k])
        segment = seq[i:i+k]
        idxs = convert_nuc_to_ind(segment)
        for j in idxs:
            for n in range(i,i+k):
                score = motif_wmm[j][n]
        scores.append(score)
    return scores

def Estep():
    return None

def Mstep():
    return None
