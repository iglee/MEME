import pandas as pd
import numpy as np
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--str_input", action="store_true", help = "indicates string inputs for comparison")
parser.add_argument("--file-input", action="store_true", help="indicates file inputs for comparison")
args = parser.parse_args()

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
        seq = lines[i+1].strip() + lines[i+2].strip()
        d[ident] = list(seq)

    return d


if args.file:
    None

if args.str_input:
    None



def makeCountMatrix():
    return None

def addPseudo():
    return None

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
