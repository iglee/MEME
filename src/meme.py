import pandas as pd
import numpy as np
from numpy import log2
import argparse
from collections import defaultdict, Counter
#import pickle as pkl
import matplotlib.pyplot as plt
import logomaker as lm



parser = argparse.ArgumentParser()
parser.add_argument("-tr", "--trainf", action="store", help="input train file name")
parser.add_argument("-tt", "--testf", action="store", help="input test file name")
parser.add_argument("-o", "--output-name", action="store", type=str, help="output pkl file name")
args = parser.parse_args()

MAXLEN = 113 # maximum length of a sequence
k = 10 # Motif width
IND_TO_NUC = dict( zip( range(4), ['A','C','G','T'] ) )
NUC_TO_IND = dict( zip( ['A','C','G','T'], range(4) ) )

def convert_nuc_to_ind(segment):
    return [NUC_TO_IND[x] for x in segment]

def convert_ind_to_nuc(segment):
    return [IND_TO_NUC[x] for x in segment]

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


Dtrain = read_data(args.trainf)
Dtest = read_data(args.testf)    



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
    l = count_matrix.shape[1]                                                                                                                                                                                                                                                          
    pseudo_vec_a = np.zeros((l)) + pseudo_count[0]
    pseudo_vec_c = np.zeros((l)) + pseudo_count[1]
    pseudo_vec_g = np.zeros((l)) + pseudo_count[2]
    pseudo_vec_t = np.zeros((l)) + pseudo_count[3]
    pseudo_matrix = np.stack((pseudo_vec_a, pseudo_vec_c, pseudo_vec_g, pseudo_vec_t))
    return count_matrix + pseudo_matrix

def makeFrequencyMatrix(count_matrix):
    normalization = np.sum(count_matrix, axis=0)
    return count_matrix/normalization

def entropy(wmm, freq_matrix):
    return sum(sum(wmm*freq_matrix))

def makeWMM(frequency_matrix, background_vec = (0.25, 0.25, 0.25, 0.25)):
    entropies = []

    for i in range(4):
        entropies.append(log2(frequency_matrix[i]/background_vec[i]))

    return np.stack(entropies)

def scanWMM(seq, motif_wmm):
    l = len(seq)
    scores = []
    for i in range(l-k+1):
        score = 0
        segment = seq[i:i+k]
        idxs = convert_nuc_to_ind(segment)
        for i in range(k):
            score += motif_wmm[idxs[i]][i]
        scores.append(score)
    return np.asarray(scores)

def Estep(seqs, motif_wmm):
    Ys = []
    for seq in seqs:
        scores = scanWMM(seq, motif_wmm)
        
        probs = 2**scores
        N = probs.sum()
        probs = probs/N
        Ys.append(probs)
    return np.asarray(Ys)

def Mstep(seqs, Ys):
    final_matrix = np.zeros((4,10))

    for seq, Y in zip(seqs,Ys):
        # make a Y weighted count matrix for each position
        substrings = [seq[i:i+k] for i in range(len(seq)-k+1)]
        cnt_matrix = np.array([makeCountMatrix(x) for x in substrings])
        Y_expanded = np.expand_dims(Y, axis=(1,2))
        weighted_cnt = (cnt_matrix * Y_expanded).sum(axis=0)

        # add pseudo count, with pseudo count 1
        weighted_cnt_pseudo = addPseudo(weighted_cnt)

        freq_matrix = makeFrequencyMatrix(weighted_cnt_pseudo)

        # add new WMM estimate for one sequence to final estimate
        final_matrix += makeWMM(freq_matrix)

    return final_matrix


#initialization step
seqs = list(Dtrain.values())
test_seqs = list(Dtest.values())
init_seq = seqs[0] # first string is used for seeding the algorithm
substrings = [init_seq[i:i+k] for i in range(0,len(init_seq)-k+1,k//2)]
cnt_matrices = [makeCountMatrix(x) for x in substrings]
cnt_matrices_pseudo = [addPseudo(x,(0.0625,0.0625,0.0625,0.0625)) for x in cnt_matrices]
freq_matrices = [makeFrequencyMatrix(x) for x in cnt_matrices_pseudo]
motif_wmms = [makeWMM(x) for x in freq_matrices]
init_wmms = motif_wmms.copy()

for t in range(3):
    for i in range(len(motif_wmms)):
        Ys = Estep(seqs, motif_wmms[i])
        motif_wmms[i] = Mstep(seqs, Ys)

# get max, mid, med motif_wmms from 3 runs; A, B, C
entropies = []
for x,y in zip(motif_wmms, freq_matrices):
    entropies.append(entropy(x,y))


entropies = np.asarray(entropies)
idx_sorted = entropies.argsort()[::-1]
max_idx, mid_idx, min_idx = idx_sorted[0], idx_sorted[len(idx_sorted)//2], idx_sorted[-1]

A, A_freq = motif_wmms[max_idx], freq_matrices[max_idx]
B, B_freq = motif_wmms[mid_idx], freq_matrices[mid_idx]
C, C_freq = motif_wmms[min_idx], freq_matrices[min_idx]
#print(A)

# proceed to 7 more runs
for t in range(7):
    for i in range(len(motif_wmms)):
        Ys = Estep(seqs, motif_wmms[i])
        motif_wmms[i] = Mstep(seqs, Ys)

entropies = []
for x,y in zip(motif_wmms, freq_matrices):
    entropies.append(entropy(x,y))

entropies = np.asarray(entropies)
final_idx = entropies.argmax()

D, D_freq = motif_wmms[final_idx], freq_matrices[final_idx]
#print(D)



####################################################
#             FORMAT RESULTS AND SAVE              #
####################################################


def score_wmms(test_seqs, wmm):
    scores = []
    #y_true = []
    positions = []

    # scores based on D
    for x in test_seqs:
        score = scanWMM(x, wmm)
        
        positions.append(score.argmax())
        scores.append(score)
        #y_seq = np.zeros(len(scores))
        #y_seq[46:59] = 1
        #y_true.append(y_seq)
    return scores, positions

# score motifs A, B, C, D
A_scores, A_positions = score_wmms(test_seqs, A)
B_scores, B_positions = score_wmms(test_seqs, B)
C_scores, C_positions = score_wmms(test_seqs, C)
D_scores, D_positions = score_wmms(test_seqs, D)


# plot histograms for A, B, C, D
for motif, positions in zip(["A", "B", "C", "D"],[A_positions, B_positions, C_positions, D_positions]):

    c = Counter(positions)
    plt.hist(positions, bins=30, color='#0504aa', alpha=0.7, rwidth=0.85)
    plt.text(73, 500, "Mode position: {},\n Count: {}".format(c.most_common(1)[0][0],c.most_common(1)[0][1]))
    plt.xlabel("Starting Positions", fontsize=13)  
    plt.ylabel("Counts", fontsize=13)
    plt.xticks(fontsize=12)  
    plt.yticks(fontsize=12)
    plt.title("Histogram of Starting Positions, Motif {}".format(motif), fontsize=14)
    
    plt.savefig("output/motif_{}_histogram.png".format(motif))
    plt.close()



# plot sequence logos with Logomaker
for freq_matrix, label in zip([ A_freq, B_freq, C_freq, D_freq ], ["A", "B", "C", "D"]):
    df= pd.DataFrame(freq_matrix.T, columns=['A','C','G','T'])
    lm.Logo(df)
    plt.title("Sequence Logo of Frequency Matrix, Motif {}".format(label), fontsize=14)
    plt.savefig("output/seq_logo_{}.png".format(label))
    plt.close()


#[y_true_test, scores_test] = format_results(test_seqs, D)
#[y_true_test, scores_test] = format_results(test_seqs, D)
#[y_true_test, scores_test] = format_results(test_seqs, D)
#[y_true_test, scores_test] = format_results(test_seqs, D)
#[y_true_train, scores_train] = format_results(seqs, D)


#with open("output/"+args.output_name+"_train.pkl", "wb") as f:
#    pkl.dump([y_true_train, scores_train], f)
#f.close()

#with open("output/"+args.output_name+"_test.pkl", "wb") as f:
#    pkl.dump([y_true_test, scores_test], f)
#f.close()