import itertools
import pandas as pd
import re
from ast import literal_eval
import numpy as np
import time
import operator


def generate_index(seq, probs, w=7, mutation_rate=0.5):
    x = 'ACTG'
    indices = [p for p in itertools.product(x, repeat=w)]
    indices = [''.join(t) for t in indices]

    df = {'index':indices, "pos":[[] for _ in range(len(indices))]}
    df = pd.DataFrame(df)
    df = df.set_index('index')

    for i in range(len(seq)-w):
        if (i%1000 == 0):
            print("Current: ",i)
        word = seq[i:i+w]
        probability = 1 
        #library[word].append(i) 
        df.loc[word,'pos'].append(i)

        for j in range(i, i+w):
            '''
            probability = probability*self.probs[j]
        # if the probability is over the treshold, then add it to the index
        if probability >= self.thresh1:
            df.loc[word, 'pos'].append(i)
        '''
            if probs[j]<mutation_rate:
                w1=seq[i:j]+"A"+seq[j+1:i+w]
                w2=seq[i:j]+"G"+seq[j+1:i+w]
                w3=seq[i:j]+"C"+seq[j+1:i+w]
                w4=seq[i:j]+"T"+seq[j+1:i+w]
                #print(word[:j],j,word[j+1:],len(w1),len(word))
                w_list=[w1,w2,w3,w4]
                #print(w_list)
                #print(word)
                for ele in w_list:
                    if i not in df.loc[ele,'pos']:
                        df.loc[ele,'pos'].append(i)

    df.to_csv('BoreoEutherian_'+ str(w) + '.csv', index=True)
    return df

if __name__ == "__main__":
    w = 9
    mutation_rate = 0.5      #default

    def read_probs(txt):
        with open(txt, 'r') as f:
            probs = f.readline()
        probs = probs.split(' ')
        probs = probs[:-1]         #issue where last char is not a float but a space
        probs = [float(i) for i in probs] 
        return probs

    def read_seq(txt):
        with open(txt, 'r') as f:
            seq = f.read()
        regex = re.compile('[^a-zA-Z]')
        regex.sub('', seq)
        return seq

    seq = read_seq("chr22.maf.ancestors.42000000.complete.boreo.fa.txt")
    probs = read_probs("chr22.maf.ancestors.42000000.complete.boreo.conf.txt")

    lib = generate_index(seq, probs, w, mutation_rate)
    print(lib.head())