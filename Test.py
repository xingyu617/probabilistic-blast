from numpy import inf
import time
import operator
import pickle
import numpy as np
import itertools
import pandas as pd
import re
import matplotlib.pyplot as plt
import random
from random import randrange
from ast import literal_eval

import BLAST_algorithm as BLAST_algorithm



class Test:
    def __init__(self,blast):
        self.prob_list=blast.probs
        self.seq_list=blast.seq
        self.blast = blast

    def generate_distribution(self,index):
        seq=self.seq_list
        prob=self.prob_list
        #index=self.index
        nucs=np.array(["A",'G','C','T'])
        max_nuc=np.where(nucs==seq[index])
        min_prob=(1-prob[index])/3
        distribution=np.full(4,min_prob)
        distribution[max_nuc]=prob[index]
        #print(nucs,distribution)
        return nucs,distribution

    def generate_random_query(self,length,index,error_rate):
        seq=self.seq_list
        prob=self.prob_list

        seq_length=len(seq)
        query=[]
        
        for i in range(length):
            nucs,distribution=self.generate_distribution(index+i)
            nuc=np.random.choice(nucs, p=distribution)
            query.append(nuc)
            #print(nuc)
        query="".join(query)
        #print("before permutation",query)
        if (error_rate == 0):
            return query, index
        query=self.permute_sequence(query,error_rate)

        #print("after permutation",query)
        return query,index


    def permute_sequence(self,sequence,error_rate):
        nucs=["A","G","C","T"]
        result = ""
        for s in sequence:
            if random.random() <= error_rate:
                if random.randint(0,1):
                        # delete
                    continue
                else:
                        # insert
                    result += nucs[random.randint(0,3)]
            else:
                result += s

        return result
    def run(self,iteration, length=[50, 100, 200, 400, 800, 1600, 3200, 5000], error_rate=[0.00,0.03, 0.05, 0.07, 0.1, 0.2]):

        hit_array=np.full((len(length),len(error_rate)),0.0)
        time_array=np.full((len(length),len(error_rate)),0.0)
        
        for i in range(len(length)):
            for k in range(len(error_rate)):
                hit=0
                total_time=0
                for j in range(iteration):
                    index=randrange(len(self.seq_list)-length[i])
                    seq,real_start_pos = self.generate_random_query(length[i],index,error_rate[k])
                    seq_len = len(seq)
                    shift=abs(seq_len-length[i])

                    #print("generate random sequence at", real_start_pos,shift)
                    end_pos = real_start_pos + seq_len

                    start_time = time.time()

                    score,align1,align2,start_pos,end_pos=self.blast.scanning_query(seq)
                    diff_time = time.time() - start_time
                    #print("TIME",time.time(),start_time,diff_time)

                    #print("Actual Position", real_start_pos)
                    #print("Returned Position", start_pos)
                    if score == -inf:
                        hit+=0
                    elif real_start_pos-shift<=start_pos<=real_start_pos+shift:
                        hit+=1
                    
                    total_time+=diff_time

                average_hit=float(hit)/iteration
                print("HITS: ",hit, "/", iteration, "AVERAGE HIT: ", average_hit)
                average_time=float(total_time)/iteration
                hit_array[i][k]=average_hit
                time_array[i][k]=average_time

        # Returns n x m array, where n = different length of sequences, and m = different error rates
        return hit_array,time_array



def test_blast_algorithm():
    # Testing Scanning Query
    START_POS = 1000
    QUERY_LENGTH = 100
    MUTATION_RATE = 0.05

    blast = BLAST_algorithm.BLAST("chr22.maf.ancestors.42000000.complete.boreo.fa.txt", 
    "chr22.maf.ancestors.42000000.complete.boreo.conf.txt", 
    db="BoreoEutherian_9.csv",
    w=9, thresh2=0.7, g=-0.75)
    test=Test(blast)

    query, index=test.generate_random_query(QUERY_LENGTH,START_POS,MUTATION_RATE)
    query2 = blast.scanning_query(query)

    print("query", query)
    print(query2[0])
    print(query2[1])
    print(query2[2])
    print("Start", query2[3])
    print("End", query2[4])
    # The "Start" position stored in query2[3] should match the START_POS value
    # The "End" position sotred in query2[4] should equal START_POS + QUERY_LENGTH


def test_accuracy():
    # Generate graphs for accuracy of different sequence lengths
    w = 9          #word length
    HSP_thresh = 0.7
    g = -0.75       #linear gap cost

    ERROR_RATE = [0, 0.05, 0.1, 0.2]
    NUM_ITERATIONS = 50
    SEQUENCE_LENGTH = [[50], [100], [200]]


    sequence = "chr22.maf.ancestors.42000000.complete.boreo.fa.txt"
    probs = "chr22.maf.ancestors.42000000.complete.boreo.conf.txt"
    db = {7:"BoreoEutherian_7.csv", 11:"BoreoEutherian_11.csv", 9:"BoreoEutherian_9.csv"}

    #print("PARAMETERS: ", w, HSP_thresh, g)
    testblast = BLAST_algorithm.BLAST(sequence, probs, db[w], w, HSP_thresh, g)
    test = Test(testblast)
    hit_array_50, time_array_50 = test.run(NUM_ITERATIONS, SEQUENCE_LENGTH[0], ERROR_RATE)
    hit_array_100, time_array_100 = test.run(NUM_ITERATIONS, SEQUENCE_LENGTH[1], ERROR_RATE)
    hit_array_200, time_array_200 = test.run(NUM_ITERATIONS, SEQUENCE_LENGTH[2], ERROR_RATE)

    return hit_array_50[0], hit_array_100[0], hit_array_200[0]


def plot_bar_graphs(hit_array_50, hit_array_100, hit_array_200):
    labels = ['0.00', '0.05', '0.10', '0.20']

    x = np.arange(len(labels))  # the label locations
    width = 0.4  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x-width , hit_array_50, width/2, label='50')
    rects2 = ax.bar(x-width/2 , hit_array_100, width/2, label='100')
    rects3= ax.bar(x , hit_array_200, width/2, label='200')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Accuracy')
    ax.set_xlabel("Error rate")
    ax.set_title('Accuracy by query length and error rate')
    ax.set_xticklabels(labels)
    ax.set_xticks([-0.20, 0.80, 1.80, 2.80])
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.legend(bbox_to_anchor=(1.05, 1.05))

    fig.tight_layout()

    plt.savefig("test_accuracy.png")
    plt.show()


if __name__ == "__main__":
    test_blast_algorithm()
    plot_bar_graphs(test_accuracy())