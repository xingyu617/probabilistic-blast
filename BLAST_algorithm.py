from numpy import inf
import time
import operator
import pickle
import itertools
import pandas as pd
import re
from ast import literal_eval
import numpy as np




class BLAST:
    def __init__(self, seq="chr22.maf.ancestors.42000000.complete.boreo.fa.txt", 
    probs="chr22.maf.ancestors.42000000.complete.boreo.conf.txt", db="BoreoEutherian_11.csv",
    w=7, HSP_thresh=0.7, g=-0.5):
        self.w = w
        self.seq = self._read_seq(seq)
        self.probs = self._read_probs(probs)
        self.HSP_thresh = HSP_thresh
        self.g = g
        #self.seq = sequence           ### FOR TESTING SEQ DIRECTLY
        #self.probs = probs            ### FOR TESTING PROBS DIRECTLY
        self.database = self._read_index(db)
        self.mutation_rate = 0.5
    
    def _read_probs(self, txt):
        with open(txt, 'r') as f:
            probs = f.readline()
        probs = probs.split(' ')
        probs = probs[:-1]         #issue where last char is not a float but a space
        probs = [float(i) for i in probs] 
        return probs
    
    def _read_seq(self, txt):
        with open(txt, 'r') as f:
	        seq = f.read()
        regex = re.compile('[^a-zA-Z]')
        regex.sub('', seq)
        return seq

    def _read_index(self, index):
        index = pd.read_csv(index)
        return index

    def processing_wmers(self,query):
        df=self.database
        word_dict={}
        all_pos=[]
        count_index={}
        pos_to_query_index = {}
        
        for i in range(len(query)-self.w):
            word=query[i:i+self.w]
            row=df.loc[df['index'] == word].pos.tolist()
            row = literal_eval(row[0])
            word_dict[word]=row
            for x in row:
                all_pos.append(x)
                count_index[x]=0
                pos_to_query_index[x] = i
                

        if all_pos==[]:
            print("no w-mer is found")
            return None
        for index in all_pos:
            i=1
            while i<=(len(query)):
                if (index+i) in all_pos:
                    #print("hit",index+i)
                    count_index[index]+=1
                i+=1

        count_index = dict((k, v) for k, v in count_index.items() if v > 1)    #remove indices with 0 neighbours
        top_indices = sorted(count_index, key=count_index.get, reverse=True)

        #top_indices=top_indices[:5]      # Return top 5 indices for faster processing

        return top_indices, pos_to_query_index        
            
            

    def scanning_query(self, query):
        print("Scanning Query...", query)
        df=self.database
        score_list=[]
        best_hsp_score=-inf
        all_indices = []
       
        all_indices, pos_to_query_index = self.processing_wmers(query)

        for pos in all_indices:
            score,seq_left,seq_right,q_left,q_right = self._HSP(query, pos_to_query_index[pos], pos)
            
            if score >=self.HSP_thresh:
                sl2 = max(0, int(seq_left-(len(query)-self.w)//2))
                sr2 = min(len(self.seq), int(seq_right+(len(query)+self.w)//2))
                left_seq = self.seq[sl2:seq_left]
                right_seq = self.seq[seq_right+1:sr2]
                
                word_query = query[q_left:q_right+1]
                left_query=query[0:q_left]
                right_query=query[q_right+1:]
                
 
                align_seq_left,align_query_left,score_left,start_left,end_left=self.alignment(len(query),left_seq,left_query,0,seq_left-1,"left")
                align_seq_right,align_query_right,score_right,start_right,end_right=self.alignment(len(query),right_seq,right_query,seq_right+1,len(self.seq)-1,"right")
                sum_score=score+score_left+score_right

                #print("left",align_query_left)
                #print("left_seq:", align_seq_left)
                #print("right",align_query_right)
                #print("right_seq:", align_seq_right)
                #print("#######sum_score##########",score_left,score_right,score,sum_score)

                if score_left==1 and score_right==1:
                    best_hsp_score=sum_score
                    best_align_seq=align_seq_left+self.seq[seq_left:seq_right+1]+align_seq_right
                    best_start=start_left
                    best_end=end_right
                    best_align_query=align_query_left+query[q_left:q_right+1]+align_query_right
                    return best_hsp_score,best_align_query,best_align_seq,best_start,best_end
                if best_hsp_score<sum_score:
                    best_hsp_score=sum_score
                    best_align_seq=align_seq_left+self.seq[seq_left:seq_right+1]+align_seq_right
                    best_start=start_left
                    best_end=end_right
                    best_align_query=align_query_left+query[q_left:q_right+1]+align_query_right

        #handling exceptions
        if best_hsp_score==-inf:
            best_align_query=""
            best_align_seq=""
            best_start=-inf
            best_end=-inf
        
        #print("Score: ", best_hsp_score)
        print(best_align_query)
        print(best_align_seq)
        
        return best_hsp_score,best_align_query,best_align_seq,best_start,best_end


    def _HSP(self, query, i, pos):
        # split query up into words of size k
        # return all left and right positions 

        q_left = i
        q_right = i + self.w-1
        seq_left = pos
        seq_right = pos + self.w-1
        HSP_step = self.w
        score = []

        while q_left > 0 and seq_left > 0 and HSP_step >= 0:
            q_left -= 1
            seq_left -= 1
            if query[q_left] == self.seq[seq_left]:
                score.append(self.probs[seq_left])
            else:
                score.append((1-self.probs[seq_left])/3)
            HSP_step -= 1

        HSP_step = self.w

        while q_right < len(query)-1 and seq_right < len(self.seq)-1 and HSP_step >= 0:
            q_right += 1
            seq_right += 1
            if query[q_right] == self.seq[seq_right]:
                score.append(self.probs[seq_right])
            else:
                score.append((1-self.probs[seq_right])/3)
            HSP_step -= 1

        score=sum(score)/len(score)
        return score,seq_left,seq_right,q_left, q_right

    
    def score_gapped(self,seq1,seq2,pos1,pos2):
        match = 1
        mismatch = -1
        if seq1[pos1] == seq2[pos2]:
            return self.probs[pos1]*match + (1-self.probs[pos1])*mismatch
        else:
            return self.probs[pos1]*mismatch + (1-self.probs[pos1])/3*match

    
    def alignment(self,qlength,query, long_seq, start,end,direction):
        if len(query)==0 and direction=="left":
            return "","",0,end+1,end+1
        if len(query)==0 and direction=="right":
            return "","",0,start-1,start-1
        if len(long_seq)==0 and direction=="left":

            return "","",0,end+1,end+1
        if len(long_seq)==0 and direction=="right":
            return "","",0,start-1,start-1
        m = len(query)
        n = len(long_seq)

        matrix = np.zeros((m+1, n+1))
        for sii in range(0, m+1):
            matrix[sii][0] = sii*self.g
        if direction=="left":   
            #don't count everything after the last nuc of query
            for sjj in range(0, n+1):
                matrix[0][sjj] = sjj*self.g
        else:
            for sjj in range(0, n+1):
                matrix[0][sjj] = 0

        for siii in range(1, m+1):
            for sjjj in range(1, n+1):
                matrix[siii][sjjj] = max(matrix[siii-1][sjjj] + self.g, matrix[siii - 1][sjjj - 1] + self.score_gapped(query,long_seq,siii-1, sjjj-1), matrix[siii][sjjj-1] + self.g)
 
        sequ1 = []
        sequ2 = []
        #select the start point to trace back
        if direction=="right":
            res = np.where(matrix[-1] == np.amax(matrix[-1]))
            n = res[0][0]
        
        while m > 0 and n > 0:
            #print("m", m,"n", n,"rows", len(matrix),"cols", len(matrix[0]))
            if (matrix[m-1][n-1]>=matrix[m-1][n] and matrix[m-1][n-1]>=matrix[m][n-1]) or query[m-1] == long_seq[n-1]:
                sequ1.append(query[m-1])
                sequ2.append(long_seq[n-1])
                m -= 1
                n -= 1
            elif matrix[m][n-1]>=matrix[m-1][n-1] and matrix[m][n-1]>=matrix[m-1][n]:
                #print("add gap 1")
                sequ1.append('-')
                sequ2.append(long_seq[n-1])
                n -= 1
            else:
                sequ1.append(query[m-1])
                #print("add gap 2")
                sequ2.append('-')
                m -= 1

        if direction=="right":
            seq_start=start
            seq_end=np.argmax(matrix[n-1])-1+start
            #seq_end=matrix[n-1].index(max(matrix[n-1]))-1+start
            if m>n:
                for i in range(m-n):
                    sequ1.append(query[m-n-1-i])
                    sequ2.append('-')
                    #print("add gap")
        else:
            seq_start=end-len(query)
            seq_end=end
        '''
        if m>n:
            for i in range(m-n):
                sequ1.append(seq1[m-n-1-i])
                sequ2.append('-')
        elif n>m:
            for i in range(n-m):
                sequ2.append(seq2[n-m-1-i])
                sequ1.append('-')
        '''

        sequ1.reverse()
        sequ2.reverse()
        result1 = ''.join(sequ1)
        result2 = ''.join(sequ2)
        align_score = 0.
        for k in range(0, len(result1)):
            if result1[k] == result2[k]:
                align_score += 1

        align_score = float(align_score)/qlength

        return result1, result2, align_score,seq_start,seq_end

if __name__ == "__main__":
    blast = BLAST()