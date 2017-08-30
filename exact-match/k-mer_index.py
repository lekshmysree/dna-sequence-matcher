#!/usr/bin/env python

"""kmer_index.py: A k-mer index for indexing a text."""

import bisect


class Index(object):
    """ Holds a substring index for a text T """

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
        
 
    def queryindex(self,p,t,index):
        k=index.k
        offsets=[]
        for i in index.query(p):
            if p[k:]==t[(i+k):(i+len(p))]:
                offsets.append(i)
        return offsets
        
    def appmatch(self,p,t,n):
        seq_length=int(round(len(p)/(n+1)))
        all_match=set()
        for i in range(n+1):
            start=i*seq_length
            end=min((i+1)*seq_length,len(p))
            #index=Index(t,k)
            match=self.queryindex(p[start:end],t,index)
            for m in match:
                if m<start or m-start+len(p) > len(t):
                    continue
                mismatch=0
                for j in range(0,start):
                    if p[j]!=t[m-start+j]:
                        mismatch+=1
                        if mismatch>n:
                            break
                for j in range(end,len(p)):
                    if p[j]!=t[m-start+j]:
                        mismatch+=1
                        if mismatch>n:
                            break
                if mismatch<=n:
                    all_match.add(m-start)
        return list(all_match)
        
        
        

def readgenome(filename):
    genome=''
    with open(filename,'r') as f:
        for line in f:
            if line[0]!='>':
                genome+=line.rstrip()
    return genome

p='GGCGCGGTGGCTCACGCCTGTAAT'
genome=readgenome('resources/chr1.GRCh38.excerpt.fasta')
index=Index(genome,8)
print len(index.query(p))
print index.queryindex(p,genome,index)
#print(index.appmatch(p,genome,2))
#print len(index.appmatch(p,genome,2))
