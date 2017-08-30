#!/usr/bin/python

def readgenome(filename):
    genome=''
    with open(filename,'r') as f:
        for line in f:
            if line[0]!='>':
                genome+=line.rstrip()
    return genome

def readfastq(filename):
    sequences=[]
    qualities=[]
    with open(filename,'r') as f:
        while True:
            f.readline()
            seqs=f.readline()
            f.readline()
            qual=f.readline()
            if len(seqs)==0:
                break
            sequences.append(seqs)
            qualities.append(qual)
    return sequences,qualities
        
def gccheck(reads):
    gc=[0]*101
    total=[0]*101
    for read in reads:
        for i in range(len(read)):
            if read[i]=='G' or read[i]=='C':
                gc[i]+=1
            total[i]+=1
    for i in range(len(gc)):
        if total[i]!=0:
            gc[i]/=float(total[i])
    return gc
        

def reversecomplement(p):
    basecomplement={'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    t=''
    for base in p:
        t=basecomplement[base]+t
    return t


def naive(p,t):
    occurences=[]
    for i in range(len(t)-len(p)+1):
        match=True
        for j in range(len(p)):
            if t[i+j]!= p[j]:
                match=False
                break
        if match:
            occurences.append(i)
    return occurences

def naive2mm(p,t):
    occurences=[]
    for i in range(len(t)-len(p)+1):
        match=True
        num_mismatch=0
        for j in range(len(p)):
            if t[i+j]!= p[j]:
                num_mismatch+=1
                if num_mismatch>2:
                   match=False
                   break
        if match:
            occurences.append(i)
    return occurences

seqs,quals=readfastq('resources/ERR037900_1.first1000.fastq')
gc=gccheck(seqs)
print "gc:{}".format(gc)
#print seqs
genome=readgenome('resources/lambda_virus.fa')
#with open('res.txt','w') as f:
    #f.write(genome)

#genome=readgenome('FA1.txt')

pattern='AGGAGGTT'
t=reversecomplement(pattern)
occurences=naive2mm(pattern,genome)

#print "pattern={}   t={}".format(repr(pattern),repr(t))
#if pattern != t:
    #reverseoccurence=naive(t,genome)
    #print "reverseoccurences is :{}".format(reverseoccurence)
    #occurences.extend(reverseoccurence)
    #print "checked the complementary strand"

   
#print "the number of matches occured is {}".format(len(occurences))
#print occurences
#print len(occurences)
#print genome
#print t
        