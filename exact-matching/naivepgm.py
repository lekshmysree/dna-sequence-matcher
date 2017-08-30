#!/usr/bin/python

def readgenome(filename):
    genome=''
    with open(filename,'r') as f:
        for line in f:
            #print line
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
            #print "hiii"
            f.readline()
            qual=f.readline()
            #print seqs
            #print qual
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
        #print "the value of i={}".format(i)
        for j in range(len(p)):
            #print "t[i+j]={} p[j]={}".format(t[i+j], p[j])
            if t[i+j]!= p[j]:
                #print "not matched"
                match=False
                break
        if match:
            #print "matched"
            occurences.append(i)
    return occurences

def naive2mm(p,t):
    occurences=[]
    for i in range(len(t)-len(p)+1):
        match=True
        num_mismatch=0
        #print "the value of i={}".format(i)
        for j in range(len(p)):
            #print "t[i+j]={} p[j]={}".format(t[i+j], p[j])
            if t[i+j]!= p[j]:
                num_mismatch+=1
                if num_mismatch>2:
                   #print "more than 2 characters mismatched"
                   match=False
                   break
        if match:
            #print "matched"
            occurences.append(i)
    return occurences

seqs,quals=readfastq('ERR037900_1.first1000.fastq')
gc=gccheck(seqs)
print "gc:{}".format(gc)
#print seqs
genome=readgenome('lambda_virus.fa')
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
        