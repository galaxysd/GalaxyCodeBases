#!/usr/bin/env python3

#InFile = 'out.kmers.1k'

#import csv
import primer3
import itertools
import re

from pprint import pprint

import math
from collections import Counter

myKmerLen = 18
myEntropyKmer = myKmerLen - 3
myMinEntropy = myKmerLen * 3

def entropy(seq, max_k=5):
    #p, lns = Counter(s), float(len(s))
    #return -sum( count/lns * math.log(count/lns, 2) for count in p.values())
    SeqLength = len(seq)
    Kmers = [ [] for x in range(max_k) ]
    for i in range(SeqLength):
        for thisK in range(1,max_k+1):
            if SeqLength-i+1 > thisK:
                this_kmer = seq[i:(i+thisK)]
                Kmers[thisK-1].append(this_kmer)
    sumSNa = 0
    sumSNe = 0
    for thisKsub1 in range(max_k):  # https://rosettacode.org/wiki/Entropy
        p = Counter(Kmers[thisKsub1])
        lns = float(len(Kmers[thisKsub1]))
        n = 4**(thisKsub1+1)
        thisSNa = -sum( count * math.log(count/lns, 2) for count in p.values())
        thisSNe = -sum( count * math.log(count/lns, n) for count in p.values())
        sumSNa = sumSNa + thisSNa
        sumSNe = sumSNe + thisSNe
    return [sumSNa,sumSNe]
#ttt = entropy('GATTACATATC',11)
#print(str(ttt))


reTri = re.compile(r"(\w)\1{2,}")

for i in itertools.product('ACTG', repeat=myKmerLen):
    if len(set(i)) < 4:
        continue
    oneKmer = ''.join(i)
    reTriFound = reTri.search(oneKmer)
    if reTriFound:
        continue
    #print(oneKmer)
    #print(reTriFound)
    thisE = entropy(oneKmer,myEntropyKmer)
    if thisE[1] < myMinEntropy:
        continue
    hp=primer3.calcHairpin(oneKmer)
    hd=primer3.calcHomodimer(oneKmer)
    primer3Flag = ['pF','dF']
    if hp.structure_found:
        primer3Flag[0] = 'pT'
    if hd.structure_found:
        primer3Flag[1] = 'dT'
    print(' '.join([oneKmer,str(thisE[0]),str(thisE[1]),primer3Flag[0],primer3Flag[1],str(round(hp.tm,2)),str(round(hd.tm,2))]))

exit()

hp=primer3.calcHairpin('CCCCCATCCGATCAGGGGG')
hd=primer3.calcHomodimer('CCCCCATCCGATCAGGGGG')

"""
grep OK out.kmers.f1 | awk -v size=2.5 '{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }'

grep OK out.kmers.f1 | awk -v size=0.5 '{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }'
"""
