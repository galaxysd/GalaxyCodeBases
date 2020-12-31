#!/usr/bin/env python3

InFile = 'out.kmers.1k'

import csv
import primer3
import itertools

import math
from collections import Counter

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



with open(InFile) as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter=' ')
    skip = next(tsvreader, [None])[0]
    for i in itertools.product('ACTG', repeat=18):
        oneKmer = ''.join(i)
        if oneKmer == skip:
            skip = next(tsvreader, [None])[0]
        else:
            #print(oneKmer)
            thisE = entropy(oneKmer,15)
            if thisE[1] < 53:
                continue
            hp=primer3.calcHairpin(oneKmer)
            hd=primer3.calcHomodimer(oneKmer)
            if hp.structure_found or hd.structure_found:
                print(' '.join(['ST',str(round(hp.tm,2)),str(round(hd.tm,2)),oneKmer,str(thisE[0]),str(thisE[1])]))
            else:
                print(' '.join(['OK','0','0',oneKmer,str(thisE[0]),str(thisE[1])]))

exit()

hp=primer3.calcHairpin('CCCCCATCCGATCAGGGGG')
hd=primer3.calcHomodimer('CCCCCATCCGATCAGGGGG')

"""
grep OK out.kmers.f1 | awk -v size=2.5 '{ b=int($5/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }'

grep OK out.kmers.f1 | awk -v size=0.5 '{ b=int($6/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }'
"""
