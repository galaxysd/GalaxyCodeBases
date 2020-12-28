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
        print(p)
        print(n)
        print(lns)
        print(thisSNa)
        print(thisSNe)
        sumSNa = sumSNa + thisSNa
        sumSNe = sumSNe + thisSNe
    return [sumSNa,sumSNe]
ttt = entropy('GGGGGCCCCCCCCC',11)
print(ttt)

exit()

import math
import numpy
import scipy.special as spc

def approximate_entropy(bin_data: str, pattern_length=10):
    """
    Note that this description is taken from the NIST documentation [1]
    [1] http://csrc.nist.gov/publications/nistpubs/800-22-rev1a/SP800-22rev1a.pdf
    As with the Serial test of Section 2.11, the focus of this test is the frequency of all possible overlapping
    m-bit patterns across the entire sequence. The purpose of the test is to compare the frequency of overlapping
    blocks of two consecutive/adjacent lengths (m and m+1) against the expected result for a random sequence.
    :param bin_data: a binary string
    :param pattern_length: the length of the pattern (m)
    :return: the P value
    """
    n = len(bin_data)
    # Add first m+1 bits to the end
    # NOTE: documentation says m-1 bits but that doesnt make sense, or work.
    bin_data += bin_data[:pattern_length + 1:]

    # Get max length one patterns for m, m-1, m-2
    max_pattern = ''
    for i in range(pattern_length + 2):
        max_pattern += '1'

    # Keep track of each pattern's frequency (how often it appears)
    vobs_one = numpy.zeros(int(max_pattern[0:pattern_length:], 2) + 1)
    vobs_two = numpy.zeros(int(max_pattern[0:pattern_length + 1:], 2) + 1)

    for i in range(n):
        # Work out what pattern is observed
        vobs_one[int(bin_data[i:i + pattern_length:], 2)] += 1
        vobs_two[int(bin_data[i:i + pattern_length + 1:], 2)] += 1

    # Calculate the test statistics and p values
    vobs = [vobs_one, vobs_two]
    sums = numpy.zeros(2)
    for i in range(2):
        for j in range(len(vobs[i])):
            if vobs[i][j] > 0:
                sums[i] += vobs[i][j] * math.log(vobs[i][j] / n)
    sums /= n
    ape = sums[0] - sums[1]
    chi_squared = 2.0 * n * (math.log(2) - ape)
    p_val = spc.gammaincc(pow(2, pattern_length-1), chi_squared/2.0)
    return p_val

test_str='GATCTA'
bin_str = ''.join(format(ord(i), 'b') for i in test_str)
print(approximate_entropy(bin_str,6))

exit()
#import csv
#import itertools

with open(InFile) as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter=' ')
    skip = next(tsvreader, [None])[0]
    for i in itertools.product('ACTG', repeat=18):
        oneKmer = ''.join(i)
        if oneKmer == skip:
            skip = next(tsvreader, [None])[0]
        else:
            #print(oneKmer)
            hp=primer3.calcHairpin(oneKmer)
            hd=primer3.calcHomodimer(oneKmer)
            if hp.structure_found or hd.structure_found:
                print(' '.join([str(round(hp.tm,2)),str(round(hd.tm,2)),oneKmer]))
            else:
                print(' '.join(['OK',oneKmer]))

exit()

hp=primer3.calcHairpin('CCCCCATCCGATCAGGGGG')
hd=primer3.calcHomodimer('CCCCCATCCGATCAGGGGG')

print(t)

