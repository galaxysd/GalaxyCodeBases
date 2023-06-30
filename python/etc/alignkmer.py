#!/usr/bin/env python3

from skbio.alignment import StripedSmithWaterman
from skbio import DNA

Adapters = {
    'adapter1': 'CAAGCAGAAGACGGCATACGAGATTCTTTCCCTACACGACGCTCTTCCGATCT',
    'adapter2': 'CCCGTTCGCAACATGTCTGGCGTCATATCTTGTGACTACAGCACCCTCGACTCTCGCAGACTTTCACCAGTCCATGAT',
    'adapter3': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGTAGATCTCGGTGGTCGCCGTATCATT'
}

Seq = 'ACATAGTCTGGCGTCATTCTTGTGTACA'
Query = StripedSmithWaterman(Seq)

rcSeq = str(DNA(Seq).reverse_complement())
rcQuery = StripedSmithWaterman(rcSeq)

Alignments = {}
Scores = {}
for tid in Adapters.keys():
    target_sequence = Adapters[tid]
    Alignment = Query(target_sequence)
    rcAlignment = rcQuery(target_sequence)
    strand = '*'
    if rcAlignment == None or Alignment.optimal_alignment_score > rcAlignment.optimal_alignment_score:
        resAli = Alignment
        strand = '+'
    else:
        resAli = rcAlignment
        strand = '-'
    Alignments[tid] = (resAli,strand)
    Scores[tid] = resAli.optimal_alignment_score
maxScores = [key for key, value in Scores.items() if value == max(Scores.values())]

for tid in maxScores:
    ali = Alignments[tid]
    prnStr = ' '.join((
        tid, ali[1], ali[0].cigar,
        str(ali[0].optimal_alignment_score),
        str(ali[0].suboptimal_alignment_score)
    ))
    print(prnStr)
    prnStr = ' -> '.join((
        ali[0].aligned_query_sequence,ali[0].aligned_target_sequence
    ))
    print(prnStr)
    print()
