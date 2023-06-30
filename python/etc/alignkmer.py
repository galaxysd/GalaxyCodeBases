#!/usr/bin/env python3

import sys

def buildQuery():
    from skbio import DNA
    from skbio.alignment import StripedSmithWaterman
    Adapters = {
        'adapter1': 'CAAGCAGAAGACGGCATACGAGATTCTTTCCCTACACGACGCTCTTCCGATCT',
        'adapter2': 'CCCGTTCGCAACATGTCTGGCGTCATATCTTGTGACTACAGCACCCTCGACTCTCGCAGACTTTCACCAGTCCATGAT',
        'adapter3': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGTAGATCTCGGTGGTCGCCGTATCATT'
    }
    Querys={}
    for qid in Adapters.keys():
        Seq = Adapters[qid]
        rcSeq = str(DNA(Seq).reverse_complement())
        Query = StripedSmithWaterman(Seq)
        rcQuery = StripedSmithWaterman(rcSeq)
        Querys[qid] = (Query,rcQuery)
    return Querys

def aliOneSeq(QuerySSWs,targetSeq):
    Alignments = {}
    Scores = {}
    for qid in QuerySSWs.keys():
        (Query,rcQuery) = QuerySSWs[qid]
        Alignment = Query(targetSeq)
        rcAlignment = rcQuery(targetSeq)
        strand = '*'
        if rcAlignment == None or Alignment.optimal_alignment_score > rcAlignment.optimal_alignment_score:
            resAli = Alignment
            strand = '+'
        else:
            resAli = rcAlignment
            strand = '-'
        Alignments[qid] = (resAli,strand)
        Scores[qid] = resAli.optimal_alignment_score
    maxScores = [key for key, value in Scores.items() if value == max(Scores.values())]
    maxAlignments = {k: Alignments[k] for k in maxScores}
    return maxAlignments

def main():
    if len(sys.argv) < 3 :
        print('Usage:',sys.argv[0],'<kmc_dump.output_file> <min count>',file=sys.stderr,flush=True);
        exit(0);

    QuerySSWs = buildQuery()
    targetSeq = 'ACATAGTCTGGCGTCATTCTTGTGTACA'
    maxAlignments = aliOneSeq(QuerySSWs,targetSeq)

    for tid in maxAlignments.keys():
        ali = maxAlignments[tid]
        prnStr = ' '.join((
            tid, ali[1], ali[0].cigar,
            str(ali[0].optimal_alignment_score),
            str(ali[0].suboptimal_alignment_score)
        ))
        print(prnStr)
        prnStr = ' -> '.join((
            ali[0].aligned_target_sequence,
            ali[0].aligned_query_sequence
        ))
        print(prnStr)
        print()

if __name__ == "__main__":
    main()
