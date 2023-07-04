#!/usr/bin/env python3

import sys
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

substitution_matrix = {'A': {'A':  1, 'C': -2, 'G': -2, 'T': -2, 'U': -2, 'N': -1},
                       'C': {'A': -2, 'C':  1, 'G': -2, 'T': -2, 'U': -2, 'N': -1},
                       'G': {'A': -2, 'C': -2, 'G':  1, 'T': -2, 'U': -2, 'N': -1},
                       'T': {'A': -2, 'C': -2, 'G': -2, 'T':  1, 'U': -2, 'N': -1},
                       'U': {'A': -2, 'C': -2, 'G': -2, 'T': -2, 'U':  1, 'N': -1},
                       'N': {'A': -1, 'C': -1, 'G': -1, 'T': -1, 'U': -1, 'N':  0}}

def getAdapters():
    IlluminaP7adapter = 'CAAGCAGAAGACGGCATACGAGAT'
    TruseqRead1primer5E = 'TCTTTCCCTACACGACGCTCTTCCGATCT'
    SpacialBarcode = 'N' * 30
    SeqScopeUMIpHR1B = 'CCCGTTCGCAACATGTCTGGCGTCATA'   # https://doi.org/10.1016/j.cell.2021.05.010
    SeqScopeHDMIRead1 = 'TCTTGTGACTACAGCACCCTCGACTCTCGC'
    pUMI = 'N' * 10
    polyT = 'T' * 30
    #mRNAdocker = 'VN'
    # mRNA
    #qUMI = 'N' * 9
    TruSeqRead2primerRC = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    IlluminaP5adapterRC = 'GTGTAGATCTCGGTGGTCGCCGTATCATT'

    AdapterL = IlluminaP7adapter + TruseqRead1primer5E + SpacialBarcode + SeqScopeUMIpHR1B + SeqScopeHDMIRead1 + pUMI + polyT
    AdapterR = TruSeqRead2primerRC + IlluminaP5adapterRC
    Adapters = {'AdapterL':DNA(AdapterL), 'AdapterR':DNA(AdapterR)}
    return Adapters

def aliOneSeq(Adapters,targetSeq):
    Alignments = {}
    Scores = {}
    for aid in Adapters.keys():
        Adapter = Adapters[aid]
        TargetSeq = DNA(targetSeq)
        rcTargetSeq = TargetSeq.reverse_complement()
        alignment, score, start_end_positions = global_pairwise_align_nucleotide(
            TargetSeq, Adapter, substitution_matrix=substitution_matrix
        )
        RCalignment, RCscore, RCstart_end_positions = global_pairwise_align_nucleotide(
            rcTargetSeq, Adapter, substitution_matrix=substitution_matrix
        )
        strand = '*'
        thisScore = 0
        if RCalignment == None or score > RCscore:
            resAli = alignment
            strand = '+'
            thisScore = score
        else:
            resAli = RCalignment
            strand = '-'
            thisScore = RCscore
        if thisScore > 5:
            Alignments[aid] = (resAli,strand,thisScore)
    return Alignments

def printAli(Alignments,seq,seqid):
    prnStr = ':'.join((seqid,seq))
    print(prnStr)
    for tid in Alignments.keys():
        (resAli,strand,thisScore) = Alignments[tid]
        tmp = []
        resAli.write(tmp,format='stockholm')
        prnStr = ' '.join((
            tid, strand,
            str(thisScore)
        ))
        print(prnStr)
        prnStr = ''.join(tmp[1:3])
        print(prnStr)
    return

def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

def main():
    if len(sys.argv) < 3 :
        print('Usage:',sys.argv[0],'<fq file> <umi length> >alignments.out',file=sys.stderr,flush=True);
        exit(0);
    fqFile = sys.argv[1]
    umiLength = int(sys.argv[2])

    Adapters = getAdapters()

    n = 4
    with open(fqFile) as fqfile:
        lines = []
        for line in fqfile:
            lines.append(line.rstrip())
            if len(lines) == n:
                record = process(lines)
                record['name'] = ' '.join((record['name'][1:],record['sequence'][:umiLength]))
                record['sequence'] = record['sequence'][umiLength:]
                lines = []
                maxAlignments = aliOneSeq(Adapters,record['sequence'])
                printAli(maxAlignments,record['sequence'],record['name'])

if __name__ == "__main__":
    main()  # time ./alignkmer.py deumi31-9.out.rnk2 1000000 2>/dev/null |head

# GGAGTTGCAAAAGGTCTGCGAGAGTCGAGAGTGCTGTAGTCACAAGATATGACGCCAGACATGTTGCGAACGGGTAAAACTACCCTACACT
