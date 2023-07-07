#!/usr/bin/env python3

import sys
import pysam
from skbio import DNA
from skbio.alignment import global_pairwise_align_nucleotide
import alignfq

def main():
    if len(sys.argv) < 2 :
        print('Usage:',sys.argv[0],'<bam file> >alignments.out',file=sys.stderr,flush=True);
        exit(0);
    bamFile = sys.argv[1]
    Adapters = alignfq.getAdapters()

    with pysam.AlignmentFile(bamFile, "rb") as inbam:
        for aln in inbam:
            # aln.query_name
            fqSeq = aln.get_forward_sequence()
            refNfo = ''.join((
                str(aln.flag),'|',
                aln.reference_name,':',
                str(aln.reference_start),'+',
                str(aln.reference_length),
                ':',aln.cigarstring
            ))
            faID = ' '.join((
                aln.query_name,
                refNfo
            ))
            """
            alnSeq = aln.seq
            if aln.is_reverse:
                alnSeq = str(DNA(alnSeq).reverse_complement())
            if fqSeq != alnSeq:
                eprint('-'.join('[x]',fqSeq,alnSeq))
            """
            maxAlignments = alignfq.aliOneSeq(Adapters,fqSeq)
            alignfq.printAli(maxAlignments,fqSeq,faID)

if __name__ == "__main__":
    main()  # time ./alignbam.py t.multi.bam 2>/dev/null >t.multi.bam.adapterali
