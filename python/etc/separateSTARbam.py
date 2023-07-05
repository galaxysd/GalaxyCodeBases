#!/usr/bin/env python3

import sys
import os
import pysam

def main():
    if len(sys.argv) < 3 :
        print('Usage:',sys.argv[0],'<bam file> <output prefix>',file=sys.stderr,flush=True);
        exit(0);
    bamFile = sys.argv[1]
    outPrefix = sys.argv[2]
    if os.path.dirname(outPrefix): os.makedirs(os.path.dirname(outPrefix), exist_ok=True)
    uniqBamName = ''.join((outPrefix,'.uniq.bam'))
    multiBamName = ''.join((outPrefix,'.multi.bam'))
    nfoName = ''.join((outPrefix,'.info'))
    readCnt = {'Uniq':0, 'Multi':0, 'Secondary':0, 'Total':0}

    with pysam.AlignmentFile(bamFile, "rb") as samIn:
        uniqBam = pysam.AlignmentFile(uniqBamName, "w", template=samIn)
        multiBam = pysam.AlignmentFile(multiBamName, "w", template=samIn)
        for aln in samIn:
            readCnt['Total'] += 1
            if aln.mapq == 255:
                readCnt['Uniq'] += 1
                uniqBam.write(aln)
            elif not aln.is_secondary:
                readCnt['Multi'] += 1
                multiBam.write(aln)
            else:
                readCnt['Secondary'] += 1
        uniqBam.close()
        multiBam.close()

    with open(nfoName,'w') as nfo:
        nfo.write(readCnt)

if __name__ == "__main__":
    main()  # ./separateSTARbam.py Aligned.sortedByCoord.out.bam Separated &
