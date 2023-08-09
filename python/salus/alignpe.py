#!/usr/bin/env python3
###!/scratch/src/mambaforge/envs/py311/bin/python3

import sys
#sys.path.append('/scratch/src/mambaforge/envs/py311/lib/python3.11/site-packages')
import io
import argparse
import pathlib
import gzip
#import pyfastx
import mappy
import tqdm

import pprint
pp = pprint.PrettyPrinter(indent=4)

def eprint(*args, **kwargs) -> None:
    print(*args, **kwargs, file=sys.stderr, flush=True)

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Align Salus PE barcode with 1-100 as forward and 101-150 as reverse.',
        epilog='Contact: <huxs@salus-bio.com>')
    parser.add_argument('outfile', type=argparse.FileType('w', encoding='UTF-8'), help='Alignment Results', nargs='?', default='-')
    parser.add_argument('-i', '--read1', type=pathlib.Path, required=True, dest='infile', help='For Lane01_fastq_R1.fq[.gz]')
    #parser.add_argument('-n', '--dryrun', '--dry-run', action='store_true', dest='dryrun')
    parser.add_argument(
        "-v", "--version", action="version",
        version=f"{parser.prog} version 1.0.0"
    )
    return parser

def main() -> None:
    parser = init_argparse()
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0);
    args = parser.parse_args()
    pp.pprint(args)
    matchNumbers = list('1234567890'*15)
    for name,seq,qual in mappy.fastx_read(args.infile.as_posix()):
        seqF = seq[:100]
        seqR = seq[100:]
        AlnRef = mappy.Aligner(seq=seqF,preset='sr',k=11)
        for hit in AlnRef.map(seqR, cs=False, MD=False): # traverse alignments
            if hit.is_primary:
                if hit.NM or hit.cigar_str != '30M':
                    print( ', '.join(map(str,[hit.ctg, hit.r_st, hit.r_en, hit.strand, hit.q_st, hit.q_en, hit.cigar_str, hit.NM, name])), file=args.outfile )
                    seqFstr = seqF[:30] + seqF[30:].lower()
                    seqRstr = seqR[:30] + seqR[30:].lower()
                    queryST = hit.q_st
                    if hit.strand == -1:
                        revR = mappy.revcomp(seqRstr)
                        queryST = len(seqR) - hit.q_st - hit.q_en
                        queryStr = revR
                    elif hit.strand == 1:
                        queryStr = seqRstr
                    (firstPadding, secondPadding) = (0,0)
                    if hit.r_st > queryST:
                        secondPadding = hit.r_st - queryST
                    elif hit.r_st < queryST:
                        firstPadding = queryST - hit.r_st
                    rStr = ' '*firstPadding + seqFstr
                    qStr = ' '*secondPadding + queryStr
                    pairs = zip( list(rStr),list(qStr) )
                    matchStr = ''
                    matchNumberIndex = 0
                    for item in pairs:
                        if item[0].lower() == item[1].lower():
                            #matchStr += '|'
                            matchStr += matchNumbers[matchNumberIndex]
                            matchNumberIndex += 1
                        else:
                            if item[0] == ' ':
                                matchStr += '_'
                            elif item[1] == ' ':
                                matchStr += '‾'
                            else:
                                matchStr += '●'
                    print( rStr,matchStr,qStr,sep='\n', file=args.outfile)
        args.outfile.flush()
    exit(0);

if __name__ == "__main__":
    main()  # time ./alignpe.py -i pe150.fq.gz
