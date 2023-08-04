#!/usr/bin/env python3

import sys
import io
import argparse
import pathlib
import gzip
import pyfastx
import pafpy

import pprint
pp = pprint.PrettyPrinter(indent=4)

def eprint(*args, **kwargs) -> None:
    print(*args, **kwargs, file=sys.stderr, flush=True)

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='intersection of FastQ Read 2 with spatial barcodes, and dump matching Read 1 with spatial coordinates',
        epilog='Contact: <huxs@salus-bio.com>')
    #parser.add_argument('-s', '--spatial', type=pathlib.Path, default='spatial.txt', metavar='file', help='For spatial.txt[.gz]')
    parser.add_argument('-1', '--read1', type=pathlib.Path, default='Unmapped.out.mate1', metavar='file', help='For Unmapped.out.mate1[.gz]')
    #parser.add_argument('-2', '--read2', type=pathlib.Path, default='Unmapped.out.mate2', metavar='file', help='For Unmapped.out.mate2[.gz]')
    parser.add_argument('-p', '--read2-paf', type=pathlib.Path, default='Unmapped.mate2.paf', metavar='file', help='For Unmapped.mate2.paf')
    #parser.add_argument('-m', '--max-mismatch', dest='mismatch', type=int, default=1, help='max allowed mismatch, default=1')
    parser.add_argument('-o', '--output', type=pathlib.Path, default='Unmapped.fq.gz', dest='outfile')
    #parser.add_argument('-z', '--gzip', action=argparse.BooleanOptionalAction, default=True, help='Output gzipped files, default on', dest='gzip')
    #parser.add_argument('-n', '--dryrun', '--dry-run', action='store_true', dest='dryrun')
    parser.add_argument(
        "-v", "--version", action="version",
        version=f"{parser.prog} version 1.0.0"
    )
    return parser

def fileOpener(filename):
    f = open(filename,'rb')
    fh = f
    if (f.read(2) == b'\x1f\x8b'):
        f.seek(0)
        fh = gzip.GzipFile(fileobj=f, mode='rb')
    else:
        f.seek(0)
    #fht = io.TextIOWrapper(fh, encoding='utf-8', line_buffering=True)
    return fh

def main() -> None:
    parser = init_argparse()
    if len(sys.argv) == 1:
        parser.print_help()
        eprint('''
Requirements:
    (01) perl -lane 'print ">",join("_",@F),"\\n$F[0]"' spatial.txt | minimap2 -k 15 -d spatial.miniref - 2>spatial.miniref.log
    (2a) seqtk trimfq -L 30 Unmapped.out.mate2 | minimap2 -x sr spatial.miniref - -k15 -w10 -N1 -t8 -QL2c --eqx --cs --sr --end-bonus 200 --for-only -A4 -B0 -o Unmapped.mate2.paf 2>Unmapped.mate2.paf.log
    (2b) seqtk trimfq -L 30 Unmapped.out.mate2 | minimap2 -x sr spatial.miniref - -k15 -w10 -N1 -t8 -QL2c --eqx --cs --sr --end-bonus 200 --for-only -A4 -B0 2>Unmapped.mate2.paf.log | gzip -1c >Unmapped.mate2.paf.gz
''')
        exit(0);
    args = parser.parse_args()
    pp.pprint(args)
    eprint('[!]Read1:[',args.read1,'], Read2.PAF:[',args.read2_paf,']. OutFile:[',args.outfile,']',sep='');
    skipped = 0
    found = 0
    with gzip.open(args.outfile, mode='wt', compresslevel=1) as fh:
        with fileOpener(args.read2_paf) as fp2:
            with pafpy.PafFile(fp2) as paf:
                for record in paf:
                    if record.is_primary():
                    #pp.pprint(record)
                        for name,seq,qual in pyfastx.Fastq(args.read1.as_posix(), build_index=False):
                            #print('|'.join((name, seq, qual)))
                            if record.qname == name:
                                print( '|'.join((
                                    str(found), str(skipped),
                                    str(record.tags['cg']),str(record.tags['cs']),
                                    name
                                )) )
                                found +=1
                                break
                            else:
                                skipped +=1

if __name__ == "__main__":
    main()  # time ./spffq.py -1 n4457360.Unmapped.out.mate1.gz -p n175410.Unmapped.mate2.paf.gz

'''
[1]+  Running                 perl -lane 'print ">",join("_",@F),"\n$F[0]"' spatial.txt | minimap2 -k 15 -d spatial.miniref - 2> spatial.miniref.log &
[2]+  Running                 seqtk trimfq -L 30 Unmapped.out.mate2 | minimap2 -x sr spatial.miniref - -k15 -w10 -N1 -t8 -QL2c --eqx --cs --sr --end-bonus 200 --for-only -A4 -B0 -o Unmapped2.paf 2> Unmapped2.err &
'''
