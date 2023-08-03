#!/usr/bin/env python3

import sys
import argparse
import pathlib
import gzip

import pprint
pp = pprint.PrettyPrinter(indent=4)

def eprint(*args, **kwargs) -> None:
    print(*args, **kwargs, file=sys.stderr, flush=True)

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='intersection of FastQ Read 2 with spatial barcodes, and dump matching Read 1 with spatial coordinates',
        epilog='Contact: <huxs@salus-bio.com>')
    parser.add_argument('-s', '--spatial', type=pathlib.Path, default='spatial.txt', metavar='file', help='For spatial.txt[.gz]')
    parser.add_argument('-1', '--read1', type=pathlib.Path, default='Unmapped.out.mate1', metavar='file', help='For Unmapped.out.mate1[.gz]')
    parser.add_argument('-2', '--read2', type=pathlib.Path, default='Unmapped.out.mate2', metavar='file', help='For Unmapped.out.mate2[.gz]')
    parser.add_argument('-m', '--max-mismatch', dest='mismatch', type=int, default=1, help='max allowed mismatch, default=1')
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
    fht = io.TextIOWrapper(fh, encoding='utf-8', line_buffering=True)
    return fht

def main() -> None:
    parser = init_argparse()
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0);
    args = parser.parse_args()
    pp.pprint(args)
    eprint('[!]maxMisMatch=[',args.mismatch,']. Spatial:[',args.spatial,'], Read1:[',args.read1,'], Read2:[',args.read2,']. OutFile:[',args.outfile,']',sep='');


if __name__ == "__main__":
    main()  # time ./spffq.py
