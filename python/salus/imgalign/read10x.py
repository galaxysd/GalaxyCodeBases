#!/usr/bin/env python3

import concurrent.futures
import sys
import os
import io
import functools
import re
import argparse
import pathlib
import gzip
import graphblas as gb
import fast_matrix_market
import tqdm
import time
#import recordclass

import dahuffman
# rANS4x16 as used in CRAM should be the best option. See <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8896640/> and <https://github.com/samtools/htscodecs> or *<https://github.com/jkbonfield/rans_static>.
# Others like Arithmetic Coding <https://github.com/fab-jul/torchac>, <https://marknelson.us/posts/2014/10/19/data-compression-with-arithmetic-coding.html>. Or <https://github.com/ahmedfgad/ArithmeticEncodingPython>. And <https://michaeldipperstein.github.io/arithmetic.html#download>.
# Range coding <https://en.wikipedia.org/wiki/Range_coding>. *<https://github.com/richgel999/sserangecoding> <https://github.com/powturbo/Turbo-Range-Coder>.

spatialDB = None
mgBoolMtx = None
#spPosition = recordclass.make_dataclass("Point", [("x",int), ("y",int)], readonly=True)
baseCodec = dahuffman.HuffmanCodec.from_frequencies({'A':812342385,'T':764495360,'C':778170645,'G':622248280,'N':100})
coordinatesCodec = dahuffman.HuffmanCodec.from_frequencies({
    ' ':49041072,'.':98082144,'0':50072308,'1':87291011,'2':85393505,'3':84616786,
    '4':65118253,'5':58134804,'6':58360232,'7':58371094,'8':59921116,'9':59890253,
})
'''
>>> baseCodec.print_code_table()
Bits Code Value Symbol
   2 00       0 'C'
   4 0100     4 _EOF
   4 0101     5 'N'
   3 011      3 'G'
   2 10       2 'T'
   2 11       3 'A'
>>> coordinatesCodec.print_code_table()
Bits Code  Value Symbol
   3 000       0 '3'
   3 001       1 '2'
   3 010       2 '1'
   3 011       3 '.'
   5 10000    16 _EOF
   5 10001    17 ' '
   4 1001      9 '0'
   4 1010     10 '5'
   4 1011     11 '6'
   4 1100     12 '7'
   4 1101     13 '9'
   4 1110     14 '8'
   4 1111     15 '4'
'''

def eprint(*args, **kwargs) -> None:
    print(*args, **kwargs, file=sys.stderr, flush=True)

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='merge scSeq data with spBarcode coordinates that gridded by given binning size ',
        epilog='Contact: <huxs@salus-bio.com>')
    #parser.add_argument('-b', '--bin', type=int, required = True, help='grid binning pixels')
    parser.add_argument('-s', '--spatial', type=pathlib.Path, default='spatial.txt', metavar='txt', dest='scSpatial', help='For spatial.txt[.gz]')
    parser.add_argument('-f', '--scseq-files', type=pathlib.Path, nargs=3, action='extend', metavar='<f>', dest='scSeqFiles', help='matrix.mtx[.gz] barcodes.tsv[.gz] (features|genes).tsv[.gz]')
    #parser.add_argument("files", nargs="*")
    parser.add_argument('-o', '--output-path', type=pathlib.Path, default='./gridded/', dest='outpath')
    parser.add_argument('-z', '--gzip', action=argparse.BooleanOptionalAction, default=True, help='Output gzipped files, default on', dest='gzip')
    parser.add_argument('-n', '--dryrun', '--dry-run', action='store_true', dest='dryrun')
    parser.add_argument('-d', '--debug', action='store_true', dest='dodebug')
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

maxBarcodeLen = 0
SpatialBarcodeRange_xXyY = [0,0,0,0]
gridRangeCnt = ()
def readSpatial(infile, db):
    global maxBarcodeLen
    global SpatialBarcodeRange_xXyY
    global GenesCnt
    global BarcodesCnt
    pbar = tqdm.tqdm(desc='Spatial', total=BarcodesCnt, ncols=70, mininterval=0.5, maxinterval=10, unit='', unit_scale=True, dynamic_ncols=True)
    with fileOpener(infile) as f:
        for index,line in enumerate(f, start=1):
            [ seq, Xpos, Ypos, *_ ] = line.split()
            seqLen = len(seq)
            if seqLen > maxBarcodeLen:
                maxBarcodeLen = seqLen
            theXpos = int(float(Xpos))
            theYpos = int(float(Ypos))
            if (SpatialBarcodeRange_xXyY[0]==0) or (SpatialBarcodeRange_xXyY[0] > theXpos):
                SpatialBarcodeRange_xXyY[0] = theXpos
            if (SpatialBarcodeRange_xXyY[1]==0) or (SpatialBarcodeRange_xXyY[1] < theXpos):
                SpatialBarcodeRange_xXyY[1] = theXpos
            if (SpatialBarcodeRange_xXyY[2]==0) or (SpatialBarcodeRange_xXyY[2] > theYpos):
                SpatialBarcodeRange_xXyY[2] = theYpos
            if (SpatialBarcodeRange_xXyY[3]==0) or (SpatialBarcodeRange_xXyY[3] < theYpos):
                SpatialBarcodeRange_xXyY[3] = theYpos
            binSeq = baseCodec.encode(seq)
            #strSeq = baseCodec.decode(binSeq)
            binCoordStr = coordinatesCodec.encode(f'{Xpos} {Ypos}')
            #db[binSeq] = spPosition(theXpos,theYpos)
            db[binSeq] = binCoordStr
            if not index % 1000:
                pbar.update(index - pbar.n)
    pbar.update(index - pbar.n)
    return index

GenesCnt = 0
BarcodesCnt = 0
mtxNNZ = 0
def checkmtx(mtxfile) -> None:
    mheader = fast_matrix_market.read_header(mtxfile)
    global GenesCnt, BarcodesCnt, mtxNNZ
    GenesCnt = mheader.nrows
    BarcodesCnt = mheader.ncols
    mtxNNZ = mheader.nnz

def write2gzip(outfile):
    fh = gzip.open(outfile, mode='wb', compresslevel=1)
    return fh

def main() -> None:
    parser = init_argparse()
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0);
    args = parser.parse_args()
    scSeqFiles = tuple( args.scSeqFiles )
    eprint('[!]Spatial=[',args.scSpatial,'], 10Xmatrix:[',','.join(map(str,scSeqFiles)),']. OutPath:[',args.outpath,']',sep='');
    #scFileNameTuple = ('matrix.mtx', 'barcodes.tsv', 'features.tsv', 'genes.tsv')
    #spFileNameList = ['spatial.txt']; spFileNameList.extend(scFileNameTuple[0:3])
    checkmtx(scSeqFiles[0].as_posix())
    eprint('[!]Matrix Size: Gene count(nrows)=',GenesCnt,', Barcode count(ncols)=',BarcodesCnt,', Values(nnz)=',mtxNNZ,'.',sep='')
    if args.dryrun: exit(0);
    global spatialDB, SpatialBarcodeRange_xXyY, gridRangeCnt, mgBoolMtx

    start = time.perf_counter()
    eprint('[!]Reading spatial file ...')
    spatialDB = {}
    lineCnt = readSpatial(args.scSpatial.as_posix(), spatialDB)

if __name__ == "__main__":
    gb.init("suitesparse", blocking=False)
    main()  # time ./read10x.py -s HY0726A1.spatial.txt.gz -f HY0726A1.matrix.mtx.gz HY0726A1.barcodes.tsv.gz HY0726A1.features.tsv.gz
