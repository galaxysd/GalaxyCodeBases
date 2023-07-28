#!/usr/bin/env python3

import sys
import os
import argparse
import pathlib
import gzip
import graphblas as gb
import dinopy
from speedict import Rdict

import pprint
pp = pprint.PrettyPrinter(indent=4)
# import gc
# gc.collect()

def eprint(*args, **kwargs) -> None:
    print(*args, **kwargs, file=sys.stderr, flush=True)

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='merge scSeq data with spBarcode coordinates that gridded by given binning size ',
        epilog='Contact: <huxs@salus-bio.com>')
    parser.add_argument('-b', '--bin', type=int, required = True, help='grid binning pixels')
    parser.add_argument('-i', '--spatial', type=pathlib.Path, default='spatial.txt', metavar='txt', help='For spatial.txt[.gz]')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-r', '--scseq-path', type=pathlib.Path, dest='scSeqPath')
    group.add_argument('-f', '--scseq-files', type=pathlib.Path, nargs=3, action='extend', metavar='<f>', dest='scSeqFiles', help='matrix.mtx[.gz] barcodes.tsv[.gz] features.tsv[.gz]')

    parser.add_argument('-s', '--split-zones', dest='zones', type=int, choices=[0,4,5], default=0, help='split to 4 or 5 zones, default 0=off')
    #parser.add_argument("files", nargs="*")
    parser.add_argument('-o', '--output-path', type=pathlib.Path, default='./gridded/', dest='outpath')
    parser.add_argument('-z', '--gzip', action=argparse.BooleanOptionalAction, default=True, help='Output gzipped files, default on', dest='gzip')
    parser.add_argument('-n', '--dryrun', '--dry-run', action='store_true', dest='dryrun')
    parser.add_argument(
        "-v", "--version", action="version",
        version=f"{parser.prog} version 1.0.0"
    )
    return parser

def checkFile(PathList, suffixStrs):
    for onePath in PathList:
        for oneExt in suffixStrs:
            thisPath = pathlib.Path(''.join((onePath.as_posix(),oneExt)))
            #print(thisPath)
            if thisPath.exists():
                return thisPath;
    return None;

def main() -> None:
    parser = init_argparse()
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0);
    args = parser.parse_args()
    #pp.pprint(args)
    eprint('[!]GridBin=[',args.bin,'], SplitZone:[',args.zones,']. OutPath:[',args.outpath,']',sep='');
    scFileNameTuple = ('matrix.mtx', 'barcodes.tsv', 'features.tsv', 'genes.tsv')
    spFileNameList = ['spatial.txt']; spFileNameList.extend(scFileNameTuple[0:3])
    #pp.pprint(spFileNameList)
    if args.scSeqPath == None:
        #args.scSeqFiles.append( args.scSeqFiles[2].with_stem('genes') )
        scSeqFiles = tuple( args.scSeqFiles )
    else:
        scSeqFiles = tuple( args.scSeqPath.joinpath(x) for x in scFileNameTuple )
    FileDotExts = ('', '.gz')
    #pp.pprint(scSeqFiles)
    spNameTuple = ('spatial', 'matrix', 'barcodes', 'features')
    spStandardNameDict = dict(zip(spNameTuple,[ '.'.join((fn,'gz')) if args.gzip else fn for fn in spFileNameList ]))
    #pp.pprint(spStandardNameDict)
    InFileDict={}
    InFileDict['spatial'] = checkFile([args.spatial], FileDotExts)
    InFileDict['matrix'] = checkFile([scSeqFiles[0]], FileDotExts)
    InFileDict['barcodes'] = checkFile([scSeqFiles[1]], FileDotExts)
    InFileDict['features'] = checkFile(scSeqFiles[2:], FileDotExts)
    #pp.pprint(inFiles)
    eprint('[!]Confirmed Input Files:[',', '.join([ str(x) if x else '<Missing>' for x in InFileDict.values() ]),'].',sep='')
    for fname in spNameTuple:
        if InFileDict[fname]==None:
            eprint('[x]The',fname,'file is missing !\n')
            exit(1)
    OutFileDict={}
    for fname in spNameTuple:
        OutFileDict[fname] = args.outpath.joinpath(spStandardNameDict[fname])
    #pp.pprint(OutFileDict)
    args.outpath.mkdir(parents=True, exist_ok=True)
    eprint('[!]Output Files:[',', '.join([ x.as_posix() for x in OutFileDict.values()]),'].',sep='')
    if args.dryrun: exit(0);
    exit(0);
    #outMtx = ''.join((outPrefix,'.mtx'))
    #matrixData = gb.io.mmread(matrixFile)

if __name__ == "__main__":
    gb.init("suitesparse", blocking=True)
    main()  # time ./splanegrid.py ...
