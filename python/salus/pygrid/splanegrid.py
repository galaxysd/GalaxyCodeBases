#!/usr/bin/env python3

import sys
import os
import argparse
import pathlib
import graphblas as gb

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
    parser.add_argument('-i', '--spatial', type=pathlib.Path, default='spatial.txt', metavar='txt', help='For spatial.txt[.gz]')
    parser.add_argument('-b', '--bin', type=int, required = True, help='grid binning pixels')
    parser.add_argument('-o', '--output-prefix', type=pathlib.Path, default='./gridded/sp', dest='outprefix')
    #parser.add_argument('-r', '--scseq-path', type=pathlib.Path, default='.', dest='scSeqPath')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-r', '--scseq-path', type=pathlib.Path, dest='scSeqPath')
    group.add_argument('-s', '--scseq-files', type=pathlib.Path, nargs=3, action='extend', metavar='<f>', dest='scSeqFiles', help='matrix.mtx[.gz] barcodes.tsv[.gz] features.tsv[.gz]')

    parser.add_argument('-n', '--dryrun', '--dry-run', action='store_true', dest='dryrun')
    parser.add_argument('-z', '--split-zones', dest='zones', type=int, choices=[0,4,5], default=0, help='split to 4 or 5 zones, 0=off')
    #parser.add_argument("files", nargs="*")
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
    pp.pprint(args)
    eprint('[!]scSeq:[',args.scSeqPath,'|',args.scSeqFiles,'], Spatial:[',args.spatial,'], GridBin=[',args.bin,'], SplitZone:[',args.zones,']. OutPut:[',args.outprefix,']',sep='');
    if args.dryrun: exit(0);
    if os.path.dirname(args.outprefix): os.makedirs(os.path.dirname(args.outprefix), exist_ok=True)
    if args.scSeqPath == None:
        args.scSeqFiles.append( args.scSeqFiles[2].with_stem('genes') )
        scSeqFiles = tuple( args.scSeqFiles )
    else:
        scSeqFiles = tuple( args.scSeqPath.joinpath(x) for x in ['matrix.mtx', 'barcodes.tsv', 'features.tsv', 'genes.tsv'] )
    FileDotExts = ('', '.gz')
    #pp.pprint(scSeqFiles)
    FileMatrix = checkFile([scSeqFiles[0]], FileDotExts)
    FileBarcodes = checkFile([scSeqFiles[1]], FileDotExts)
    FileFeatures = checkFile([scSeqFiles[2],scSeqFiles[3]], FileDotExts)
    FileSpatial = checkFile([args.spatial], FileDotExts)
    pp.pprint([FileMatrix, FileBarcodes, FileFeatures, FileSpatial])
    if FileSpatial==None or FileMatrix==None or FileBarcodes==None:
        exit(1)
    exit(0);
    outMtx = ''.join((outPrefix,'.mtx'))
    matrixData = gb.io.mmread(matrixFile)

if __name__ == "__main__":
    gb.init("suitesparse", blocking=True)
    main()  # time ./splanegrid.py ...
