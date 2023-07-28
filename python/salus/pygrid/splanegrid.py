#!/usr/bin/env python3

import sys
import os
import io
import argparse
import pathlib
import gzip
import graphblas as gb
import dinopy
import speedict

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

def db_options():
    opt = speedict.Options(raw_mode=False)
    # create table
    opt.create_if_missing(True)
    # config to more jobs
    opt.set_max_background_jobs(os.cpu_count())
    # configure mem-table to a large value (256 MB)
    opt.set_write_buffer_size(0x10000000)
    opt.set_level_zero_file_num_compaction_trigger(4)
    # configure l0 and l1 size, let them have the same size (1 GB)
    opt.set_max_bytes_for_level_base(0x40000000)
    # 256 MB file size
    opt.set_target_file_size_base(0x10000000)
    # use a smaller compaction multiplier
    opt.set_max_bytes_for_level_multiplier(4.0)
    # use 8-byte prefix (2 ^ 64 is far enough for transaction counts)
    opt.set_prefix_extractor(speedict.SliceTransform.create_max_len_prefix(8))
    # set to plain-table
    opt.set_plain_table_factory(speedict.PlainTableFactoryOptions())
    # by Galaxy
    opt.set_compaction_style(speedict.DBCompactionStyle.level())
    opt.optimize_level_style_compaction(0x20000000) # 512 MB
    opt.increase_parallelism(os.cpu_count())
    opt.set_compression_type(speedict.DBCompressionType.snappy())
    return opt

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
def readSpatial(infile, db):
    global maxBarcodeLen
    global SpatialBarcodeRange_xXyY
    with fileOpener(infile) as f:
        for index,line in enumerate(f, start=1):
            [ seq, Xpos, Ypos, *_ ] = line.split()
            seqLen = len(seq)
            if seqLen > maxBarcodeLen:
                maxBarcodeLen = seqLen
            theXpos = int(float(Xpos))
            theYpos = int(float(Ypos))
            if (not SpatialBarcodeRange_xXyY[0]) or (SpatialBarcodeRange_xXyY[0] > theXpos):
                SpatialBarcodeRange_xXyY[0] = theXpos
            if (not SpatialBarcodeRange_xXyY[1]) or (SpatialBarcodeRange_xXyY[1] < theXpos):
                SpatialBarcodeRange_xXyY[1] = theXpos
            if (not SpatialBarcodeRange_xXyY[2]) or (SpatialBarcodeRange_xXyY[2] > theYpos):
                SpatialBarcodeRange_xXyY[2] = theYpos
            if (not SpatialBarcodeRange_xXyY[3]) or (SpatialBarcodeRange_xXyY[3] < theYpos):
                SpatialBarcodeRange_xXyY[3] = theYpos
            intSeq = dinopy.conversion.encode_twobit(seq)
            #strSeq = dinopy.conversion.decode_twobit(intSeq, maxBarcodeLen, str)
            #pp.pprint([seq, Xpos, Ypos, f'{intSeq:b}', strSeq])
            db[intSeq] = [ theXpos, theYpos, 0 ]
    return index

def updateBarcodesID(infile, db):
    missingCnt = 0
    with fileOpener(infile) as f:
        for index,line in enumerate(f, start=1):
            seq = line.strip()
            intSeq = dinopy.conversion.encode_twobit(seq)
            if db.key_may_exist(intSeq):
                thisValue = db[intSeq]
                thisValue[2] = index
                db[intSeq] = thisValue
            else:
                ++missingCnt
    return missingCnt

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
    OutFileDict['Rdict'] = args.outpath.joinpath('_rdict').as_posix()
    #pp.pprint(OutFileDict)
    args.outpath.mkdir(parents=True, exist_ok=True)
    eprint('[!]Output Files:[',', '.join([ OutFileDict[x].as_posix() for x in spNameTuple]),'].',sep='')
    if args.dryrun: exit(0);

    eprint('[!]Reading spatial file ...', end='')
    spatialDB = speedict.Rdict(OutFileDict['Rdict'],db_options())
    lineCnt = readSpatial(InFileDict['spatial'], spatialDB)
    eprint('\b\b\b\b. Finished with [',lineCnt,'] records. X∈[',','.join(map(str,SpatialBarcodeRange_xXyY[0:2])),'], Y∈[',','.join(map(str,SpatialBarcodeRange_xXyY[2:4])),'].',sep='')
    #pp.pprint(SpatialBarcodeRange_xXyY)
    eprint('[!]Reading barcodes file ...', end='')
    missingCnt = updateBarcodesID(InFileDict['barcodes'], spatialDB)
    eprint('\b\b\b\b. Finished with [',missingCnt,'] missing barcodes.')
    spatialDB.close()
    exit(0);
    spatialDB.destroy(OutFileDict['Rdict'])
    exit(0);
    #outMtx = ''.join((outPrefix,'.mtx'))
    #matrixData = gb.io.mmread(matrixFile)

if __name__ == "__main__":
    gb.init("suitesparse", blocking=True)
    main()  # time ./splanegrid.py ...
