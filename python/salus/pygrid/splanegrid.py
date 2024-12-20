#!/usr/bin/env python3
# pip3 install python-graphblas speedict dinopy fast-matrix-market tqdm

#from numba import jit
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
import dinopy
import speedict
import fast_matrix_market
import tqdm
#from collections import defaultdict
import time

import pprint
pp = pprint.PrettyPrinter(indent=4)
# import gc
# gc.collect()

spatialDB = None
mgBoolMtx = None

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

def cmpGridID(a, b):
    print("comparing ", a, " and ", b)
    global spatialDB, args
    #pp.pprint(gridRangeCnt)
    Va = spatialDB[a]
    aXgrid = Va[0] // args.bin
    aYgrid = Va[1] // args.bin
    agridID = aXgrid * gridRangeY + aYgrid
    Vb = spatialDB[b]
    bXgrid = Vb[0] // args.bin
    bYgrid = Vb[1] // args.bin
    bgridID = bXgrid * gridRangeY + bYgrid
    return cmp(agridID,bgridID)

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
            db[intSeq] = [ theXpos, theYpos, None, None ]
            if not index % 1000:
                pbar.update(index - pbar.n)
    pbar.update(index - pbar.n)
    return index

def updateBarcodesID(infile, db, binPixels):
    missingCnt = 0
    global gridRangeCnt
    global mgBoolMtx
    #eprint(len(mtxBar2sp))
    (gridRangeX, gridRangeY, gridCnt) = gridRangeCnt
    pbar = tqdm.tqdm(desc='Barcodes', total=BarcodesCnt, ncols=70, mininterval=0.5, maxinterval=10, unit='', unit_scale=True, dynamic_ncols=True)
    RePattern = re.compile("[-_|,./\\:;`\'!~!@#$%^&*()+= \t\r\n]+")
    with fileOpener(infile) as f:
        for index,line in enumerate(f, start=0):
            [ seq, *_ ] = RePattern.split(line)
            #seq = line.strip()
            intSeq = dinopy.conversion.encode_twobit(seq)
            if db.key_may_exist(intSeq):
            #if intSeq in db:
                thisValue = db[intSeq]
                Xgrid = thisValue[0] // binPixels
                Ygrid = thisValue[1] // binPixels
                gridID = Xgrid * gridRangeY + Ygrid
                #thisValue[2] = index
                #thisValue[3] = gridID
                #db[intSeq] = thisValue
                #mtxBar2sp[index] = intSeq
                #eprint("Pos:",str(index),', ',str(gridID),'.')
                mgBoolMtx[index,gridID] << True
            else:
                ++missingCnt
            pbar.update(index - pbar.n + 1)
    return missingCnt

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

def mkcopy(fromFile, toFile):
    if toFile.exists():
        if toFile.samefile(fromFile):
            return 0
        else:
            return 1
    else:
        try:
            toFile.hardlink_to(fromFile)
            return 0
        except OSError as error :
            eprint(error)
            try:
                toFile.symlink_to(fromFile)
                return 0
            except OSError as error :
                eprint(error)
                return 1

def mkGridSpatial(spFile, scBarcodeFile, gridRangeCnt):
    spFh = gzip.open(spFile, mode='wt', compresslevel=1)
    scFh = gzip.open(scBarcodeFile, mode='wt', compresslevel=1)
    numLen = len(str(gridRangeCnt[2]))
    for i in range(gridRangeCnt[2]):
        barcodeStr = "Barcode{:0{}d}".format(i,numLen)
        ### gridID = Xgrid * gridRangeY + Ygrid
        Ygrid = i // gridRangeCnt[1]
        Xgrid = i - (Ygrid * gridRangeCnt[1])
        #eprint(barcodeStr,str(Xgrid),str(Ygrid))
        print(barcodeStr, file=scFh)
        print(barcodeStr,str(Xgrid),str(Ygrid), file=spFh)
    spFh.close()
    scFh.close()

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
    OutFileDict['mgBoolMtx'] = args.outpath.joinpath('mgBoolMtx.mtx.gz').as_posix()
    #pp.pprint(OutFileDict)
    args.outpath.mkdir(parents=True, exist_ok=True)
    eprint('[!]Output Files:[',', '.join([ OutFileDict[x].as_posix() for x in spNameTuple]),'].',sep='')
    checkmtx(InFileDict['matrix'])
    eprint('[!]Matrix Size: Gene count(nrows)=',GenesCnt,', Barcode count(ncols)=',BarcodesCnt,', Values(nnz)=',mtxNNZ,'.',sep='')
    if args.dryrun: exit(0);
    global spatialDB, SpatialBarcodeRange_xXyY, gridRangeCnt, mgBoolMtx
    #mtxBar2sp = [None] * BarcodesCnt
    mkcopy(InFileDict['features'], OutFileDict['features'])

    start = time.perf_counter()
    eprint('[!]Reading spatial file ...')
    spatialDB = speedict.Rdict(OutFileDict['Rdict'],db_options())
    #spatialDB = {}
    lineCnt = readSpatial(InFileDict['spatial'], spatialDB)
    eprint('[!]Finished with [',lineCnt,'] records. X∈[',','.join(map(str,SpatialBarcodeRange_xXyY[0:2])),'], Y∈[',','.join(map(str,SpatialBarcodeRange_xXyY[2:4])),'].',sep='') # X∈[8000,38000], Y∈[9000,39000]
    #pp.pprint(SpatialBarcodeRange_xXyY)
    SpatialGridRange_xXyY = [ (x // args.bin) for x in SpatialBarcodeRange_xXyY ]
    #gridRangeX = 1 + SpatialGridRange_xXyY[1] - SpatialGridRange_xXyY[0]
    #gridRangeY = 1 + SpatialGridRange_xXyY[3] - SpatialGridRange_xXyY[2]
    (gridRangeX, gridRangeY) = (1+SpatialGridRange_xXyY[1], 1+SpatialGridRange_xXyY[3])
    gridRangeCnt = (gridRangeX, gridRangeY, gridRangeX * gridRangeY)
    eprint('[!]Gridded by Bin [',args.bin,'], GridSize=','×'.join(map(str,(gridRangeX,gridRangeY))),'=',str(gridRangeCnt[2]),'.',sep='' )
    mgBoolMtx = gb.Matrix(gb.dtypes.BOOL, BarcodesCnt, gridRangeCnt[2])
    end1p = time.perf_counter()
    eprint("\tElapsed {}s".format((end1p - start)))

    executor = concurrent.futures.ProcessPoolExecutor(max_workers=4)
    executor.submit( mkGridSpatial, OutFileDict['spatial'],OutFileDict['barcodes'], gridRangeCnt )

    eprint('[!]Reading barcodes file ...')
    #cmpGridID(1,2)
    missingCnt = updateBarcodesID(InFileDict['barcodes'], spatialDB, args.bin)
    eprint('[!]Finished with [',missingCnt,'] missing barcodes.',sep='')
    fh = write2gzip(OutFileDict['mgBoolMtx'])
    gb.io.mmwrite(target=fh, matrix=mgBoolMtx)
    fh.close()
    spatialDB.close()
    end2p = time.perf_counter()
    eprint("\tElapsed {}s".format((end2p - end1p)))

    eprint('[!]Reading Matrix file ...')
    scMtx = gb.io.mmread(source=InFileDict['matrix'], engine='fmm')
    outGrid = scMtx.mxm(mgBoolMtx)  # lazy
    end3p = time.perf_counter()
    eprint("\tElapsed {}s".format((end3p - end2p)))

    eprint('[!]Calculating Matrix file ...')
    outGridResult = outGrid.new()
    end4p = time.perf_counter()
    eprint("\tElapsed {}s".format((end4p - end3p)))

    eprint('[!]Writing Matrix file ...')
    fh = write2gzip(OutFileDict['matrix'])
    gb.io.mmwrite(target=fh, matrix=outGridResult)
    fh.close()
    end5p = time.perf_counter()
    eprint("\tElapsed {}s".format((end5p - end4p)))

    executor.shutdown(wait=True)
    eprint('[!]All done !')
    exit(0);
    #spatialDB.destroy(OutFileDict['Rdict'])    # It is better to keep db file to enable supporting restore running.
    exit(0);
    #outMtx = ''.join((outPrefix,'.mtx'))
    #matrixData = gb.io.mmread(matrixFile)

if __name__ == "__main__":
    gb.init("suitesparse", blocking=False)
    main()  # time ./splanegrid.py -b20 -f matrix2.mtx.gz barcodes.tsv.gz features.tsv.gz -i spatial.txt.gz

# ./splanegrid.py -b20 -i GSE166635_RAW/GSM5076750_HCC2.barcodes.spatial.txt -f GSE166635_RAW/GSM5076750_HCC2.matrix.mtx.gz GSE166635_RAW/GSM5076750_HCC2.barcodes.tsv.gz GSE166635_RAW/GSM5076750_HCC2.features.tsv.gz
