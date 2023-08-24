#!/usr/bin/env python3
###!/share/src/third_party/miniconda/miniconda_py39_23/bin/python3
# pip3 install python-graphblas dinopy fast-matrix-market tqdm

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
import fast_matrix_market
import tqdm
import time
import recordclass

spatialDB = None
mgBoolMtx = None
spPosition = recordclass.make_dataclass("Point", [("x",int), ("y",int)], readonly=True)

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
    parser.add_argument('-d', '--debug', action='store_true', dest='dodebug')
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
            if (SpatialBarcodeRange_xXyY[0]==0) or (SpatialBarcodeRange_xXyY[0] > theXpos):
                SpatialBarcodeRange_xXyY[0] = theXpos
            if (SpatialBarcodeRange_xXyY[1]==0) or (SpatialBarcodeRange_xXyY[1] < theXpos):
                SpatialBarcodeRange_xXyY[1] = theXpos
            if (SpatialBarcodeRange_xXyY[2]==0) or (SpatialBarcodeRange_xXyY[2] > theYpos):
                SpatialBarcodeRange_xXyY[2] = theYpos
            if (SpatialBarcodeRange_xXyY[3]==0) or (SpatialBarcodeRange_xXyY[3] < theYpos):
                SpatialBarcodeRange_xXyY[3] = theYpos
            intSeq = dinopy.conversion.encode_twobit(seq)
            #strSeq = dinopy.conversion.decode_twobit(intSeq, maxBarcodeLen, str)
            db[intSeq] = spPosition(theXpos,theYpos)
            #db[intSeq] = [ theXpos, theYpos, None, None ]
            if not index % 1000:
                pbar.update(index - pbar.n)
    pbar.update(index - pbar.n)
    return index

def updateBarcodesID(infile, db, binPixels):
    missingCnt = 0
    global gridRangeCnt
    global mgBoolMtx
    #eprint(len(mtxBar2sp))
    (gridRangeX, gridRangeY, gridCnt, minX, minY) = gridRangeCnt
    pbar = tqdm.tqdm(desc='Barcodes', total=BarcodesCnt, ncols=70, mininterval=0.5, maxinterval=10, unit='', unit_scale=True, dynamic_ncols=True)
    RePattern = re.compile("[-_|,./\\:;`\'!~!@#$%^&*()+= \t\r\n]+")
    with fileOpener(infile) as f:
        for index,line in enumerate(f, start=0):
            [ seq, *_ ] = RePattern.split(line)
            #seq = line.strip()
            intSeq = dinopy.conversion.encode_twobit(seq)
            if intSeq in db:
                thisValue = db[intSeq]
                Xgrid = (thisValue.x - minX) // binPixels
                Ygrid = (thisValue.y - minY) // binPixels
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
            if sys.version_info < (3, 10):
                fromFile.link_to(toFile)
            else:
                toFile.hardlink_to(fromFile)
            return 0
        except OSError as error :
            #eprint(error)
            try:
                toFile.symlink_to(fromFile.absolute())
                return 0
            except OSError as error :
                eprint(error)
                return 1
        except AttributeError as error :
            toFile.symlink_to(fromFile.absolute())

def mkGridSpatial(spFile, scBarcodeFile, gridRangeCnt, binSize):
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
        print(barcodeStr,str(gridRangeCnt[3]+Xgrid*binSize),str(gridRangeCnt[4]+Ygrid*binSize), file=spFh)
    spFh.close()
    scFh.close()

def main() -> None:
    parser = init_argparse()
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0);
    args = parser.parse_args()
    eprint('[!]GridBin=[',args.bin,'], SplitZone:[',args.zones,']. OutPath:[',args.outpath,']',sep='');
    scFileNameTuple = ('matrix.mtx', 'barcodes.tsv', 'features.tsv', 'genes.tsv')
    spFileNameList = ['spatial.txt']; spFileNameList.extend(scFileNameTuple[0:3])
    if args.scSeqPath == None:
        #args.scSeqFiles.append( args.scSeqFiles[2].with_stem('genes') )
        scSeqFiles = tuple( args.scSeqFiles )
    else:
        scSeqFiles = tuple( args.scSeqPath.joinpath(x) for x in scFileNameTuple )
    FileDotExts = ('', '.gz')
    spNameTuple = ('spatial', 'matrix', 'barcodes', 'features')
    spStandardNameDict = dict(zip(spNameTuple,[ '.'.join((fn,'gz')) if args.gzip else fn for fn in spFileNameList ]))
    InFileDict={}
    InFileDict['spatial'] = checkFile([args.spatial], FileDotExts)
    InFileDict['matrix'] = checkFile([scSeqFiles[0]], FileDotExts)
    InFileDict['barcodes'] = checkFile([scSeqFiles[1]], FileDotExts)
    InFileDict['features'] = checkFile(scSeqFiles[2:], FileDotExts)
    spStandardNameDict['features'] = InFileDict['features'].name
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
    spatialDB = {}
    lineCnt = readSpatial(InFileDict['spatial'], spatialDB)
    eprint('[!]Finished with [',lineCnt,'] records. X∈[',','.join(map(str,SpatialBarcodeRange_xXyY[0:2])),'], Y∈[',','.join(map(str,SpatialBarcodeRange_xXyY[2:4])),'].',sep='') # X∈[8000,38000], Y∈[9000,39000]
    SpatialGridRange_xXyY = [ (x // args.bin) for x in SpatialBarcodeRange_xXyY ]
    #gridRangeX = 1 + SpatialGridRange_xXyY[1] - SpatialGridRange_xXyY[0]
    #gridRangeY = 1 + SpatialGridRange_xXyY[3] - SpatialGridRange_xXyY[2]
    (gridRangeX, gridRangeY) = (1+SpatialGridRange_xXyY[1]-SpatialGridRange_xXyY[0], 1+SpatialGridRange_xXyY[3]-SpatialGridRange_xXyY[2])
    gridRangeCnt = (gridRangeX, gridRangeY, gridRangeX * gridRangeY, SpatialBarcodeRange_xXyY[0], SpatialBarcodeRange_xXyY[2])
    eprint('[!]Gridded by Bin [',args.bin,'], GridSize=','×'.join(map(str,(gridRangeX,gridRangeY))),'=',str(gridRangeCnt[2]),'.',sep='' )
    mgBoolMtx = gb.Matrix(gb.dtypes.BOOL, BarcodesCnt, gridRangeCnt[2])
    end1p = time.perf_counter()
    eprint("\tElapsed {}s".format((end1p - start)))

    executor = concurrent.futures.ProcessPoolExecutor(max_workers=4)
    executor.submit( mkGridSpatial, OutFileDict['spatial'],OutFileDict['barcodes'], gridRangeCnt, args.bin )

    eprint('[!]Reading barcodes file ...')
    #cmpGridID(1,2)
    missingCnt = updateBarcodesID(InFileDict['barcodes'], spatialDB, args.bin)
    eprint('[!]Finished with [',missingCnt,'] missing barcodes.',sep='')
    if args.dodebug:
        fh = write2gzip(OutFileDict['mgBoolMtx'])
        gb.io.mmwrite(target=fh, matrix=mgBoolMtx)
        fh.close()
    #spatialDB.close()
    end2p = time.perf_counter()
    eprint("\tElapsed {}s".format((end2p - end1p)))

    eprint('[!]Reading Matrix file ...')
    scMtx = gb.io.mmread(source=InFileDict['matrix'], engine='fmm')
    #outGrid = gb.Matrix(gb.dtypes.INT64, scMtx.nrows, mgBoolMtx.ncols)
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

if __name__ == "__main__":
    gb.init("suitesparse", blocking=False)
    main()  # time ./splanegridict.py -b20 -f matrix2.mtx.gz barcodes.tsv.gz features.tsv.gz -i spatial.txt.gz

# ./splanegridict.py -b100 -i GSE166635_RAW/GSM5076749_HCC1.barcodes.spatial.txt -f GSE166635_RAW/GSM5076749_HCC1.matrix.mtx.gz GSE166635_RAW/GSM5076749_HCC1.barcodes.tsv.gz GSE166635_RAW/GSM5076749_HCC1.features.tsv.gz -o gridded
