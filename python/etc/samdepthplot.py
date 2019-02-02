#!/usr/bin/env pypy3
import sys

SamplesList = ('D3B_Crick', 'D3B_Watson', 'Normal_Crick', 'Normal_Watson', 'D3B_WGS', 'Normal_WGS')
from collections import defaultdict
DepthCnt = dict.fromkeys(SamplesList, defaultdict(int))
#nDepthCnt = defaultdict(lambda: defaultdict(int))

def main():
    if len(sys.argv) < 3 :
        print('Usage:',sys.argv[0],'<samtools.depth.gz> <outprefix> [verbose=0]',file=sys.stderr,flush=True)
        exit(0)
    try:
        verbose = int(sys.argv[3])
    except: # `except IndexError:` and `except ValueError:`
        verbose = 0

    inDepthFile = sys.argv[1]
    outPrefix = sys.argv[2]
    print('From:[{}], To:[{}.out].\nVerbose: [{}].'.format(inDepthFile,outPrefix,verbose),file=sys.stderr,flush=True)
    RecordCnt,MaxDepth = inStat(inDepthFile,verbose)
    print('{}\t{}'.format('#depth','\t'.join(SamplesList)))
    for depth in range(0,MaxDepth):
        print( '{}\t{}'.format(depth,'\t'.join(str(DepthCnt[col][depth]) for col in SamplesList)) )
    plotDepth(RecordCnt,MaxDepth)

def inStat(inDepthFile,verbose):
    import gzip
    import csv
    RecordCnt = 0
    MaxDepth = 0
    with gzip.open(inDepthFile, 'rt') as tsvin:
        tsvin = csv.DictReader(tsvin, delimiter='\t', fieldnames=('ChrID','Pos')+SamplesList )
        #headers = tsvin.fieldnames
        #print(headers)
        try:
            for row in tsvin:
                #print(', '.join(row[col] for col in SamplesList))
                RecordCnt += 1
                for k in SamplesList:
                    if int(row[k]) > MaxDepth:
                        MaxDepth = int(row[k])
                    DepthCnt[k][int(row[k])] += 1
                    #nDepthCnt[int(row[k])][k] += 1
                #print(MaxDepth,DepthCnt)
        except KeyboardInterrupt:
            print('\n[!]Ctrl+C pressed.',file=sys.stderr,flush=True)
            pass
        print('[!]Lines Read:[{}], MaxDepth is [{}].'.format(RecordCnt,MaxDepth),file=sys.stderr,flush=True)
    return RecordCnt,MaxDepth

def plotDepth(RecordCnt,MaxDepth):
    import matplotlib.pyplot as plt
    return

if __name__ == "__main__":
    main()  # time python3 ./samdepthplot.py t.tsv.gz 1

    import platform
    get_implementation_name = platform.python_implementation()
    if get_implementation_name == 'CPython':
        try:
            # get_size only works under Python, not PyPy.
            mylibfile = "../lib/pysize/pysize.py"
            with open(mylibfile) as f:
                code = compile(f.read(), mylibfile, 'exec')
                exec(code, globals(), locals())
            print('[!]DepthCnt:',get_size(DepthCnt),file=sys.stderr,flush=True)   # smaller
            #print('nDepthCnt:',get_size(nDepthCnt))
        except FileNotFoundError:
            pass
