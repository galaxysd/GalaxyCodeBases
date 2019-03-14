#!/usr/bin/env pypy3
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import concurrent
import multiprocessing
import csv
#import time
#import random

SamplesList = ('D3B_Crick', 'D3B_Watson', 'Normal_Crick', 'Normal_Watson', 'D3B_WGS', 'Normal_WGS')

ChunkSize = 1024 * 128
verbose = 0
Nworkers = 16

def main():
    import math

    if len(sys.argv) < 3 :
        print('Usage:',sys.argv[0],'<samtools.depth.gz> <out.tsv> [verbose=0]',file=sys.stderr,flush=True)
        exit(0)
    try:
        verbose = int(sys.argv[3])
    except: # `except IndexError:` and `except ValueError:`
        verbose = 0

    inDepthFile = sys.argv[1]
    outFile = sys.argv[2]
    print('From:[{}], To:[{}].\nVerbose: [{}].'.format(inDepthFile,outFile,verbose),file=sys.stderr,flush=True)
    RecordCnt,MaxDepth,cDepthCnt,cDepthStat = CallStat(inDepthFile)
    if RecordCnt == 0:
        RecordCnt = 1
        print('[!]Empty Record !',file=sys.stderr,flush=True)
    for k in SamplesList:
        cDepthStat[k][2] = cDepthStat[k][0] / RecordCnt # E(X)
        cDepthStat[k][3] = cDepthStat[k][1] / RecordCnt # E(X^2)
        cDepthStat[k][4] = math.sqrt(cDepthStat[k][3] - cDepthStat[k][2]*cDepthStat[k][2])   # E(X^2)-E(X)^2
    tsvout = open(outFile, 'wt')
    print('#{}\t{}'.format('Depth','\t'.join(SamplesList)),file=tsvout)
    #RecordCntLength = len(str(RecordCnt))
    print( '#N={},SD:\t{}'.format(RecordCnt,'\t'.join(str(round(cDepthStat[col][4],1)) for col in SamplesList)),file=tsvout)
    for depth in range(0,MaxDepth+1):
        print( '{}\t{}'.format(depth,'\t'.join(str(cDepthCnt[col][depth]) for col in SamplesList)),file=tsvout)
        #pass
    #print('#MaxDepth={}'.format(MaxDepth),file=tsvout)
    tsvout.close()
    pass

def CallStat(inDepthFile):
    import gzip
    import itertools
    import functools
    import signal
    RecordCnt = 0
    MaxDepth = 0
    cDepthCnt = {key:Counter() for key in SamplesList}
    cDepthStat = {key:[0,0,0,0,0] for key in SamplesList} # x and x^2
    #lines_queue = Queue()
    oldSignal=signal.signal(signal.SIGINT, signal.SIG_IGN)
    manager = multiprocessing.Manager()
    lines_queue = manager.Queue()
    executor = ProcessPoolExecutor(max_workers=Nworkers)
    fStator = functools.partial(iStator,inQueue=lines_queue,inSamplesList=SamplesList)
    workers_pool = [executor.submit(fStator) for i in range(Nworkers)]
    signal.signal(signal.SIGINT, oldSignal)
    try:
        with gzip.open(inDepthFile, 'rt') as tsvfin:
            while True:
                lines = tsvfin.readlines(ChunkSize)
                lines_queue.put(lines)
                if not lines:
                    for i in range(Nworkers):
                        lines_queue.put(b'\n\n')
                    #stater_pool.close()
                    break
    except KeyboardInterrupt:
        print('\n[!]Ctrl+C pressed.',file=sys.stderr,flush=True)
        for i in range(Nworkers):
            lines_queue.put(b'\n\n')
        stater_pool.close()
        pass
    for future in concurrent.futures.as_completed(workers_pool):
            try:
                (iRecordCnt,iMaxDepth,icDepthCnt,icDepthStat) = future.result()
            except Exception as exc:
                print('Generated an exception: %s' % (exc))
            else:
                RecordCnt += iRecordCnt
                if iMaxDepth > MaxDepth:
                    MaxDepth = iMaxDepth
                for k in SamplesList:
                    cDepthCnt[k].update(icDepthCnt[k])
                    cDepthStat[k][0] += icDepthStat[k][0]
                    cDepthStat[k][1] += icDepthStat[k][1]
    return RecordCnt,MaxDepth,cDepthCnt,cDepthStat

def iStator(inQueue,inSamplesList):
    # Looking up things in global scope takes longer then looking up stuff in local scope. <https://stackoverflow.com/a/54645851/159695>
    cDepthCnt = {key:Counter() for key in inSamplesList}
    cDepthStat = {key:[0,0] for key in inSamplesList} # x and x^2
    RecordCnt = 0
    MaxDepth = 0
    for lines in iter(inQueue.get, b'\n\n'):
        tsvin = csv.DictReader(lines, delimiter='\t', fieldnames=('ChrID','Pos')+inSamplesList )
        for row in tsvin:
            #print(', '.join(row[col] for col in inSamplesList))
            RecordCnt += 1
            for k in inSamplesList:
                theValue = int(row[k])
                if theValue > MaxDepth:
                    MaxDepth = theValue
                cDepthCnt[k][theValue] += 1
                cDepthStat[k][0] += theValue
                cDepthStat[k][1] += theValue * theValue
            #print(MaxDepth,DepthCnt)
        #print('[!]{} Lines Read:[{}], MaxDepth is [{}].'.format(multiprocessing.current_process().name,RecordCnt,MaxDepth),file=sys.stderr,flush=True)
    return RecordCnt,MaxDepth,cDepthCnt,cDepthStat

if __name__ == "__main__":
    #multiprocessing.freeze_support()
    main()  # time python3 ./samdepthstat.pool.py t.tsv.gz out.pool

    import platform
    get_implementation_name = platform.python_implementation()
    if get_implementation_name == 'CPython':
        try:
            # get_size only works under Python, not PyPy.
            mylibfile = "../lib/pysize/pysize.py"
            with open(mylibfile) as f:
                code = compile(f.read(), mylibfile, 'exec')
                exec(code, globals(), locals())
            print('[!]cDepthCnt:',get_size(cDepthCnt),file=sys.stderr,flush=True)
        except FileNotFoundError:
            pass
        except NameError:
            pass

'''
[!]cDepthCnt: 994207

# https://stackoverflow.com/questions/11041683/python-multiprocess-profiling
time pypy3 -m cProfile -o pool.cprof ./samdepthstat.pool.py t.tsv.gz t.out
pyprof2calltree -i pool.cprof -k

(time python3 ./samdepthstat.pool.py t.tsv.gz out.pool) 2>time.pool
'''
