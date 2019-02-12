#!/usr/bin/env pypy3
import sys

SamplesList = ('D3B_Crick', 'D3B_Watson', 'Normal_Crick', 'Normal_Watson', 'D3B_WGS', 'Normal_WGS')

#from collections import defaultdict
#DepthCnt = {key:defaultdict(int) for key in SamplesList}
#yDepthCnt = defaultdict(lambda: defaultdict(int))

from collections import Counter
cDepthCnt = {key:Counter() for key in SamplesList}

cDepthStat = {key:[0,0,0,0,0] for key in SamplesList} # x and x^2

def main():
    import math

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
    for k in SamplesList:
        cDepthStat[k][2] = cDepthStat[k][0] / RecordCnt # E(X)
        cDepthStat[k][3] = cDepthStat[k][1] / RecordCnt # E(X^2)
        cDepthStat[k][4] = math.sqrt(cDepthStat[k][3] - cDepthStat[k][2]*cDepthStat[k][2])   # E(X^2)-E(X)^2
    print('#{}\t{}'.format('Depth','\t'.join(SamplesList)))
    #RecordCntLength = len(str(RecordCnt))
    print( '#N={},SD:\t{}'.format(RecordCnt,'\t'.join(str(round(cDepthStat[col][4],1)) for col in SamplesList)) )
    for depth in range(0,MaxDepth+1):
        #print( '{}\t{}'.format(depth,'\t'.join(str(DepthCnt[col][depth]) for col in SamplesList)) )
        #print( '{}\t{}'.format(depth,'\t'.join(str(yDepthCnt[depth][col]) for col in SamplesList)) )
        print( '{}\t{}'.format(depth,'\t'.join(str(cDepthCnt[col][depth]) for col in SamplesList)) )
        #pass
    pass

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
                    theValue = int(row[k])
                    if theValue > MaxDepth:
                        MaxDepth = theValue
                    #DepthCnt[k][theValue] += 1  # PyPy3:30.54 ns, Python3:22.23 ns
                    #yDepthCnt[theValue][k] += 1 # PyPy3:30.47 ns, Python3:21.50 ns
                    cDepthCnt[k][theValue] += 1 # PyPy3:29.82 ns, Python3:30.61 ns
                    cDepthStat[k][0] += theValue
                    cDepthStat[k][1] += theValue * theValue
                #print(MaxDepth,DepthCnt)
        except KeyboardInterrupt:
            print('\n[!]Ctrl+C pressed.',file=sys.stderr,flush=True)
            pass
        print('[!]Lines Read:[{}], MaxDepth is [{}].'.format(RecordCnt,MaxDepth),file=sys.stderr,flush=True)
    return RecordCnt,MaxDepth

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
            print('[!]cDepthCnt:',get_size(cDepthCnt),file=sys.stderr,flush=True)
            print('[!]DepthCnt: ',get_size(DepthCnt),file=sys.stderr,flush=True)   # smaller
            print('[!]yDepthCnt:',get_size(yDepthCnt),file=sys.stderr,flush=True)
        except FileNotFoundError:
            pass
        except NameError:
            pass

'''
[!]DepthCnt:  993487
[!]yDepthCnt: 1953307
[!]cDepthCnt: 994207

time pypy3 -m cProfile -o t2.cprof ./samdepthplot.py t.tsv.gz 1
pyprof2calltree -i t2.cprof -k
'''
