#!/usr/bin/env pypy3
import sys

SamplesList = ('D3B-Crick', 'D3B-Watson', 'Normal-Crick', 'Normal-Watson', 'D3B-WGS', 'Normal-WGS')
from collections import defaultdict
DepthCnt = dict.fromkeys(SamplesList, defaultdict(int))

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
    inStat(inDepthFile,verbose)

def inStat(inDepthFile,verbose):
    import gzip
    import csv
    RecordCnt = 0
    MaxDepth = 0
    with gzip.open(inDepthFile, 'rt') as tsvin:
        tsvin = csv.DictReader(tsvin, delimiter='\t', fieldnames=('ChrID','Pos')+SamplesList )
        headers = tsvin.fieldnames
        print(headers)
        for row in tsvin:
            #print(', '.join(row[col] for col in SamplesList))
            RecordCnt += 1
            for k in SamplesList:
                if int(row[k]) > MaxDepth:
                    MaxDepth = int(row[k])
                DepthCnt[k][int(row[k])] += 1
            #print(MaxDepth,DepthCnt['Normal-WGS'])


if __name__ == "__main__":
    main()
