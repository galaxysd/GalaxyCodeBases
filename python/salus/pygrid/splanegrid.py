#!/usr/bin/env python3

import sys

# import gc
# gc.collect()

def eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr, flush=True)

def main():
    if len(sys.argv) != 3 :
        eprint('Usage:',sys.argv[0],'<matrix.mtx.gz> <bin size>');
        exit(0);
    matrixFile = sys.argv[1]
    binSize = int(sys.argv[2])
    eprint('[!]From:[',matrixFile,'], Bin:[',binSize,'].\n', sep='')

if __name__ == "__main__":
    main()  # time ./splanegrid.py deumi31-9.out.rnk2 1000000 |head
