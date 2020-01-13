#!/usr/bin/env python3

import sys
import os
from pathlib import Path
from pysam import VariantFile

import pprint

def readVcfDat(fname):
    vcf_in = VariantFile(fname)
    for rec in vcf_in.fetch():
        pprint.pprint(rec)
        print ([rec.chrom,rec.pos,rec.id,rec.alleles])
        if rec.id.startswith('rs'):
            # SNP
        else:
            # STR
        for key,value in rec.samples.iteritems():
            print((key, value['GT']))
            

def main():
    try:
        p = Path('unphased_all_vcf')
        for f in p.iterdir():
            if f.is_file():
                if '.vcf' in f.suffixes:
                    print((f,f.suffixes))
                    readVcfDat(f)

        print ('[!]Done', file=sys.stderr)
        sys.stdout.flush()
    except BrokenPipeError: # https://stackoverflow.com/a/58517082
        # https://docs.python.org/3/library/signal.html#note-on-sigpipe
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE

if __name__ == '__main__':
    main()
