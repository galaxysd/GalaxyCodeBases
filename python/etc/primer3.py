#!/usr/bin/env python3
import argparse
import logging
import primer3  # https://github.com/libnano/primer3-py

# https://brcaexchange.org/variants and click "Show All Public Data", then click "Download" to get `variants.tsv`.
# gzcat ~/Downloads/variants.tsv.gz|head -n30|awk -F'\t' '{if (length($124)+length($125)==2 || NR==1) print $43,$108,$122,$123,$124,$125}'|column -t
# DBID_LOVD     Genomic_Coordinate_hg38  Chr  Pos       Ref  Alt
# BRCA1_001574  chr17:g.43093299:A>G     17   43093299  A    G
# BRCA1_003516  chr17:g.43078507:T>G     17   43078507  T    G
# BRCA1_004379  chr17:g.43103085:A>G     17   43103085  A    G

from os.path import expanduser
InFile: str = r'~/Downloads/variants.tsv.gz'
InFile = expanduser(InFile)
#InFile: str = r'/Users/galaxy/Downloads/variants.tsv.gz'

#InColNames = ['DBID_LOVD','Chr','Pos','Ref','Alt','Genomic_Coordinate_hg38']
InColNames = ['Chr','Pos','Ref','Alt']
#import numpy
#import pandas as pd
#pd.read_table(InFile,compression='gzip',sep='\t')

import gzip
import csv
Total: int = 0
skipped: int = 0
with gzip.open(InFile, 'rt') as tsvin:
    tsvin = csv.DictReader(tsvin, delimiter='\t')
    #headers = tsvin.fieldnames
    #print(headers)
    for row in tsvin:
        #print(', '.join(row[col] for col in InColNames))
        Total += 1
        if len(row['Ref']) > 1 or len(row['Alt']) > 1 :
            #print(', '.join(row[col] for col in ['Chr','Pos','Ref','Alt']))
            skipped += 1
        else :
            print(', '.join(row[col] for col in InColNames))
#print(b'[!] %(skipped)d InDels skipped in %(Total)d items.' % {b'skipped': skipped, b'Total': Total})
