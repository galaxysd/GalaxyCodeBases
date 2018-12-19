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
InFile: str = r'~/tmp/variants.tsv.0.gz'
InFile = expanduser(InFile)

InRef: str = expanduser(r'~/tmp/GRCh38.fa')
from pyfaidx import Fasta
RefSeqs = Fasta(InRef)
#print(RefSeqs['chr1'])
# Random access of BGZip is not supported now, see https://github.com/mdshw5/pyfaidx/issues/126

#InColNames = ['DBID_LOVD','Chr','Pos','Ref','Alt','Genomic_Coordinate_hg38']
InColNames = ['Chr','Pos','Ref','Alt']
#import numpy
#import pandas as pd
#pd.read_table(InFile,compression='gzip',sep='\t')

import gzip
import csv
Total: int = 0
Skipped: int = 0
from typing import Dict, List, Tuple
InData: Dict[str,Dict[int,Tuple[str,List[str]]]] = {}

'''
第一种引物，上游引物3‘端设一个，下游距离300bp-400bp设置
第二种，目标点上游100bp设置上游引物，不要覆盖目标点，下游，200-300，
只考虑一对引物中间的部分，引物本身不考虑。
Tm 参考范围55-62
'''

thePara: Dict[str,int] = dict(MaxAmpLen=400, MinAmpLen=300, P5Up1=0, P5Up2=100,
    TmMax=63, TmMin=55, TmDeltra=5,
    PrimerLenMin=25, PrimerLenMax=36, Mode2LeftMax=100
    )

with gzip.open(InFile, 'rt') as tsvin:
    tsvin = csv.DictReader(tsvin, delimiter='\t')
    #headers = tsvin.fieldnames
    #print(headers)
    for row in tsvin:
        #print(', '.join(row[col] for col in InColNames))
        Total += 1
        if len(row['Ref']) > 1 or len(row['Alt']) > 1 :
            #print(', '.join(row[col] for col in ['Chr','Pos','Ref','Alt']))
            Skipped += 1
        else :
            print(', '.join(row[col] for col in InColNames))
            row['Pos'] = int(row['Pos'])
            if row['Chr'] in InData :
                if row['Pos'] in InData[row['Chr']] :
                    InData[row['Chr']][row['Pos']][1].append(row['Alt'])
                    #print(InData[row['Chr']][row['Pos']])
                else :
                    InData[row['Chr']][row['Pos']] = (row['Ref'],[row['Alt']])
            else :
                InData[row['Chr']] = { row['Pos'] : (row['Ref'],[row['Alt']]) }

Primer3GlobalArgs: Dict = {
    'PRIMER_OPT_SIZE': 2+thePara['PrimerLenMin'],
    'PRIMER_PICK_INTERNAL_OLIGO': 1,
    'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': thePara['PrimerLenMin'],
    'PRIMER_MAX_SIZE': thePara['PrimerLenMax'],
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': thePara['TmMin'],
    'PRIMER_MAX_TM': thePara['TmMax'],
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 10,
    'PRIMER_INTERNAL_MAX_POLY_X': 10,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
    'PRIMER_PRODUCT_SIZE_RANGE': [[thePara['MinAmpLen']-thePara['PrimerLenMax'],thePara['MaxAmpLen']+thePara['PrimerLenMax']]],
    'PRIMER_TASK': 'generic',
    'PRIMER_PICK_LEFT_PRIMER': 1,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_PICK_RIGHT_PRIMER': 1,
    'PRIMER_PAIR_MAX_DIFF_TM': thePara['TmDeltra'],
}
primer3.bindings.setP3Globals(Primer3GlobalArgs)

for ChrID in InData.keys() :
    for thePos in InData[ChrID].keys() :
        FulChrID: str = ''.join(['chr',ChrID])
        # Start attributes are 1-based
        Left: int = thePos - thePara['Mode2LeftMax'] - thePara['PrimerLenMax'] -1
        if Left < 0 : Left = 0
        #Left = thePos-1
        # End attributes are 0-based
        Right: int = thePos + thePara['MaxAmpLen'] + thePara['PrimerLenMax']
        if Right > len(RefSeqs[FulChrID]) : Right = len(RefSeqs[FulChrID])
        theSeq: str = RefSeqs[FulChrID][Left:Right]
        print(':'.join([ChrID,str(thePos),FulChrID,str(theSeq),str(InData[ChrID][thePos]) ]))
        Primer3Ret: Dict = primer3.bindings.designPrimers({
            'SEQUENCE_ID': theSeq.fancy_name,
            'SEQUENCE_TEMPLATE': str(theSeq),
            'SEQUENCE_INCLUDED_REGION': [ thePara['PrimerLenMax'],thePara['MaxAmpLen'] ],
        })
        print(Primer3Ret)


print(b'[!] %(skipped)d InDels skipped in %(Total)d items.' % {b'skipped': Skipped, b'Total': Total})
