#!/usr/bin/env python3

import os
import csv

infile = 'sites.txt'
bases = ['A','C','G','T']
header = ['rsID','Chr','Pos']
headerstr = '\t'.join(header+bases+['sum','dbg'])
print(headerstr)

with open( infile, 'r' ) as theFile:
    reader = csv.DictReader(theFile,dialect=csv.excel_tab)
    for line in reader:
        # line is { 'workers': 'w0', 'constant': 7.334, 'age': -1.406, ... }
        # e.g. print( line[ 'workers' ] ) yields 'w0'
        # print(line)
        freqs = dict.fromkeys(bases, 0)
        freqs[line['Alt1']] = line['Freq1']
        freqs[line['Alt2']] = line['Freq2']
        freqs[line['Alt3']] = line['Freq3']
        dbgstr = ','.join([line['Alt1'],line['Freq1'],line['Alt2'],line['Freq2'],line['Alt3'],line['Freq3']])
        outlist = [line['rsID'],line['Chr'],line['Pos']]
        sumfreq = 0.0
        for e in bases:
            sumfreq += float(freqs[e])
            outlist.append(str(freqs[e]))
        outlist.append(str(sumfreq))
        outlist.append(dbgstr)
        outstr = '\t'.join(outlist)
        if sumfreq > 0.1:
            print(outstr)
