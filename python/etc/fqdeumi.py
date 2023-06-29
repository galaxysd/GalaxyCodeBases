#!/usr/bin/env python3

import sys
import os

UMI_LENGTH = 9

def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

try:
    fn = sys.argv[1]
except IndexError as ie:
    raise SystemError("Error: Specify file name\n")

if not os.path.exists(fn):
    raise SystemError("Error: File does not exist\n")

n = 4
with open(fn, 'r') as fh:
    lines = []
    for line in fh:
        lines.append(line.rstrip())
        if len(lines) == n:
            record = process(lines)
            #record['umi'] = record['sequence'][:9]
            record['name'] = ' '.join((record['name'],record['sequence'][:UMI_LENGTH]))
            record['sequence'] = record['sequence'][UMI_LENGTH:]
            sys.stderr.write("Record: %s\n" % (str(record)))
            lines = []
