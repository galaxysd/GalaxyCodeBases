#!/usr/bin/env python3

InFile = 'out.kmers.1k'

import csv
import primer3
from itertools import product

for i in product('ACTG', repeat = 18):
    oneKmer = ''.join(i)
    print(oneKmer)

exit()

with open(InFile) as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter=" ")
    for line in tsvreader:
        print(line[0])
        print(line[1])


t=primer3.calcHairpin('CCCCCATCCGATCAGGGGG')

print(t)
