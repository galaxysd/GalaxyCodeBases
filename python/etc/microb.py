#!/usr/bin/env python3

import vcf
import gzip,csv

f = gzip.open('snpcandidatforpcr.out.gz', 'rt')
tsvin = csv.reader(f, delimiter='\t')
for row in tsvin:
    if len(row)==0: continue
    print(row)
    ChrID = row[0].split('.')[1]
    thePos = row[2]
    theRef = row[3]
    Alleles = row[4].split(',')
    Allele1 = Alleles[0].split('|')
    Allele2 = Alleles[1].split('|')
    print((ChrID,thePos,theRef))
    print(Allele1)
    print(Allele2)
    ...
f.close()

vcf_reader = vcf.Reader(filename='microb.vcf.gz',compressed=True)
for record in vcf_reader:
    print(record)
    print(record.genotype('ERR589860'))
    for sample in record.samples:
        print(sample)
