#!/usr/bin/env python3

import vcf
import gzip
f = gzip.open('snpcandidatforpcr.out.gz', 'rb')
for line in f.readlines():
    ...
f.close()

vcf_reader = vcf.Reader(filename='microb.vcf.gz',compressed=True)
for record in vcf_reader:
    print(record)
    print(record.genotype('ERR589860'))
    for sample in record.samples:
        print(sample)
