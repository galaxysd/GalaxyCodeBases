#!/usr/bin/env python3

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

import zstandard as zstd
import io

path = 'PRJNA380127.gbk.zst'
with open(path, 'rb') as fh:
    dctx = zstd.ZstdDecompressor()
    stream_reader = dctx.stream_reader(fh)
    text_stream = io.TextIOWrapper(stream_reader, encoding='utf-8')
    for record in SeqIO.parse(text_stream, "genbank"):
        print(record.id)

exit()

# https://www.ncbi.nlm.nih.gov/nuccore?term=380127%5BBioProject%5D

# get all sequence records for the specified genbank file
recs = [rec for rec in SeqIO.parse("PRJNA380127.gbk", "genbank")]

# print the number of sequence records that were extracted
print(len(recs))

# print annotations for each sequence record
for rec in recs:
	print(rec.annotations)

# print the CDS sequence feature summary information for each feature in each
# sequence record
for rec in recs:
    feats = [feat for feat in rec.features if feat.type == "CDS"]
    for feat in feats:
        print(feat)
