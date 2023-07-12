#!/usr/bin/env python3

import sys
import io
import os
import tqdm
import functools
import gtfparse
from collections import defaultdict

import logging
#logging.getLogger("gtfparse").setLevel(logging.WARNING)
#gtfparse.logger.setLevel(logging.WARNING)
for handler in logging.root.handlers:
    handler.setLevel(logging.WARNING)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parse_with_polars_lazy(
        filepath_or_buffer,
        split_attributes=True,
        features=None,
        fix_quotes_columns=["attribute"]):
    # use a global string cache so that all strings get intern'd into
    # a single numbering system
    polars.toggle_string_cache(True)
    kwargs = dict(
        has_header=False,
        separator="\t",
        comment_char="#",
        null_values=".",
        dtypes={
            "seqname": polars.Categorical, 
            "source": polars.Categorical, 
            "start": polars.Int64,
            "end": polars.Int64,
            "score": polars.Float32,
            "feature": polars.Categorical, 
            "strand": polars.Categorical, 
            "frame": polars.UInt32,
        })
    try:
        if type(filepath_or_buffer) is StringIO:
            df = polars.read_csv(
                filepath_or_buffer,
                new_columns=REQUIRED_COLUMNS,
                **kwargs).lazy()
        elif filepath_or_buffer.endswith(".gz") or filepath_or_buffer.endswith(".gzip"):
            with gzip.open(filepath_or_buffer) as f:
                df = polars.read_csv(
                    f,
                    new_columns=REQUIRED_COLUMNS,
                    **kwargs).lazy()
        else:
            df = polars.scan_csv(
                filepath_or_buffer, 
                with_column_names=lambda cols: REQUIRED_COLUMNS,
                **kwargs).lazy()
    except polars.ShapeError:
        raise ParsingError("Wrong number of columns")
    df = df.with_columns([
        polars.col("frame").fill_null(0),
        polars.col("attribute").str.replace_all('"', "'")
    ])
    for fix_quotes_column in fix_quotes_columns:
        # Catch mistaken semicolons by replacing "xyz;" with "xyz"
        # Required to do this since the Ensembl GTF for Ensembl
        # release 78 has mistakes such as:
        #   gene_name = "PRAMEF6;" transcript_name = "PRAMEF6;-201"
        df = df.with_columns([
            polars.col(fix_quotes_column).str.replace(';\"', '\"').str.replace(";-", "-")
        ])
    if features is not None:
        features = sorted(set(features))
        df = df.filter(polars.col("feature").is_in(features))
    if split_attributes:
        df = df.with_columns([
            polars.col("attribute").str.split(";").alias("attribute_split")
        ])
    return df
# https://github.com/openvax/gtfparse/pull/35
gtfparse.parse_with_polars_lazy = parse_with_polars_lazy

def main():
    if len(sys.argv) < 2 :
        print('Usage:',sys.argv[0],'<gtf file> >geneRange.out',file=sys.stderr,flush=True);
        exit(0);
    gtfFile = sys.argv[1]   # 'GCF_000001405.40_GRCh38.p14_genomic.gtf', 'h200.p14_genomic.gtf'
    Genes = defaultdict(functools.partial(defaultdict, list))
    with tqdm.tqdm(total=os.path.getsize(gtfFile)) as pbar:
        i = 0
        with open(gtfFile,'r') as gtfileh:
            #for line in gtfileh:
            while line:= gtfileh.readline():
            #for i, line in enumerate(gtfileh):
                #i += len(line)
                i += 1
                if not i % 1000:
                    pbar.update(gtfileh.tell() - pbar.n)
                    #pbar.update(i - pbar.n)
                    #eprint(i)
                gtfline = gtfparse.parse_gtf_and_expand_attributes(io.StringIO(line), restrict_attribute_columns=['gene_name'])
                if 'gene_name' in gtfline:
                    for k in ('start','end'):
                        Genes['\t'.join((gtfline['gene_name'][0],gtfline['strand'][0],gtfline['seqname'][0]))][k].append(gtfline[k][0])
                else:
                    print(str(gtfline))
                pbar.update(len(line))
        for k in sorted(Genes.keys()):
            Genes[k]['MinStart'] = min(Genes[k]['start'])
            Genes[k]['MaxEnd'] = max(Genes[k]['end'])
            #print(Genes[k])
            print('\t'.join((k,str(Genes[k]['MinStart']),str(Genes[k]['MaxEnd']))))
            print('\t'.join(('#',str(Genes[k]['start']),str(Genes[k]['end']))))
            sys.stdout.flush()
        #print(Genes)

if __name__ == "__main__":
    main()  # time ./gtfGeneRange.py h200.p14_genomic.gtf
