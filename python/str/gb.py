#!/usr/bin/env python3

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

import zstandard as zstd
import io

from collections import defaultdict
seqDat = defaultdict(list)

path = 'PRJNA380127.gbk.zst'
with open(path, 'rb') as fh:
    dctx = zstd.ZstdDecompressor()
    stream_reader = dctx.stream_reader(fh)
    text_stream = io.TextIOWrapper(stream_reader, encoding='utf-8')
    for record in SeqIO.parse(text_stream, "genbank"):
        print(record.id)
        print(record.description)
        datHumanSTR = record.annotations['structured_comment']['HumanSTR']
        print(datHumanSTR['STR locus name'])
        print(datHumanSTR['Length-based allele'])
        print(record.annotations['structured_comment']['HumanSTR']['Bracketed repeat'])
        print(datHumanSTR['RefSeq Accession'])
        print(datHumanSTR['Chrom. Location'])
        print(record.seq)
        print('---')
        seqDat[datHumanSTR['STR locus name']].append(['|'.join([record.id, record.description.replace('Homo sapiens microsatellite ','').replace(' ','_')]),record.seq])

print(seqDat)
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

"""
for record in SeqIO.parse("t.gb","genbank"):
    print(record.id)



>>> print(record.annotations['structured_comment']['HumanSTR'])
OrderedDict([('STR locus name', 'TPOX'), ('Length-based allele', '6'), ('Bracketed repeat', '[AATG]6'), ('Sequencing technology', 'ForenSeq, MiSeq FGx; PowerSeq Auto, MiSeq'), ('Coverage', '>30X'), ('Length-based tech.', 'PowerPlex Fusion, ABI3500xl'), ('Assembly', 'GRCh38 (GCF_000001405)'), ('Chromosome', '2'), ('RefSeq Accession', 'NC_000002.12'), ('Chrom. Location', '1489532..1489698'), ('Repeat Location', '1489653..1489684'), ('Cytogenetic Location', '2p25.3')])

>>> print(record)
ID: MF044246.1
Name: MF044246
Description: Homo sapiens microsatellite TPOX 6 [AATG]6 sequence
Database cross-references: BioProject:PRJNA380554
Number of features: 4
/molecule_type=DNA
/topology=linear
/data_file_division=PRI
/date=04-SEP-2018
/accessions=['MF044246']
/sequence_version=1
/keywords=['STRSeq, STR, TPOX']
/source=Homo sapiens (human)
/organism=Homo sapiens
/taxonomy=['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini', 'Catarrhini', 'Hominidae', 'Homo']
/references=[Reference(title='STRSeq: A catalog of sequence diversity at human identification Short Tandem Repeat loci', ...), Reference(title='Direct Submission', ...)]
/comment=Annotation ('bracketing') of the repeat region is consistent with
the guidance of the ISFG (International Society of Forensic
Genetics), PMID: 26844919.  Lower case letters in the 'Bracketed
repeat' region below denote uncounted bases.   The given
length-based allele value was determined using the designated
length-based technology.  Variation in the length-based allele
between individuals or assays can result from indels in flanking
regions.  The length of reported sequence is dependent on the assay
and the quality of the flanking sequence. This information is
provided as part of the STR Sequencing Project (STRseq), a
collaborative effort of the international forensic DNA community.
The purpose of this project is to facilitate the description of
sequence-based STR alleles.  Additional resources can be found at
strseq.nist.gov.  For questions or feedback, please contact
strseq@nist.gov.  Allele frequency data can be accessed in the
strider.online database.
/structured_comment=OrderedDict([('HumanSTR', OrderedDict([('STR locus name', 'TPOX'), ('Length-based allele', '6'), ('Bracketed repeat', '[AATG]6'), ('Sequencing technology', 'ForenSeq, MiSeq FGx; PowerSeq Auto, MiSeq'), ('Coverage', '>30X'), ('Length-based tech.', 'PowerPlex Fusion, ABI3500xl'), ('Assembly', 'GRCh38 (GCF_000001405)'), ('Chromosome', '2'), ('RefSeq Accession', 'NC_000002.12'), ('Chrom. Location', '1489532..1489698'), ('Repeat Location', '1489653..1489684'), ('Cytogenetic Location', '2p25.3')]))])
Seq('TGGCCTGTGGGTCCCCCCATAGATCGTAAGCCCAGGAGGAAGGGCTGTGTTTCA...AAA', IUPACAmbiguousDNA())

LOCUS       MF044246                 159 bp    DNA     linear   PRI 04-SEP-2018
DEFINITION  Homo sapiens microsatellite TPOX 6 [AATG]6 sequence.
ACCESSION   MF044246
VERSION     MF044246.1
DBLINK      BioProject: PRJNA380554
KEYWORDS    STRSeq, STR, TPOX.
SOURCE      Homo sapiens (human)
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;
            Catarrhini; Hominidae; Homo.
REFERENCE   1  (bases 1 to 159)
  AUTHORS   Gettings,K.B., Borsuk,L.A., Ballard,D., Bodner,M., Budowle,B.,
            Devesse,L., King,J., Parson,W., Phillips,C. and Vallone,P.M.
  TITLE     STRSeq: A catalog of sequence diversity at human identification
            Short Tandem Repeat loci
  JOURNAL   Forensic Sci Int Genet 31, 111-117 (2017)
   PUBMED   28888135
REFERENCE   2  (bases 1 to 159)
  AUTHORS   NIST,A.G.G.
  TITLE     Direct Submission
  JOURNAL   Submitted (04-MAY-2017) Applied Genetics Group, National Institute
            of Standards and Technology, 100 Bureau Drive, MS-8314,
            Gaithersburg, MD 20899, USA
COMMENT     Annotation ('bracketing') of the repeat region is consistent with
            the guidance of the ISFG (International Society of Forensic
            Genetics), PMID: 26844919.  Lower case letters in the 'Bracketed
            repeat' region below denote uncounted bases.   The given
            length-based allele value was determined using the designated
            length-based technology.  Variation in the length-based allele
            between individuals or assays can result from indels in flanking
            regions.  The length of reported sequence is dependent on the assay
            and the quality of the flanking sequence. This information is
            provided as part of the STR Sequencing Project (STRseq), a
            collaborative effort of the international forensic DNA community.
            The purpose of this project is to facilitate the description of
            sequence-based STR alleles.  Additional resources can be found at
            strseq.nist.gov.  For questions or feedback, please contact
            strseq@nist.gov.  Allele frequency data can be accessed in the
            strider.online database.

            ##HumanSTR-START##
            STR locus name        :: TPOX
            Length-based allele   :: 6
            Bracketed repeat      :: [AATG]6
            Sequencing technology :: ForenSeq, MiSeq FGx; PowerSeq Auto, MiSeq
            Coverage              :: >30X
            Length-based tech.    :: PowerPlex Fusion, ABI3500xl
            Assembly              :: GRCh38 (GCF_000001405)
            Chromosome            :: 2
            RefSeq Accession      :: NC_000002.12
            Chrom. Location       :: 1489532..1489698
            Repeat Location       :: 1489653..1489684
            Cytogenetic Location  :: 2p25.3
            ##HumanSTR-END##
FEATURES             Location/Qualifiers
     source          1..159
                     /organism="Homo sapiens"
                     /mol_type="genomic DNA"
                     /db_xref="taxon:9606"
     misc_feature    1..159
                     /note="Promega PowerSeq Sequence"
     misc_feature    120..150
                     /note="Illumina ForenSeq Sequence"
     repeat_region   122..145
                     /rpt_type=tandem
                     /satellite="microsatellite:TPOX"
ORIGIN
        1 tggcctgtgg gtccccccat agatcgtaag cccaggagga agggctgtgt ttcagggctg
       61 tgatcactag cacccagaac cgtcgactgg cacagaacag gcacttaggg aaccctcact
      121 gaatgaatga atgaatgaat gaatgtttgg gcaaataaa
//
"""
