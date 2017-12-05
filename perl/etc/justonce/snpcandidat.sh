#!/bin/bash
# https://www.ncbi.nlm.nih.gov/books/NBK25500/
# Based on https://www.biostars.org/p/244479/#244495

mkdir t
set -u

#TAXLIST=`gzcat snpcandidatforpcr.out.gz|awk 'NF > 0 {split($0,a,".");print a[1]}'|uniq`

TAXLIST=(213810 245014 349741 411459 411460 411461 411462 411463 411469 411474 411477 411479 411481 411483 411484 411486 411489 411902 411903 428125 445970 445972 451640 457412 457421 469586 469587 470145 470146 471875 483216 483218 484018 500632 511680 515619 515620 518636 518637 537007 537011 537012 545696 547042 548480 552396 563191 563193 566550 585057 592028 622312 657313 657314 657315 657319 657321 657322 657323 679189 702443 702446 717959 717960 717962 718252)

for TAX in "${TAXLIST[@]}" ; do
	echo getting genome for: $TAX
	GENOME=$(esearch -db genome -query txid"$TAX"[Organism:noexp] | efetch -format docsum | tee "./t/${TAX}.genome.esearch.docsum")
	ACC=`echo $GENOME | xtract -pattern DocumentSummary  -element Assembly_Accession`
	NAME=`echo $GENOME |  xtract -pattern DocumentSummary -element Assembly_Name`
	echo authoritative genome: $ACC $NAME
	RESULT=$(esearch -db assembly -query "$ACC" |
		efetch -format docsum | tee "${TAX}.assembly.esearch.docsum")
	FTPP=`echo $RESULT | xtract -pattern DocumentSummary  -element FtpPath_GenBank`
	TAXID=`echo $RESULT | xtract -pattern DocumentSummary  -element Taxid`
	echo FtpPath: $FTPP
	BASENAME=`basename $FTPP`
	FTPPATHG=$FTPP/$BASENAME'_genomic.fna.gz'
	FTPPATHP=$FTPP/$BASENAME'_genomic.gff.gz'
	# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/307/495/GCF_000307495.1_Para_merd_CL09T00C40_V1/GCF_000307495.1_Para_merd_CL09T00C40_V1_genomic.fna.gz
	# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/307/495/GCF_000307495.1_Para_merd_CL09T00C40_V1/GCF_000307495.1_Para_merd_CL09T00C40_V1_genomic.gff.gz
	echo Downloading $FTPPATHG ...

	## get genome data
	wget $FTPPATHG
	BASENAME=`basename $FTPPATHG`
	#gunzip -f $BASENAME
	#BASENAME=`echo $BASENAME | sed s/.gz//`
	#makeblastdb -in $BASENAME -dbtype nucl -parse_seqids -taxid $TAXID -title "$TAX $NAME genomic"
	echo Downloading $FTPPATHP ...
	## get protein data
	wget $FTPPATHP
	BASENAME=`basename $FTPPATHP`
	#gunzip -f $BASENAME
	#BASENAME=`echo $BASENAME | sed s/.gz//`
	#makeblastdb -in $BASENAME -dbtype prot -parse_seqids -taxid $TAXID -title "$TAX $NAME proteins"
done
