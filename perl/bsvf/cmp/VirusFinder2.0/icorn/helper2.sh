#!/bin/bash

root=$1
awk '{ print $8}'  $root.ins > Ins.Reads.fofn
echo "$PILEUP_HOME/ssaha_pileup/other_codes/get_seqreads/get_seqreads Ins.Reads.fofn $root.fastq Ins.Reads.fastq"
$PILEUP_HOME/ssaha_pileup/other_codes/get_seqreads/get_seqreads Ins.Reads.fofn $root.fastq Ins.Reads.fastq
