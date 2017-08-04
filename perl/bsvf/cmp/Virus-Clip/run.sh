#!/bin/bash

#Usage: bash /home/danielho/virus_clip.sh $1 $2 $3 $4
#Example: bash /home/danielho/virus_clip.sh /home/danielho/seq /home/danielho/virus_clip 100T 2

#tools
bwa=/ifs2/FGI/brew/bin/bwa
samtools=/ifs2/FGI/brew/bin/samtools
virus_clip=/share/FGI2017B/users/huxs/bsvf/cmp/Virus-Clip/virus_clip.pl
blastn=/share/FGI2017B/users/huxs/bsvf/cmp/pkg/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastn
annovar=/share/FGI2017B/users/huxs/bsvf/cmp/annovar/annotate_variation.pl

#directory
SEQdir=/share/FGI2017B/users/huxs/bsvf/cmp/Virus-Clip/fq
OUTdir=/share/FGI2017B/users/huxs/bsvf/cmp/Virus-Clip/out

#resources
refFasta=/share/FGI2017B/users/huxs/bsvf/cmp/Virus-Clip/ref/hbv.fa
blastdb=/share/FGI2017B/users/huxs/bsvf/cmp/Virus-Clip/ref/hg38.fa
annovardb=/share/FGI2017B/users/huxs/bsvf/cmp/annovar/humandb

#parameter
seq_file_suffix=.fq.gz

#procedures
file=WH1705004020
library=2	#1=single-end, 2=paired-end

#alignment
$bwa mem ${refFasta} ${SEQdir}/${file}_R1${seq_file_suffix} ${SEQdir}/${file}_R2${seq_file_suffix} > ${OUTdir}/align/${file}.bwa.sam

#extract softclip reads
$samtools view -Sh ${OUTdir}/align/${file}.bwa.sam | awk '$6 ~ /S/ {print $2"\t"$4"\t"$5"\t"$6"\t"$10}' > ${OUTdir}/align/${file}.softclip.txt

#detect human and virus integration breakpoint and sequence
mkdir ${OUTdir}/result/${file}
cd ${OUTdir}/result/${file}
perl ${virus_clip} ${OUTdir}/align/${file}.softclip.txt ${blastn} ${blastdb} ${annovar} ${annovardb}
