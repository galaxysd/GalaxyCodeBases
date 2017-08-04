#!/bin/bash

#Usage: bash /home/danielho/virus_clip.sh $1 $2 $3 $4
#Example: bash /home/danielho/virus_clip.sh /home/danielho/seq /home/danielho/virus_clip 100T 2

#tools
bwa=/software/sequencing/bwa_latest/bwa
samtools=/software/sequencing/samtools_latest/samtools
virus_clip=/home/danielho/virus_clip.pl
blastn=/home/danielho/blast/ncbi-blast-2.2.30+/bin/blastn
annovar=/pathool01/disk1/data/danielho/ireneng/ANNOVAR/annovar/annotate_variation.pl

#directory
SEQdir=$1
OUTdir=$2

#resources
refFasta=/pathool01/disk1/data/danielho/ireneng/ref_annotation/HBV/HBV.fa
blastdb=/pathool01/disk1/software/ref/hg19
annovardb=/pathool01/disk1/data/danielho/ireneng/ANNOVAR/annovar/humandb

#parameter
seq_file_suffix=.fastq.gz

#procedures
file=$3
library=$4	#1=single-end, 2=paired-end

#alignment
if [ "$library" -eq 1 ]; then	#single-end data
$bwa mem ${refFasta} ${SEQdir}/${file}${seq_file_suffix} > ${OUTdir}/align/${file}.bwa.sam
else				#paired-end data
$bwa mem ${refFasta} ${SEQdir}/${file}_1${seq_file_suffix} ${SEQdir}/${file}_2${seq_file_suffix} > ${OUTdir}/align/${file}.bwa.sam
fi

#extract softclip reads
$samtools view -Sh ${OUTdir}/align/${file}.bwa.sam | awk '$6 ~ /S/ {print $2"\t"$4"\t"$5"\t"$6"\t"$10}' > ${OUTdir}/align/${file}.softclip.txt

#detect human and virus integration breakpoint and sequence
mkdir ${OUTdir}/result/${file}
cd ${OUTdir}/result/${file}
perl ${virus_clip} ${OUTdir}/align/${file}.softclip.txt ${blastn} ${blastdb} ${annovar} ${annovardb}
