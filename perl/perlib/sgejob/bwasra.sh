#!/bin/sh
#$ -S /bin/bash

# --       project & hard_queue_list       --
#$ -q fgi.q -P fgimem

# --             our name                  --
#$ -N PnBwa

# --        virtual_free request           --
#$ -l vf=12G,num_proc=8

#$ -cwd
# rerun if  the  job  was aborted without leaving a consistent exit state.  (This is typically the case if  the  node  on  which  the  job  is  running crashes).
#$ -r y

# redefines the environment variables to be exported to the execution context of the job.
# cannot use ['] or [\'], or all after become 1 value.
#$ -v PERL5LIB,PATH,LD_LIBRARY_PATH,ENVA=TESTtest,ENVB="test\n rest"

#$ -t 1-81

INFILE=`sed -n "${SGE_TASK_ID}p" wgt.lst`
read -ra INDAT <<<"$INFILE"

mkdir -p bam/${INDAT[0]}
bwa mem -t 16 -Y ../ref/PanTigGCFv1.fna.gz -R "@RG\tID:${INDAT[2]}\tSM:${INDAT[0]}" ex/${INDAT[1]}/${INDAT[2]}_1.fastq.gz ex/${INDAT[1]}/${INDAT[2]}_2.fastq.gz 2>bam/${INDAT[0]}/${INDAT[2]}.log | samtools view -bS - | samtools sort -m 2G -T bam/${INDAT[0]}/${INDAT[2]}.tmp -o bam/${INDAT[0]}/${INDAT[2]}.bam
#samtools index bam/${INDAT[0]}/${INDAT[2]}.bam

