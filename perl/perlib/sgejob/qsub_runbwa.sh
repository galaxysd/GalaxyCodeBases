#!/bin/bash
#$ -S /bin/bash
#$ -q fgi.q -P fgi
#$ -l vf=8G,num_proc=9
#$ -binding linear:12
#$ -cwd -r y
#$ -v PERL5LIB,PATH,LD_LIBRARY_PATH
#$ -e /dev/null -o /dev/null

#$ -N p2bwa -hold_jid p1cutmgiadapter

SEEDFILE=fqprefix.lst
SEEDLINES=$(/usr/bin/wc -l $SEEDFILE|/bin/awk '{print $1}')

ARC=lx-amd64
QSUB=$SGE_ROOT/bin/$ARC/qsub
SLEEP=20

cmd="$QSUB -t 1-$SEEDLINES $0 $0 $@"

BWAREF="./ref/Homo_sapiens_assembly38.fasta"

# started by SGE or manually
if [ "$JOB_ID" = "" ]; then
	echo "submitting $0 ..."
	$cmd
	while [ "x$?" != "x0" ]; do
	echo "$0: qsub failed - retrying .." >&2
		sleep $SLEEP
		$cmd
	done
else
	if [ "$1" = "" ]; then
		MAIN=${JOB_NAME}_${JOB_ID}${SGE_TASK_ID}
	else
		MAIN=$1
	fi
	echo \#Running @ Host:$HOSTNAME as Job:[$JOB_ID],Task:[$SGE_TASK_ID] >${MAIN}.err
	echo \#Begin @ `date` >>${MAIN}.err
# jobs here
	uname -a > ${MAIN}.out 2>>${MAIN}.err
	sleep 2
	ARG1=$(/bin/sed -n -e "${SGE_TASK_ID} p" $SEEDFILE|/bin/awk -F"\t" '{print $1}')
	bwa mem -t 12 -Y ${BWAREF} -R "@RG\tID:${ARG1}\tSM:${ARG1}" -p fq/${ARG1}.fq.gz 2>bam/${ARG1}.log | samtools view -bS - | samtools sort -m 2G -T bam/${ARG1}.tmp -o bam/${ARG1}.bam
# jobs end
	ENDVALUE=$?
	cat <<EOFSTAT >> ${MAIN}.err
#  End @ `date`
# Used $SECONDS sec.
#Job ended as $ENDVALUE
#\$PWD is [$PWD]
#\$SGE_O_WORKDIR is [$SGE_O_WORKDIR]

ENVIRONMENT is [$ENVIRONMENT] BATCH or not
JOB_NAME is [$JOB_NAME] and REQUEST is [$REQUEST]
PATH is [$PATH]
QUEUE is [$QUEUE]
RESTARTED is [$RESTARTED] 0 or 1
TMPDIR is [$TMPDIR] and TMP is [$TMP]

PERL5LIB is [$PERL5LIB]
PYTHONPATH is [$PYTHONPATH]
HOME is [$HOME] 
`which perl`	is works after -v PATH

#I am [$0]
#I was [$@]
#All done !
EOFSTAT
fi
