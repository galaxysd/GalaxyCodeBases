#!/bin/bash
#$ -S /bin/bash
#$ -q fgi.q -P fgi
#$ -l vf=3G,num_proc=3
#$ -binding linear:4
#$ -cwd -r y
#$ -v PERL5LIB,PATH,LD_LIBRARY_PATH
#$ -e /dev/null -o /dev/null

#$ -N p1cutmgiadapter -hold_jid p0xxx

SEEDFILE=fqprefix.lst
SEEDLINES=$(/usr/bin/wc -l $SEEDFILE|/bin/awk '{print $1}')

ARC=lx-amd64
QSUB=$SGE_ROOT/bin/$ARC/qsub
SLEEP=20

cmd="$QSUB -t 1-$SEEDLINES $0 $0 $@"

SEADAPTERSTR="-g GAACGACATGGCTACGATCCGACTT -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAANNNNNNNNNNCAACTCCTTGGCTCACAGAACGACATGGCTACGATCCGACTT"
PEADAPTERSTR="${SEADAPTERSTR} -A AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTGNNNNNNNNNNTTGTCTTCCTAAGACCGCTTGGCCTCCGACTT -G TGTGAGCCAAGGAGTTGNNNNNNNNNNTTGTCTTCCTAAGACCGCTTGGCCTCCGACTT"
CUTADAPTARG="-e 0.1 -n 2 -m 1 -O 5 -j 4"

INPATH="./son"

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
	#ARG2=$(/bin/sed -n -e "${SGE_TASK_ID} p" $SEEDFILE|/bin/awk -F"\t" '{print $2}')
	cutadapt ${CUTADAPTARG} ${PEADAPTERSTR} --interleaved -o fq/${ARG1}.fq.gz ${INPATH}/${ARG1}_1.fq.gz ${INPATH}/${ARG1}_2.fq.gz >fq/${ARG1}.log
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
