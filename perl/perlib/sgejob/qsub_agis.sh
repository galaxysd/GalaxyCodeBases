#!/bin/sh
#$ -S /bin/bash

#$ -N getMatrix
#$ -l vf=3G -t 1-50

#$ -cwd
#$ -r y
#$ -v PERL5LIB,PATH,LD_LIBRARY_PATH

#$ -e /dev/null
#$ -o /dev/null

ARC=lx-amd64
QSUB=$SGE_ROOT/bin/$ARC/qsub
SLEEP=20

cmd="$QSUB $0 $0 $@"

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
	if [ "$RESTARTED" = "0" ]; then
		echo \# Begin  @ `date` >${MAIN}.err
		cat /dev/null > ${MAIN}.out
	else
		echo >>${MAIN}.err
		echo >>${MAIN}.out
		echo \#Restart @ `date` [$RESTARTED]>>${MAIN}.err
	fi
	echo \#$ENVIRONMENT $JOB_NAME of $QUEUE @ Host:$HOSTNAME as Job:[$JOB_ID],Task:[$SGE_TASK_ID] >>${MAIN}.err
# jobs here
	uname -a >>${MAIN}.out 2>>${MAIN}.err
	#perl -e 'for (1..20000000) {$a{$_}=$_;print $a{$_},"\n";};exit 1;'
SEEDFILE=id.lst
SEED=$(/bin/sed -n -e "${SGE_TASK_ID} p" $SEEDFILE)
echo [${SEED}]
perl ./new.Matrix_calculator.pl -b -r ref.fa -i bam/${SEED}.sort.bam -s vcf/${SEED}.SNPs.filter.vcf -o matrix/${SEED}.matrix 2>matrix/${SEED}.matrix.err
# jobs end
	ENDVALUE=$?
	cat <<EOFSTAT >> ${MAIN}.err
#  End  @ `date`
# Used $SECONDS sec.
#Job ended as $ENDVALUE
#\$PWD is $PWD

#\$SGE_O_WORKDIR is [$SGE_O_WORKDIR]
#PATH is [$PATH]
#LD_LIBRARY_PATH is [$LD_LIBRARY_PATH]
#PERL5LIB is [$PERL5LIB]
#PYTHONPATH is [$PYTHONPATH]

#`qstat -j $JOB_ID | grep usage`
#I am [$0]
#I was [$@]
#All done !
EOFSTAT
fi
