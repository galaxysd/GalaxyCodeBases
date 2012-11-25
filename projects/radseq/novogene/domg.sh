#!/bin/bash
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#$ -cwd -r y -l vf=1g,p=1
#$ -e /dev/null -o /dev/null
#$ -S /bin/bash -t 1-19

MY_HOST=`/bin/hostname`
MY_DATE=`/bin/date`
env
echo "================================================================"
echo JOB_NAME=$JOB_NAME
echo JOB_ID=$JOB_ID
echo SGE_TASK_ID=$SGE_TASK_ID
echo SGE_TASK_FIRST=$SGE_TASK_FIRST
echo SGE_TASK_LAST=$SGE_TASK_LAST
echo NSLOTS=$NSLOTS
echo QUEUE=$QUEUE
echo SGE_CWD_PATH=$SGE_CWD_PATH
echo PATH=$PATH
echo SGE_STDIN_PATH=$SGE_STDIN_PATH
echo SGE_STDOUT_PATH=$SGE_STDOUT_PATH
echo SGE_STDERR_PATH=$SGE_STDERR_PATH
echo SGE_O_HOST=$SGE_O_HOST
echo SGE_O_PATH=$SGE_O_PATH
echo "================================================================"

echo "Running job JOB_NAME=$JOB_NAME task SGE_TASK_ID=$SGE_TASK_ID on $MY_HOST at $MY_DATE"

SEEDFILE=ng/mg.lst
SEED=$(/bin/sed -n -e "${SGE_TASK_ID} p" $SEEDFILE)
PATH=$(/bin/sed -n -e "${SGE_TASK_ID} p" $SEEDFILE|/bin/awk -F"\t" '{print $1}')
MAIN=$(/bin/sed -n -e "${SGE_TASK_ID} p" $SEEDFILE|/bin/awk -F"\t" '{print $2}')
TYPE=$(/bin/sed -n -e "${SGE_TASK_ID} p" $SEEDFILE|/bin/awk -F"\t" '{print $3}')

echo "Running job JOB_NAME=$JOB_NAME task SGE_TASK_ID=$SGE_TASK_ID on $MY_HOST at $MY_DATE" >$PATH/$MAIN.bamlog

echo $PATH/*$MAIN*.sort.bam >>$PATH/$MAIN.bamlog
echo /usr/local/bin/samtools merge $PATH/${TYPE}_$MAIN.merge.bam $PATH/*$MAIN*.sort.bam >>$PATH/$MAIN.bamlog
/bin/date >> $PATH/$MAIN.bamlog
echo /usr/local/bin/samtools rmdup $PATH/${TYPE}_$MAIN.merge.bam $PATH/${TYPE}_$MAIN.rmdup.bam >>$PATH/$MAIN.bamlog
/bin/date >> $PATH/$MAIN.bamlog

