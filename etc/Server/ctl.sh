#!/bin/bash

MYPS=`ps -u huxs|grep ffmpeg|awk '{print $1}'`

#echo ${MYPS}

if [ "$1" = "s" ]; then
	kill -TSTP ${MYPS}
else
	kill -CONT ${MYPS}
fi

# $ crontab -l
# 24 8 * * *	/bak/seqdata/files/FL/ctl.sh s
# 12 0 * * *	/bak/seqdata/files/FL/ctl.sh go

