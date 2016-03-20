#!/bin/bash

DIR=$(/usr/bin/dirname "${0}")

mystart ()
{
	cd ${DIR}
	php -q bgprocess.php </dev/null >logfile.log &
	echo $! > bg.pid
}

case "$1" in
	start)
		mystart
	;;
	stop)
		cd ${DIR}
		mypid=`cat bg.pid`
		var2=`ps xu | grep $mypid | wc -l`
		if [ "$var2" -gt 0 ]
		then
			kill $mypid
		fi
		;;
	check)
		# check if the process is running.
		cd ${DIR}
		mypid=`cat bg.pid`
		var2=`ps xu | grep $mypid | wc -l`
		if [ "$var2" -lt 1 ]
		then
			mystart
		fi
		;;
esac

exit 0;
