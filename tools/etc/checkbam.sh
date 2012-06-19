#!/bin/bash

if [ x"$1" = x"" ]; then
   echo "Usage: $0 <outfile> [find_path(.)]"
   exit
fi

if [ x"$2" = x"" ]; then
   TOFIND="."
else
   TOFIND="$2"
fi

date > $1
echo "From [${TOFIND}] to [$1]:"
echo "Find in [${TOFIND}]" >> $1

find ${TOFIND} -name '*.bam' | while read a;
 do echo -n "Checking [$a] ...";
    RES=`gzip -t "$a" 2>&1`;
    if [ "$?" -ne "0" ]; then
       echo -e "\b\b\b\b: Error=${?}.${RES}"
       echo "[e${?}] ${a}${RES}" >> $1
    else
       echo -e "\b\b\b\b: OK."
    fi;
 done

