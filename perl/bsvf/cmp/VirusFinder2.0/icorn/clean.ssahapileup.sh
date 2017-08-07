#!/bin/bash

root=$1
num=$2

#dleete tmp.*
rm tmp.*

# delete fastq files
for ((x=1;$x<=$num;x++)) do
   rm $root.$x;
   if [ ! -s "out.$x.e" ] ; then 
	   rm out.$x.e; 
	   rm out.$x.o; 
   else 
	   echo "Error in $x";
	   cat out.$x.e;
	   cat out.$x.o;
   fi 
done;
rm $root.fastq.names
if [ -f "$root.allDistance.txt" ] ; then 
 rm $root.allDistance.txt
fi

# zip all the stuff
gzip $root.cigar

### file used in Acorn
if [ -f "$root.cigarX" ] ; then
	gzip $root.cigarX
fi

gzip $root.AllMapped.cigar
gzip $root.pileup

if [ -f "$root.fastq" ] ; then
	rm $root.fastq
fi

rm ../*.tmp
rm ../*.qry.fa
rm ../*.ref.fa