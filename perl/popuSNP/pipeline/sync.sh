#!/bin/env bash
if [ "$1" = "" ]; then
	s="/share/raid007/solexa-work/Index/PAN*"
else s="$@"
fi
echo Source from $s
mkdir -p ./0rawfq
date > ./0rawfq/_sync.log
echo Source from $s >> ./0rawfq/_sync.log
cp -sra $s ./0rawfq/ 2>> ./0rawfq/_sync.log
date >> ./0rawfq/_sync.log
find `readlink -nf ./0rawfq` -name "*.fq" > rawfq.lst
echo done !
