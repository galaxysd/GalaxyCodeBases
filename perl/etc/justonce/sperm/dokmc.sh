#!/bin/bash
mkdir -p tmp/$1
./kmc -m5 -ci2 -cs200000000 -k25 $1 out/$1.k25 tmp/$1
rmdir tmp/$1
./kmc_dump out/$1.k25 /dev/stdout | gzip -9c > stat/$1.k25.gz

./getkmhist.pl stat/$1.k25.gz

