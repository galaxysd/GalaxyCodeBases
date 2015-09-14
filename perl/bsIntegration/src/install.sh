#!/bin/sh

cd bwa; make
cp -av bwa ../../bin/
cd ..

cd samtools; make
cp -av samtools ../../bin/
cd ..

cd idba
aclocal
autoconf
automake --add-missing
./configure
#make clean
make -j
cp -av bin/idba_hybrid ../../bin/
cd ..

cd methlyAln/bin; cmake ..
cp -av alnmethly ../../../bin/
cd ../..

ls -l ../bin/
