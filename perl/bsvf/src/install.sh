#!/bin/sh

if [ ! -d "./htslib" ]; then
	./download.sh
fi

cd bwa; make
cp -av bwa ../../bin/
cd ..

cd samtools; make
cp -av samtools ../../bin/
cd ..

uname_S=`uname -s 2>/dev/null || echo not`
has_cmake=`which cmake || echo not`

cd analyser
make
cp -av analyser/bsanalyser ../../bin/
cd ..

if [[ "$uname_S" == "Darwin" ]]; then
	cd methlyAln/src; make
	cp -av ../bin/alnmethly ../../../bin/
else
	cd methlyAln/bin
	if [[ "$has_cmake" == "not" ]]; then
		make -f Makefile.Linux
	else
		cmake ..
	fi
	cp -av alnmethly ../../../bin/
fi
cd ../..

ls -l ../bin/

echo "[!]To clean: rm -fr idba bwa samtools htslib"
