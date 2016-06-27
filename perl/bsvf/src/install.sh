#!/bin/sh

if [ ! -d "./idba" ]; then
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

cd idba
aclocal
autoconf
automake --add-missing
if [[ "$uname_S" == "Darwin" ]]; then
	gcc5_S=`which g++-5 || echo not`
	if [[ "$gcc5_S" == "not" ]]; then
		brew install --without-multilib gcc
	fi
	CXX=g++-5 ./configure
else
	./configure
fi
#make clean
make -j
cp -av bin/idba_hybrid ../../bin/
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
