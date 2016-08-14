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

cd emboss
#cd EMBOSS-6.6.0/
mv configure.in configure.ac
cat ../patches/*.patch | patch -p1

libtoolize --install --copy --force --automake
aclocal -I m4
autoconf --force
autoheader
automake --add-missing --copy --force-missing

# https://github.com/ebi-pf-team/interproscan/wiki/CompilingBinaries . Also Gentoo for `--enable-systemlibs --enable-static`
./configure --disable-debug --enable-64 --with-thread --without-x --without-java --without-mysql --without-postgresql --disable-shared --without-hpdf --without-pngdriver --enable-systemlibs --enable-static
make -j
cp -av emboss/water ../../bin/
cd ..

# if [[ "$uname_S" == "Darwin" ]]; then
# 	cd methlyAln/src; make
# 	cp -av ../bin/alnmethly ../../../bin/
# else
# 	cd methlyAln/bin
# 	if [[ "$has_cmake" == "not" ]]; then
# 		make -f Makefile.Linux
# 	else
# 		cmake ..
# 	fi
# 	cp -av alnmethly ../../../bin/
# fi
# cd ../..

ls -l ../bin/

#echo "[!]To clean: rm -fr idba bwa samtools htslib"
