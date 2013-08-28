#!/usr/bin/make
#
# hash - combined shs and md5 makefile
#
# @(#) $Revision: 13.5 $
# @(#) $Id: Makefile,v 13.5 2010/10/12 21:08:53 chongo Exp $
# @(#) $Source: /usr/local/src/cmd/hash/RCS/Makefile,v $
#
# This file was written by Landon Curt Noll.
#
# This makefile has been placed in the public domain.  Please do not
# copyright this makefile.
#
# LANDON CURT NOLL DISCLAIMS ALL WARRANTIES WITH  REGARD  TO
# THIS  SOFTWARE,  INCLUDING  ALL IMPLIED WARRANTIES OF MER-
# CHANTABILITY AND FITNESS.  IN NO EVENT SHALL  LANDON  CURT
# NOLL  BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
# DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM  LOSS  OF
# USE,  DATA  OR  PROFITS, WHETHER IN AN ACTION OF CONTRACT,
# NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR  IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
# chongo (was here) /\oo/\
# http://www.isthe.com/chongo/index.html
#
# Share and enjoy! :-)
#
# See shsdrvr.c, shs1drvr.c and md5drvr.c for version and modification history.

# standard tools
#
SHELL= /bin/sh
TR= tr
DIFF= diff
INSTALL= install

# where to install binary files
#
BINDIR= /usr/local/bin
#BINDIR= /usr/bin
#BINDIR= /usr/contrib/bin

# where to install the shs, shs1 and md5 test file dirs
#
# ${TOPDIR} is the directory under which the shs and md5 dirs will be placed.
# ${SHSLIB} is where the readonly shs and shs1 test files are kept
# ${MD5LIB} is where the readonly md5 test files are kept
#
TOPDIR= /usr/local/lib
#TOPDIR= /usr/lib
#TOPDIR= /usr/libdata

SHSLIB= ${TOPDIR}/shs
MD5LIB= ${TOPDIR}/md5

# where man pages are installed
#
# Use MANDIR= to disable installation of the calc man (source) page.
#
#MANDIR=
MANDIR=/usr/local/man/man1
#MANDIR=/usr/man/man1
#MANDIR=/usr/share/man/man1
#MANDIR=/usr/man/u_man/man1
#MANDIR=/usr/contrib/man/man1

# extenstion to add on to the calc man page filename
#
# This is ignored if MANDIR is empty.
#
MANEXT= 1
#MANEXT= l

# select the proper compiler
#
# If you use a non-ANSI cc, you may have to edit md5.c to change the
# UL constant suffixes in MD5Transform() to just L and see if test works.
#
# Under BSD/OS 2.0, we cannot compile with -ansi because <sys/resource.h>
# uses struct timeval which is not defined by <sys/time.h> under -ansi.
#
CC= cc
#CC= gcc -ansi
#CC= gcc -ansi -fcaller-saves
#CC= gcc2 -ansi -fcaller-saves
#CC= gcc
#CC= gcc -fcaller-saves
#CC= gcc2 -fcaller-saves

# select the proper optimization or debug level
#
#OPT_DEBUG=
#OPT_DEBUG= -DDEBUG
#OPT_DEBUG= -g3
#OPT_DEBUG= -g3 -DDEBUG
#OPT_DEBUG= -g3 -DDEBUG -DDETAILED_DEBUG
#OPT_DEBUG= -O
#OPT_DEBUG= -O -DDEBUG
#OPT_DEBUG= -O -g3
#OPT_DEBUG= -O -g3 -DDEBUG
#OPT_DEBUG= -O2
#OPT_DEBUG= -O2 -DDEBUG
#OPT_DEBUG= -O2 -g3
#OPT_DEBUG= -O2 -g3 -DDEBUG
#OPT_DEBUG= -O3
#OPT_DEBUG= -O3 -DDEBUG
OPT_DEBUG= -O3 -g3
#OPT_DEBUG= -O3 -g3 -DDEBUG

CFLAGS= ${OPT_DEBUG} -Wall -W
#CFLAGS= ${OPT_DEBUG} -Wall -W -Werror

# Some architectures such as Sparc do not allow one to access long that
# are not alligned.  If make test dumpds core or fails to produce no
# output, use only -DMUST_ALIGN.
#
# ALIGN=		let must_align.c figure out if alignment is needed
# ALIGN=-DMUST_ALIGN 	force alignment (at the cost of performance}
# ALIGN=-UMUST_ALIGN	allow non-aligment (usually at a performance gain}
#
ALIGN=
#ALIGN=-DMUST_ALIGN
#ALIGN=-UMUST_ALIGN

all: shs sha shs1 sha1 md5

shs: shsdrvr.o shsio.o shsdual.o shs.o
	${CC} ${CFLAGS} -o shs shsdrvr.o shsio.o shsdual.o shs.o

sha: shs
	rm -f $@
	ln $? $@

shs1: shs1drvr.o shs1io.o shs1dual.o shs1.o
	${CC} ${CFLAGS} -o shs1 shs1drvr.o shs1io.o shs1dual.o shs1.o

sha1: shs1
	rm -f $@
	ln $? $@

md5: md5drvr.o md5io.o md5dual.o md5.o
	${CC} ${CFLAGS} -o md5 md5drvr.o md5io.o md5dual.o md5.o

# perform all tests
test: shs shs1 md5 shstest shs1test md5test

# perform the shs1 test
shs1test: shs1 file1 file2 file3 shs1.data shs1.ddata mega Makefile
	@echo This shs1 test should produce no output
	@./shs1 -x | ${DIFF} -u - shs1.data
	@./shs1 shs1.data | \
	    sed '/^407a65c1b4c95c5dd19da3c22b5d7357eadc2aeb shs1.data$$/d'
	@./mega | ./shs1 | \
	    sed '/^34aa973cd4c4daa4f61eeb2bdbad27316534016f$$/d'
	@./shs1 -q -c file3 | \
	    sed '/^0xa2ecae1924928cfebf7bbb01ee0974063e845824$$/d'
	@./shs1 -p abc -c file2 | \
	    sed '/^0xa2ecae1924928cfebf7bbb01ee0974063e845824 file2$$/d'
	@./shs1 -P file1 file2 | \
	    sed '/^a2ecae1924928cfebf7bbb01ee0974063e845824 file2$$/d'
	@./shs1 shs1.ddata | \
	    sed '/^a32cb7757f15c922f5466c279f12803d5848c21b shs1.ddata$$/d'
	@./shs1 -d -x -q | ${DIFF} -u - shs1.ddata
	@echo End of shs1 test

# perform the old extended shs digest test
shstest: shs file1 file2 file3 shs.data shs.ddata Makefile
	@echo This old shs test should produce no output
	@./shs -x | ${DIFF} -u - shs.data
	@./shs shs.data | \
	    sed '/^6789cc089b94c3a950745071a8ecf5b5b61bfc3c shs.data$$/d'
	@./shs -q -c file3 | \
	    sed '/^0xa9f2c8b901a7c8628133694a105fc760bdccad2c$$/d'
	@./shs -p abc -c file2 | \
	    sed '/^0xa9f2c8b901a7c8628133694a105fc760bdccad2c file2$$/d'
	@./shs -P file1 file2 | \
	    sed '/^a9f2c8b901a7c8628133694a105fc760bdccad2c file2$$/d'
	@./shs shs.ddata | \
	    sed '/^c0e05434119fec847d91fe981dafda58a7a97e5d shs.ddata$$/d'
	@./shs -d -x -q | ${DIFF} -u - shs.ddata
	echo End of old shs test

# perform the extended md5 digest test
md5test: md5 file1 file2 file3 md5.data md5.ddata Makefile
	@echo This md5 test should produce no output
	@./md5 -x | ${DIFF} -u - md5.data
	@./md5 md5.data | \
	    sed '/^9638425fd6b565e9fddf642fa4889125 md5.data$$/d'
	@./md5 -q -c file3 | \
	    sed '/^0x68cfe1a8e5242c8e1e8152d0bc710631$$/d'
	@./md5 -p abc -c file2 | \
	    sed '/^0x68cfe1a8e5242c8e1e8152d0bc710631 file2$$/d'
	@./md5 -P file1 file2 | \
	    sed '/^68cfe1a8e5242c8e1e8152d0bc710631 file2$$/d'
	@./md5 md5.ddata | \
	    sed '/^e1ac664dddab5c1b0d987a3d6ab0be9b md5.ddata$$/d'
	@./md5 -d -x -q | ${DIFF} -u - md5.ddata
	@echo End of md5 test

file1: Makefile
	rm -f file1
	echo abc | ${TR} -d '\012' > file1

file3: file1 file2
	rm -f file3
	cat file1 file2 > file3

shsdrvr.o: shsdrvr.c shs.h Makefile
	${CC} ${CFLAGS} -DTLIB=\"${SHSLIB}\" -c shsdrvr.c

shs1drvr.o: shs1drvr.c shs1.h Makefile
	${CC} ${CFLAGS} -DTLIB=\"${SHSLIB}\" -c shs1drvr.c

shsio.o: shsio.c shs.h
	${CC} ${CFLAGS} -c shsio.c

shs1io.o: shs1io.c shs1.h
	${CC} ${CFLAGS} -c shs1io.c

shsdual.o: shsdual.c shs.h Makefile
	${CC} ${CFLAGS} -DTLIB=\"${SHSLIB}\" -c shsdual.c

shs1dual.o: shs1dual.c shs1.h Makefile
	${CC} ${CFLAGS} -DTLIB=\"${SHSLIB}\" -c shs1dual.c

shs.o: shs.c shs.h align.h endian.h
	${CC} ${CFLAGS} -c shs.c

shs1.o: shs1.c shs1.h align.h endian.h
	${CC} ${CFLAGS} -c shs1.c

md5drvr.o: md5drvr.c md5.h Makefile
	${CC} ${CFLAGS} -DTLIB=\"${MD5LIB}\" -c md5drvr.c

md5io.o: md5io.c md5.h
	${CC} ${CFLAGS} -c md5io.c

md5dual.o: md5dual.c md5.h Makefile
	${CC} ${CFLAGS} -DTLIB=\"${MD5LIB}\" -c md5dual.c

md5.o: md5.c md5.h align.h endian.h
	${CC} ${CFLAGS} -c md5.c

align.h: must_align
	rm -f align.h
	-./must_align > align.h 2>/dev/null

must_align: must_align.c Makefile
	${CC} ${ALIGN} must_align.c -o must_align

endian.h: endian
	rm -f endian.h
	./endian > endian.h

endian: endian.c
	${CC} ${CFLAGS} endian.c -o endian

mega: mega.c
	${CC} ${CFLAGS} mega.c -o mega

install: all file1 file2 shs1.data shs.data shs1.ddata shs.ddata \
    md5.data shs.1 md5.1
	-rm -f ${BINDIR}/shs
	${INSTALL} -c -m 0555 shs ${BINDIR}
	-rm -f ${BINDIR}/sha
	${INSTALL} -c -m 0555 sha ${BINDIR}
	-rm -f ${BINDIR}/shs1
	${INSTALL} -c -m 0555 shs1 ${BINDIR}
	-rm -f ${BINDIR}/sha1
	${INSTALL} -c -m 0555 sha1 ${BINDIR}
	-rm -f ${BINDIR}/md5
	${INSTALL} -c -m 0555 md5 ${BINDIR}
	-@if [ ! -d "${TOPDIR}" ]; then \
	    echo "mkdir ${TOPDIR}"; \
	    mkdir "${TOPDIR}"; \
	fi
	-@if [ ! -d "${SHSLIB}" ]; then \
	    echo "mkdir ${SHSLIB}"; \
	    mkdir "${SHSLIB}"; \
	fi
	${INSTALL} -c -m 0444 file1 ${SHSLIB}
	${INSTALL} -c -m 0444 file2 ${SHSLIB}
	${INSTALL} -c -m 0444 shs1.data shs.data ${SHSLIB}
	${INSTALL} -c -m 0444 shs1.ddata shs.ddata ${SHSLIB}
	-@if [ ! -d "${MD5LIB}" ]; then \
	    echo "mkdir ${MD5LIB}"; \
	    mkdir "${MD5LIB}"; \
	fi
	${INSTALL} -c -m 0444 file1 ${MD5LIB}
	${INSTALL} -c -m 0444 file2 ${MD5LIB}
	${INSTALL} -c -m 0444 md5.data ${MD5LIB}
	-@if [ -z "${MANDIR}" ]; then \
	    echo "man pages are not installed, $${MANDIR} is empty"; \
	else \
	    echo "rm -f ${MANDIR}/md5.${MANEXT}"; \
	    rm -f ${MANDIR}/md5.${MANEXT}; \
	    echo "rm -f ${MANDIR}/shs.${MANEXT} ${MANDIR}/shs1.${MANEXT}"; \
	    rm -f ${MANDIR}/shs.${MANEXT} ${MANDIR}/shs1.${MANEXT}; \
	    echo "rm -f ${MANDIR}/sha.${MANEXT} ${MANDIR}/sha1.${MANEXT}"; \
	    rm -f ${MANDIR}/sha.${MANEXT} ${MANDIR}/sha1.${MANEXT}; \
	    echo "cp shs.1 ${MANDIR}/shs.${MANEXT}"; \
	    cp shs.1 ${MANDIR}/shs.${MANEXT}; \
	    echo "cp shs.1 ${MANDIR}/shs1.${MANEXT}"; \
	    cp shs.1 ${MANDIR}/shs1.${MANEXT}; \
	    echo "cp shs.1 ${MANDIR}/sha.${MANEXT}"; \
	    cp shs.1 ${MANDIR}/sha.${MANEXT}; \
	    echo "cp shs.1 ${MANDIR}/sha1.${MANEXT}"; \
	    cp shs.1 ${MANDIR}/sha1.${MANEXT}; \
	    echo "cp md5.1 ${MANDIR}/md5.${MANEXT}"; \
	    cp md5.1 ${MANDIR}/md5.${MANEXT}; \
	    echo "chmod 0444 ${MANDIR}/md5.${MANEXT}"; \
	    chmod 0444 ${MANDIR}/md5.${MANEXT}; \
	    echo "chmod 0444 ${MANDIR}/shs.${MANEXT}"; \
	    chmod 0444 ${MANDIR}/shs.${MANEXT}; \
	    echo "chmod 0444 ${MANDIR}/shs1.${MANEXT}"; \
	    chmod 0444 ${MANDIR}/shs1.${MANEXT}; \
	    echo "chmod 0444 ${MANDIR}/sha.${MANEXT}"; \
	    chmod 0444 ${MANDIR}/sha.${MANEXT}; \
	    echo "chmod 0444 ${MANDIR}/sha1.${MANEXT}"; \
	    chmod 0444 ${MANDIR}/sha1.${MANEXT}; \
	fi

clean:
	rm -f shs1.o shs1drvr.o shs1dual.o mega mega.o
	rm -f shs.o shsdrvr.o shsdual.o
	rm -f md5.o md5drvr.o md5dual.o
	rm -f shsio.o shs1io.o md5io.o
	rm -f syscrypt.o
	rm -f file1 file3
	rm -f endian.h endian.o endian
	rm -f align.h must_align must_align.o
	rm -f core core.must_align core.shs core.md5

clobber: clean
	rm -f shs md5 shs1 sha sha1
	rm -rf endian.dSYM mega.dSYM
