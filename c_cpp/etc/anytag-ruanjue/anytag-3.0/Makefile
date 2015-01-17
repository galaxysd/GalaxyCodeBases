VERSION=3.0
MINOR_VER=20121023
CC=gcc
CFLAGS=-g -W -Wall -O2 -finline-functions -D_FILE_OFFSET_BITS=64
GLIBS=-lm -lpthread
GENERIC_SRC=Makefile string.h vector.h file_reader.h file_reader.c bitvec.h bitsvec.h hashset.h sort.h list.h heap.h dna.h stdaln.h stdaln.c thread.h timer.h aln_cigar.h sr_aln.h sr_aln.c

ifeq (0, ${MAKELEVEL})
UNIQ_ID=$(shell date +"%s")
endif
all: anytag-$(VERSION) fillin


anytag-$(VERSION): $(GENERIC_SRC) ps_assembly.h ps_assembly.c anytag_aux.h anytag.c main.c 
	$(CC) $(CFLAGS) $(GLIBS) -DVERSION="\"$(VERSION) $(MINOR_VER)\"" -DUID="\"$(UNIQ_ID)\"" -o anytag-$(VERSION) file_reader.c stdaln.c sr_aln.c ps_assembly.c anytag.c main.c

fillin: $(GENERIC_SRC) ps_assembly.h ps_assembly.c fillin.c
	$(CC) $(CFLAGS) $(GLIBS) -o fillin file_reader.c stdaln.c sr_aln.c ps_assembly.c fillin.c

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out anytag-$(VERSION) fillin

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out
