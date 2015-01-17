VERSION=2.5.2
MINOR_VER=20120618
CC=gcc
CFLAGS=-g -W -Wall -O2 -finline-functions -D_FILE_OFFSET_BITS=64
GLIBS=-lm -lpthread
GENERIC_SRC=Makefile string.h vector.h file_reader.h file_reader.c bitvec.h hashset.h sort.h list.h heap.h dna.h stdaln.h stdaln.c thread.h timer.h aln_cigar.h sr_aln.h sr_aln.c

ifeq (0, ${MAKELEVEL})
UNIQ_ID=$(shell date +"%s")
endif
all: anytag-$(VERSION) sralign


anytag-$(VERSION): $(GENERIC_SRC) local_assembly.h local_assembly.c anytag_aux.h anytag_all.c main.c 
	$(CC) $(CFLAGS) $(GLIBS) -DVERSION="\"$(VERSION) $(MINOR_VER)\"" -DUID="\"$(UNIQ_ID)\"" -o anytag-$(VERSION) file_reader.c stdaln.c sr_aln.c local_assembly.c anytag_all.c main.c

sralign: $(GENERIC_SRC) sralign.c error_correction.h
	$(CC) $(CFLAGS) $(GLIBS) -o sralign sralign.c sr_aln.c file_reader.c stdaln.c

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out anytag-$(VERSION) sralign

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out
