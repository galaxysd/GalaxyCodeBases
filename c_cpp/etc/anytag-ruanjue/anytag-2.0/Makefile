CC=gcc
CFLAGS=-g -W -Wall -O2 -finline-functions -D_FILE_OFFSET_BITS=64
GLIBS=-lm -lpthread
GENERIC_SRC=Makefile string.h vector.h file_reader.h file_reader.c hashset.h sort.h list.h heap.h dna.h stdaln.h stdaln.c thread.h timer.h aln_cigar.h sr_aln.h sr_aln.c

all: anytag sra

anytag: $(GENERIC_SRC) local_assembly.h local_assembly.c anytag_aux.h anytag_all.c anytag_aln.c anytag_asm.c anytag_lnk.c anytag.c 
	$(CC) $(CFLAGS) $(GLIBS) -o anytag file_reader.c stdaln.c sr_aln.c local_assembly.c anytag_all.c anytag_aln.c anytag_asm.c anytag_lnk.c anytag.c

sra: $(GENERIC_SRC) sra.c
	$(CC) $(CFLAGS) $(GLIBS) -o sra file_reader.c stdaln.c sr_aln.c sra.c

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out anytag sra

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out
