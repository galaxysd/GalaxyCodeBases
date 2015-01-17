CC=gcc
CFLAGS=-g -W -Wall -O2 -finline-functions -D_FILE_OFFSET_BITS=64
GLIBS=-lm -lpthread
GENERIC_SRC=Makefile string.h vector.h file_reader.h file_reader.c hashset.h sort.h list.h heap.h dna.h stdaln.h stdaln.c thread.h timer.h

all: count_tag cluster_tag asm_tag sra anytag

count_tag: $(GENERIC_SRC) tag.h count_tag.c
	$(CC) $(CFLAGS) $(GLIBS) -o count_tag file_reader.c count_tag.c

cluster_tag: $(GENERIC_SRC) tag.h cluster_tag.c
	$(CC) $(CFLAGS) $(GLIBS) -o cluster_tag file_reader.c cluster_tag.c

asm_tag: $(GENERIC_SRC) aln_cigar.h sr_aln.h sr_aln.c simp_asm.h asm_tag.c
	$(CC) $(CFLAGS) $(GLIBS) -o asm_tag file_reader.c stdaln.c sr_aln.c asm_tag.c

sra: $(GENERIC_SRC) aln_cigar.h sr_aln.h sr_aln.c sra.c
	$(CC) $(CFLAGS) $(GLIBS) -o sra file_reader.c stdaln.c sr_aln.c sra.c

anytag: $(GENERIC_SRC) aln_cigar.h sr_aln.h sr_aln.c all_path.h all_path.c anytag.c 
	$(CC) $(CFLAGS) $(GLIBS) -o anytag file_reader.c stdaln.c sr_aln.c all_path.c anytag.c

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out count_tag cluster_tag asm_tag sra anytag

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out
