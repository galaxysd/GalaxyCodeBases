# Makefile for ssaha_pileup 
CC= gcc
CFLAGS= -O2
LFLAGS= -lm

SSAHAPILEUPO= fast.o ssaha_pileup.o 
SSAHAPAIRSO=ssaha_pairs.o 
SSAHASOLEXA=ssaha_solexa.o 
SSAHACLEANO=ssaha_clean.o 
SSAHAMATESO=ssaha_mates.o 
SSAHACIGARO=ssaha_cigar.o 
SSAHAINDELO= fast.o ssaha_indel.o 
SSAHAVIEWER= fast.o view_pileup.o 
SSAHAREADSO= fast.o ssaha_reads.o 
SSAHAMERGEO= fast.o ssaha_merge.o 
SSAHACHECKC= fast.o ssaha_check-cigar.o 

all  : ssaha_pileup ssaha_solexa ssaha_pairs ssaha_clean ssaha_mates ssaha_cigar ssaha_indel view_pileup ssaha_reads ssaha_merge ssaha_check-cigar

ssaha_pileup: makefile $(SSAHAPILEUPO)
	$(CC) $(CFLAGS) -o $@ $(SSAHAPILEUPO) $(LFLAGS) 
	chmod o-r ssaha_pileup 

ssaha_pairs: makefile $(SSAHAPAIRSO)
	$(CC) $(CFLAGS) -o $@ $(SSAHAPAIRSO) $(LFLAGS) 
	chmod o-r ssaha_pairs 

ssaha_clean: makefile $(SSAHACLEANO)
	$(CC) $(CFLAGS) -o $@ $(SSAHACLEANO) $(LFLAGS) 
	chmod o-r ssaha_clean 

ssaha_solexa: makefile $(SSAHASOLEXA)
	$(CC) $(CFLAGS) -o $@ $(SSAHASOLEXA) $(LFLAGS) 
	chmod o-r ssaha_solexa 

ssaha_mates: makefile $(SSAHAMATESO)
	$(CC) $(CFLAGS) -o $@ $(SSAHAMATESO) $(LFLAGS) 
	chmod o-r ssaha_mates 

ssaha_cigar: makefile $(SSAHACIGARO)
	$(CC) $(CFLAGS) -o $@ $(SSAHACIGARO) $(LFLAGS) 
	chmod o-r ssaha_cigar 

ssaha_check-cigar: makefile $(SSAHACHECKC)
	$(CC) $(CFLAGS) -o $@ $(SSAHACHECKC) $(LFLAGS) 
	chmod o-r ssaha_check-cigar 

ssaha_indel: makefile $(SSAHAINDELO)
	$(CC) $(CFLAGS) -o $@ $(SSAHAINDELO) $(LFLAGS) 
	chmod o-r ssaha_indel 

view_pileup: makefile $(SSAHAVIEWER)
	$(CC) $(CFLAGS) -o $@ $(SSAHAVIEWER) $(LFLAGS) 
	chmod o-r view_pileup 

ssaha_reads: makefile $(SSAHAREADSO)
	$(CC) $(CFLAGS) -o $@ $(SSAHAREADSO) $(LFLAGS) 
	chmod o-r ssaha_reads 

ssaha_merge: makefile $(SSAHAMERGEO)
	$(CC) $(CFLAGS) -o $@ $(SSAHAMERGEO) $(LFLAGS) 
	chmod o-r ssaha_merge 


