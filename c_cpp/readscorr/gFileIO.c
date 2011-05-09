#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>	// open, above 3
#include <unistd.h>	// read, close
#include <err.h>
#include <string.h>	// memcmp
#include <stdlib.h>	//calloc
#include <stdio.h>	// puts
#include <stdint.h>	// uint_fast8_t
//#include <zlib.h>	// already in gkseq.h
#include "gkseq.h"	// No more `kseq.h` with tons of Macros !!!
#include "gFileIO.h"
//KSEQ_INIT(gzFile, gzread)	// [kseq.h] Just like include, to inline some static inline functions.

// binmode: 1=base char, 2=base 2bit
SeqFileObj * inSeqFinit(const char * const filename, unsigned char binmode) {
	int fd;
	ssize_t fdstat;
	char FileID[G_HEADER_LENGTH];
	fd = open(filename, O_RDONLY);
	if (fd == -1) {warn("\n[x]Cannot open Sequence_file [%s]", filename); return NULL;}
	fdstat = read(fd, FileID, G_HEADER_LENGTH);
	if ( fdstat != G_HEADER_LENGTH ) {
		warn("[!]Too short (%zd) for sequence_file [%s].", fdstat, filename);
		return NULL;
	}
	fdstat = close(fd);
	if ( fdstat )
		warn("[!]ErrNo:[%zd] while closing file [%s].", fdstat, filename);

	SeqFileObj * const seqObj = malloc(sizeof(SeqFileObj));	// no more calloc(), just {.datePos[1]=0} is needed for kseq way.
	//CANNOT use `&(SeqFileObj) {.datePos={0}};`, which is of a short lifetime !

	if ( memcmp(G_TYPE_FAQC, FileID, G_HEADER_LENGTH) ) {	// Not G_TYPE_FAQC, thus to kseq.h
		gzFile fp;
		kseq_t *seq;

		fp = gzopen(filename, "r");
		if (! fp) err(EXIT_FAILURE, "\n[x]Cannot open with zlib for [%s]", filename);
		seq = kseq_init(fp);	// calloc, thus safe
		seqObj->datePos[0] = -1;	// seeking not available here.
		seqObj->datePos[1] = 0;
		//seqObj->hasQ = 0;
		seqObj->name = &seq->name.s;
		seqObj->comment = &seq->comment.s;
		seqObj->seq = &seq->seq.s;
		seqObj->qual = &seq->qual.s;
		seqObj->readlength = &seq->seq.l;
		seqObj->fh = seq;
		seqObj->getNextSeq = (G_int_oneIN) kseq_read;	// (int (*)(void*))
		seqObj->diBseq = NULL;	// later ...
		seqObj->hexBQ = NULL;
	//int seqlen;	// TEST ONLY !
	//seqlen = (*seqObj->getNextSeq)(seqObj->fh);	//seqlen = kseq_read(seq); // TEST ONLY !
		
	} else {	// is G_TYPE_FAQC
		errx(EXIT_FAILURE, "\n[:(]FAQC file not supported now.");
	}

	return seqObj;
}

/*
int_fast8_t read_kseq_no2bit()
int_fast8_t read_kseq_with2bit()
int_fast8_t read_faqc_nobasechar()
int_fast8_t read_faqc_withbasechar()
*/

ssize_t inSeqFreadNext(SeqFileObj * const seqObj) {
	ssize_t seqlen;	// in fact size_t, but minus values are meanful.
	seqlen = (*seqObj->getNextSeq)(seqObj->fh);
	//seqObj->hasQ = (*seqObj->qual != 0);
	return seqlen;
}

void inSeqFdestroy(SeqFileObj * const seqObj) {
	kseq_destroy(seqObj->fh);
	free(seqObj);
}

