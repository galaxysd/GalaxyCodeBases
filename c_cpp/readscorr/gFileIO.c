#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>	// open, above 3
#include <unistd.h>	// read, close
#include <err.h>
#include <string.h>	// memcmp
#include <stdlib.h>	//calloc
#include <stdio.h>	// puts
#include <zlib.h>
#include "kseq.h"
#include "gFileIO.h"
KSEQ_INIT(gzFile, gzread)	// Just like include, to inline some static inline functions.

SeqFileObj * inSeqFinit(const char * const filename) {
	int fd;
	ssize_t fdstat;
	char FileID[G_HEADER_LENGTH];
	fd = open(filename, O_RDONLY);
	if (fd == -1) err(EXIT_FAILURE, "Cannot open Sequence_file [%s]", filename);
	fdstat = read(fd, FileID, G_HEADER_LENGTH);
	if ( fdstat != G_HEADER_LENGTH )
		err(EXIT_FAILURE, "Too short (%zd) for sequence_file [%s].", fdstat, filename);
	fdstat = close(fd);
	if ( fdstat )
		warn("ErrNo:[%zd] while closing file [%s].", fdstat, filename);

	SeqFileObj * const seqObj = calloc(1,sizeof(SeqFileObj));	// in fact, just {.datePos[1]=0} is needed for kseq way.
	//CANNOT use `&(SeqFileObj) {.datePos={0}};`, which is of a short lifetime !

	if ( memcmp(G_TYPE_FAQC, FileID, G_HEADER_LENGTH) ) {	// Not G_TYPE_FAQC, thus to kseq.h
		gzFile fp;
		kseq_t *seq;

		fp = gzopen(filename, "r");
		if (! fp) err(EXIT_FAILURE, "Cannot open with zlib for [%s]", filename);
		seq = kseq_init(fp);	// calloc, thus safe
		seqObj->datePos[0] = -1;	// seeking not available here.
		seqObj->name = &seq->name;
		seqObj->comment = &seq->comment;
		seqObj->seq = &seq->seq;
		seqObj->qual = &seq->qual;
		seqObj->fh = seq;
		seqObj->getNextSeq = (G_int_oneIN) kseq_read;	// (int (*)(void*))
		seqObj->diBseq = NULL;	// later ...
		seqObj->hexBQ = NULL;
	int seqlen;	// TEST ONLY !
	seqlen = (*seqObj->getNextSeq)(seqObj->fh);	//seqlen = kseq_read(seq); // TEST ONLY !
		
	} else {	// is G_TYPE_FAQC
		errx(EXIT_FAILURE, "FAQC file not supported now.");
	}

	return seqObj;
}
