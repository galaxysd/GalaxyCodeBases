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

KSEQ_INIT(gzFile, gzread)

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

	//SeqFileObj * const seqObj = NULL;
	SeqFileObj * const seqObj = calloc(1,sizeof(SeqFileObj));	//&(SeqFileObj) {.datePos={0}}; thread safe sine not static

	if ( memcmp(G_TYPE_FAQC, FileID, G_HEADER_LENGTH) ) {	// Not G_TYPE_FAQC, thus to kseq.h
		gzFile fp;
		kseq_t *seq;

		int seqlen;
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
	(*seqObj->getNextSeq)(seqObj->fh);	//seqlen = kseq_read(seq);	// TEST ONLY !
		
	} else {	// is G_TYPE_FAQC
		errx(EXIT_FAILURE, "FAQC file not supported now.");
	}

	return seqObj;
}

/*
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l;
	if (argc == 1) {
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		return 1;
	}
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		printf("name: %s\n", seq->name.s);
		if (seq->comment.l) printf("comment: %s\n", 
seq->comment.s);
		printf("seq: %s\n", seq->seq.s);
		if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
	}
	printf("return value: %d\n", l);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
*/

