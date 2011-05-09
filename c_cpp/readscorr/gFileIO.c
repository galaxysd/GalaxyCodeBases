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

ssize_t read_kseq_with2bit(SeqFileObj * const seqObj) {
    size_t seqlen;	// in fact size_t, but minus values are meanful.
    int_fast8_t rvalue = kseq_read(seqObj->fh);
    if (rvalue>0) {
        uint_fast8_t type = rvalue; // 1 or 3. No need to &3 now.
        seqlen = ((kseq_t*)(seqObj->fh))->seq.l;
        if (rvalue&2) { // withQ
            //encodeQ;
            type |= 8;
        } 
        size_t needtomallocDW = (seqlen+3)>>2;  // 1 DWord = 4 Bytes. Well, I just mean the 4-time relationship.
        if (needtomallocDW > seqObj->binMallocedDWord) {
            KROUNDUP32(needtomallocDW);
            seqObj->binMallocedDWord = needtomallocDW;
            seqObj->diBseq = realloc(seqObj->diBseq,needtomallocDW);
            seqObj->hexBQ = realloc(seqObj->hexBQ,needtomallocDW<<2);
        }
        seqObj->readlength = seqlen;
        seqObj->type = type;
        return seqlen;
    } else return rvalue;
}

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
		//seqObj->readlength = &seq->seq.l;
		seqObj->fh = seq;
		seqObj->getNextSeq = read_kseq_with2bit;	// (int (*)(void*))
		seqObj->diBseq = NULL;	// We need NULL to free ...
		seqObj->hexBQ = NULL;
		seqObj->binMallocedDWord = 0;
		//seqObj->readlength = 0;
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
/*
ssize_t inSeqFreadNext(SeqFileObj * const seqObj) {
	ssize_t seqlen;	// in fact size_t, but minus values are meanful.
	seqlen = (*seqObj->getNextSeq)(seqObj);
	//seqObj->hasQ = (*seqObj->qual != 0);
	return seqlen;
}
*/
void inSeqFdestroy(SeqFileObj * const seqObj) {
	kseq_destroy(seqObj->fh);
	free(seqObj->diBseq);
	free(seqObj->hexBQ);
	free(seqObj);
}

