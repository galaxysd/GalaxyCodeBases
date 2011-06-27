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
#include "2bitseq.h"
#include "chrseq.h"

ssize_t read_kseq_no2bit(SeqFileObj * const seqObj) {
    size_t seqlen;	// in fact size_t, but minus values are meanful.
    int_fast8_t rvalue = kseq_read(seqObj->fobj);
    if (rvalue>0) {
        uint_fast8_t type = rvalue; // 1 or 3. No need to &3 now.
        kseq_t *kseq;
        kseq = seqObj->fobj;
        seqlen = kseq->seq.l;
		seqObj->name = kseq->name.s;
		seqObj->comment = kseq->comment.s;
		seqObj->seq = kseq->seq.s;
		NormalizeChrSeq(kseq->seq.s);   // to /[ATCGN]*/
        if (rvalue&2) { // withQ
            //encodeQ;
            type |= 8u;
            seqObj->qual = kseq->qual.s;
        } else {
            seqObj->qual = NULL;
        }
        seqObj->readlength = seqlen;
        seqObj->type = type;
        return seqlen;
    } else return rvalue;
}
ssize_t read_kseq_with2bit(SeqFileObj * const seqObj) {
    size_t seqlen;	// in fact size_t, but minus values are meanful.
    int_fast8_t rvalue = kseq_read(seqObj->fobj);
    if (rvalue>0) {
        uint_fast8_t type = rvalue; // 1 or 3. No need to &3 now.
        kseq_t *kseq;
        kseq = seqObj->fobj;
        seqlen = kseq->seq.l;
		seqObj->name = kseq->name.s;
		seqObj->comment = kseq->comment.s;
		seqObj->seq = kseq->seq.s;
        if (rvalue&2) { // withQ
            //encodeQ;
            type |= 8u;
            seqObj->qual = kseq->qual.s;
        } else {
            seqObj->qual = NULL;
        }
        size_t needtomallocQQW = (seqlen+31u)>>5;  // 1 "QQWord" = 4 QWord = 32 bp. Well, I say there is QQW.
        if (needtomallocQQW > seqObj->binMallocedQQWord) {
            KROUNDUP32(needtomallocQQW);
            seqObj->binMallocedQQWord = needtomallocQQW;
            seqObj->diBseq = realloc(seqObj->diBseq,needtomallocQQW<<3);	// 2^3=8
            seqObj->hexBQ = realloc(seqObj->hexBQ,needtomallocQQW<<5);	// 4*2^3=32
        }
        seqObj->binNcount = base2dbit(seqlen, kseq->seq.s, seqObj->qual, seqObj->diBseq, seqObj->hexBQ);
// printf("-[%s]<%s><%zx>[%s]-\n",kseq->seq.s, qstr, seqObj->diBseq[0], unit2basechr(seqObj->diBseq[0]));
    // Well, how to deal with smallcase masking ? Not using this information yet.
        NormalizeChrSeq(kseq->seq.s);   // to /[ATCGN]*/
        seqObj->readlength = seqlen;
        seqObj->type = type;
        return seqlen;
    } else return rvalue;
}
void close_kseq(SeqFileObj * const seqObj){
    gzFile fp = ((kseq_t*)seqObj->fobj)->f->f;
    gzclose(fp);
    kseq_destroy(seqObj->fobj);
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
		seqObj->fobj = seq;
		if (binmode & GFIODIBBASE) {
		    seqObj->getNextSeq = read_kseq_with2bit;	// (ssize_t (*)(void*))
		} else {
		    seqObj->getNextSeq = read_kseq_no2bit;
		}
		seqObj->closefh = close_kseq;
		seqObj->diBseq = NULL;	// We need NULL to free ...
		seqObj->hexBQ = NULL;
		seqObj->binMallocedQQWord = 0;
		seqObj->binNcount = 0;
	} else {	// is G_TYPE_FAQC
		errx(EXIT_FAILURE, "\n[:(]FAQC file not supported now.");
	}

	return seqObj;
}

void inSeqFdestroy(SeqFileObj * const seqObj) {
	(*seqObj->closefh)(seqObj);
	free(seqObj->diBseq);
	free(seqObj->hexBQ);
	free(seqObj);
}

