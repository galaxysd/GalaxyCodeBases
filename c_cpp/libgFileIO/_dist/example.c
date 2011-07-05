#define _GNU_SOURCE
#include <stdlib.h> //EXIT_FAILURE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <err.h>
#include "libgFileIO.h"

int main (void) {
	char *filename="test.seq.gz";

	ssize_t readlength;
	fprintf(stderr, "<%s>:\n", filename);
	SeqFileObj *seqobj = inSeqFinit(filename,GFIOCHRBASE);
	if (seqobj) {
		while ( (readlength = (*seqobj->getNextSeq)(seqobj)) >= 0 ) {
			printf("-ID:[%s,%s] Seqlen:%zu %zu ReturnValue:%zd\nSeq:[%s]\n  Q:[%s]\n*%zx,%u\n",
			seqobj->name,seqobj->comment,seqobj->readlength,seqobj->binMallocedQQWord,readlength,
			seqobj->seq,seqobj->qual,(size_t)seqobj->seq,seqobj->type);
		}
	} 
	fputs("\b\b\b\b, done !\n", stderr);
	inSeqFdestroy(seqobj);

	exit(EXIT_SUCCESS);
}

