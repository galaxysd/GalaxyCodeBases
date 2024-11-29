#include <stdio.h>  // fprintf

#include "kseq.h"
#if ZLIB_ID == 1
#include <zlib-ng.h>
KSEQ_INIT(gzFile, zng_gzread)
#define _GZ_OPEN_ zng_gzopen
#define _GZ_BUFFER_ zng_gzbuffer
#define _GZ_CLOSE_ zng_gzclose
#else
#if ZLIB_ID == 2
#include "izlib.h"
#else
#include <zlib.h>
#define _GZ_BUFFER_ gzbuffer
#endif
#define _GZ_OPEN_ gzopen
#define _GZ_CLOSE_ gzclose
KSEQ_INIT(gzFile, gzread)
#endif

#include "common.h"

void fqReader_init(void) {
	int rc;
	const char *pattern = "R[0-9]{3}(:)?C[0-9]{3}";
	if ((rc = regcomp(&Parameters.regex, pattern, REG_EXTENDED)) != 0) {
		regerror(rc, &Parameters.regex, Parameters.buffer, PARAMETERS_BUFFER_SIZE);
		printf("[x]=%d=regcomp(\"%s\") failed with '%s'.\n", rc, pattern, Parameters.buffer);
		exit(EXIT_FAILURE);
	}
	Parameters.ksfp = _GZ_OPEN_(Parameters.inFastqFilename, "r");
	if (unlikely(Parameters.ksfp == NULL)) {
		fprintf(stderr, "[x]gzopen error on opening [%s]: %s.\n", Parameters.inFastqFilename, strerror(errno));
		exit(1);
	}
#if defined(_GZ_BUFFER_)
	rc = _GZ_BUFFER_(Parameters.ksfp, GZBUFSIZE);
	if (unlikely(rc != 0)) {
		fprintf(stderr, "[x]gzbuffer error: %s.\n", strerror(errno));
		exit(1);
	}
#endif
	Parameters.kseq = kseq_init(Parameters.ksfp);
	// Parameters.worksQuene = (workerArray_t *)calloc(JOBQUEUESIZE, sizeof(workerArray_t));
	Parameters.worksQuene = (workerArray_t *)calloc(4, sizeof(workerArray_t));
}

void fqReader_destroy(void) {
	regfree(&Parameters.regex);
	kseq_destroy(Parameters.kseq);
	_GZ_CLOSE_(Parameters.ksfp);
	free(Parameters.worksQuene);
}

void fill_worker(int_least16_t worker_id) {
	kseq_t *seq = Parameters.kseq;
	workerArray_t *worker = &Parameters.worksQuene[worker_id];
	uint64_t index = 0;
	for (index = 0; index < JOBITEMSIZE; index++) {
		fstBCdata_t *fstBCdata_p = &worker->jobDatArray[index];
		if (kseq_read(seq) >= MINBARCODELEN) {
			// strncpy((char *)fstBCdata_p->name, seq->name.s, sizeof(fstBCdata_p->name));
			STRcpyARRAY(fstBCdata_p->name, seq->name.s);
			/* seq->comment.s is discarded */
			STRcpyARRAY(fstBCdata_p->seq, seq->seq.s);
			STRcpyARRAY(fstBCdata_p->qual, seq->qual.s);
#ifndef RELEASE
			if (unlikely(seq->comment.l > 0)) {
				size_t oldSize = MALLOCSIZE(fstBCdata_p->comment);
				if (unlikely(1 + seq->comment.l > oldSize)) {
					free(fstBCdata_p->comment);
					fstBCdata_p->comment = strdup(seq->comment.s);
				} else {
					STRcpySTR(fstBCdata_p->comment, seq->comment.s);
				}
			} else {
				if (unlikely(fstBCdata_p->comment != NULL)) {
					free(fstBCdata_p->comment);
					fstBCdata_p->comment = NULL;
				}
			}
#endif
#ifdef DEBUG
			fprintf(stderr, "- %llu -\n", index);
			ARRAYcpySTR(Parameters.buffer, fstBCdata_p->name);
			// snprintf(Parameters.buffer, 1 + sizeof(fstBCdata_p->name), "%s", fstBCdata_p->name);
			fprintf(stderr, "->Name:[%s] Comment:[%s]\n", Parameters.buffer, fstBCdata_p->comment);
			ARRAYcpySTR(Parameters.buffer, fstBCdata_p->seq);
			fprintf(stderr, "->Sequ:[%s]\n", Parameters.buffer);
			ARRAYcpySTR(Parameters.buffer, fstBCdata_p->qual);
			fprintf(stderr, "->Qual:[%s]\n", Parameters.buffer);
#endif
		} else {
			fstBCdata_p->name[0] = '\0';
			fstBCdata_p->seq[0] = '\0';
			fstBCdata_p->qual[0] = '\0';
#ifndef RELEASE
			free(fstBCdata_p->comment);
			fstBCdata_p->comment = NULL;
#endif
		}
		// continue;
	}
	// Update the worker's flag to indicate the processing is complete
	atomic_store(&worker->flag, 1);
}
