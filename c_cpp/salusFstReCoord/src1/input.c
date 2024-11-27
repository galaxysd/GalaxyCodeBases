#include <stdio.h>  // fprintf
#include <zlib-ng.h>

#include "kseq.h"
KSEQ_INIT(gzFile, zng_gzread)
#include "common.h"

void fqReader_init(void) {
	int rc;
	const char *pattern = "R[0-9]{3}(:)?C[0-9]{3}";
	if ((rc = regcomp(&Parameters.regex, pattern, REG_EXTENDED)) != 0) {
		regerror(rc, &Parameters.regex, Parameters.buffer, PARAMETERS_BUFFER_SIZE);
		printf("[x]=%d=regcomp(\"%s\") failed with '%s'.\n", rc, pattern, Parameters.buffer);
		exit(EXIT_FAILURE);
	}
	Parameters.ksfp = zng_gzopen(Parameters.inFastqFilename, "r");
	if (unlikely(Parameters.ksfp == NULL)) {
		fprintf(stderr, "[x]gzopen error on opening [%s]: %s.\n", Parameters.inFastqFilename, strerror(errno));
		exit(1);
	}
	rc = zng_gzbuffer(Parameters.ksfp, GZBUFSIZE);
	if (unlikely(rc != 0)) {
		fprintf(stderr, "[x]gzbuffer error: %s.\n", strerror(errno));
		exit(1);
	}
	Parameters.kseq = kseq_init(Parameters.ksfp);
}

void fqReader_destroy(void) {
	regfree(&Parameters.regex);
	kseq_destroy(Parameters.kseq);
	zng_gzclose(Parameters.ksfp);
}

void fill_worker(int_least16_t worker_id) {
	kseq_t *seq = Parameters.kseq;
	workerArray_t *worker = &Parameters.workerArray[worker_id];
	uint64_t index = 0;
	regmatch_t matches[2];

	for (uint64_t index = 0; index < JOBITEMSIZE; index++) {
		fstBCdata_t *fstBCdata_p = &worker->input_array[index];
		if (kseq_read(seq) >= MINBARCODELEN) {
			// strncpy((char *)fstBCdata_p->name, seq->name.s, sizeof(fstBCdata_p->name));
			STRcpyARRAY(fstBCdata_p->name, seq->name.s);
			/* seq->comment.s is discarded */
			STRcpyARRAY(fstBCdata_p->seq, seq->seq.s);
			STRcpyARRAY(fstBCdata_p->qual, seq->qual.s);
#ifndef RELEASE
			fprintf(stderr, "- %llu -\n", index);
			ARRAYcpySTR(Parameters.buffer, fstBCdata_p->name);
			// snprintf(Parameters.buffer, 1 + sizeof(fstBCdata_p->name), "%s", fstBCdata_p->name);
			fprintf(stderr, "->Name:[%s]\n", Parameters.buffer);
			ARRAYcpySTR(Parameters.buffer, fstBCdata_p->seq);
			fprintf(stderr, "->Sequ:[%s]\n", Parameters.buffer);
			ARRAYcpySTR(Parameters.buffer, fstBCdata_p->qual);
			fprintf(stderr, "->Qual:[%s]\n", Parameters.buffer);
#endif
		} else {
			fstBCdata_p->name[0] = '\0';
			fstBCdata_p->seq[0] = '\0';
			fstBCdata_p->qual[0] = '\0';
		}
		// continue;
	}
	// Update the worker's flag to indicate the processing is complete
	atomic_store(&worker->flag, 1);
}
