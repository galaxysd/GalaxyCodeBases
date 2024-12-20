#include <stdio.h>  // fprintf

#ifdef __linux__
#include <sys/mman.h>  // madvise
#endif

#include "kseq.h"
#ifdef USE_ZLIBNG
#include <zlib-ng.h>
KSEQ_INIT(gzFile, zng_gzread)
#define _GZ_OPEN_ zng_gzopen
#define _GZ_BUFFER_ zng_gzbuffer
#define _GZ_CLOSE_ zng_gzclose
#else
#ifdef USE_ZLIB
#include <zlib.h>
#define _GZ_BUFFER_ gzbuffer
#else /* USE_LIBISAL */
#include "izlib.h"
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
		fprintf(stderr, "[x]=%d=regcomp(\"%s\") failed with '%s'.\n", rc, pattern, Parameters.buffer);
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
	Parameters.ksflag = 1;
	// Parameters.worksQuene = (workerArray_t *)calloc(JOBQUEUESIZE, sizeof(workerArray_t));
	size_t HPAGE_shift = 21;  // 1 << 21 = 2M
	size_t HPAGE_alignment = 1 << HPAGE_shift;
	size_t worksQuene_size = JOBQUEUESIZE * sizeof(workerArray_t);
	size_t worksQuene_aligned_cnt = worksQuene_size >> HPAGE_shift;
	if (worksQuene_size % HPAGE_alignment != 0) {
		worksQuene_aligned_cnt++;
	}
	size_t worksQuene_aligned_size = HPAGE_alignment * worksQuene_aligned_cnt;
	fprintf(stderr, "[!]aligned_alloc(%zu, %zu) for [%zu].\n", HPAGE_alignment, worksQuene_aligned_size, worksQuene_size);
	Parameters.worksQuene = (workerArray_t *)aligned_alloc(HPAGE_alignment, worksQuene_aligned_size);
	if (Parameters.worksQuene != NULL) {
#ifdef __linux__
		madvise(Parameters.worksQuene, worksQuene_aligned_size, MADV_HUGEPAGE);
#endif
		memset(Parameters.worksQuene, 0, worksQuene_aligned_size);
	} else {
		perror("[x]aligned_alloc failed.");
		exit(EXIT_FAILURE);
	}
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
	// if (atomic_load(&worker->flag) == 1) return;
	uint64_t index = 0;
	int kseq_ret = 0;
	for (index = 0; index < JOBITEMSIZE; index++) {
		fstBCdata_t *fstBCdata_p = &worker->jobDatArray[index];
		if ((kseq_ret = kseq_read(seq)) >= MINBARCODELEN) {
			// strncpy((char *)fstBCdata_p->name, seq->name.s, sizeof(fstBCdata_p->name));
			STRcpyARRAY(fstBCdata_p->name, seq->name.s);
			/* seq->comment.s is discarded */
			STRcpyARRAY(fstBCdata_p->seq, seq->seq.s);
			STRcpyARRAY(fstBCdata_p->qual, seq->qual.s);
#ifdef DEBUG
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
			if (kseq_ret < 0) {  // -1 for FEOF
				Parameters.ksflag = kseq_ret;
			}
		}
		// continue;
	}
	// Update the workerArray's flag to indicate it is filled, thus ready for worker.
	atomic_store(&worker->flag, 1);
}
