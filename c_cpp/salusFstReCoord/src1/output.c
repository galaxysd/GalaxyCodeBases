#include "common.h"

void output_worker(int_least16_t worker_id) {
	workerArray_t *worker = &Parameters.worksQuene[worker_id];
	//if (atomic_load(&worker->flag) == 0) return;
#ifndef RELEASE
	char readName[MAXFQIDLEN + 1] = {0};
	char readSeq[BARCODELEN + 1] = {0};
	char readQual[BARCODELEN + 1] = {0};
#endif
	for (uint64_t index = 0; index < JOBITEMSIZE; index++) {
		fstBCdata_t *fstBCdata_p = &worker->jobDatArray[index];
		if (fstBCdata_p->name[0] == 0) {
			continue;
		}
#ifndef RELEASE
		ARRAYcpySTR(readName, fstBCdata_p->name);
		assert(readName[sizeof(fstBCdata_p->name)] == '\0');
		fprintf(stderr, "  ID:[%s]\n", readName);
		ARRAYcpySTR(readSeq, fstBCdata_p->seq);
		assert(readSeq[sizeof(fstBCdata_p->seq)] == '\0');
		fprintf(stderr, " Seq:[%s]\n", readSeq);
		ARRAYcpySTR(readQual, fstBCdata_p->qual);
		assert(readQual[sizeof(fstBCdata_p->qual)] == '\0');
		fprintf(stderr, "Qual:[%s]\n", readQual);
		fprintf(stderr, "R%03uC%03u XY:[%.2f],[%.2f]\n\n", fstBCdata_p->fov_row, fstBCdata_p->fov_column, fstBCdata_p->newXY[0], fstBCdata_p->newXY[1]);
#endif
		// clang-format off
		fprintf(stdout, "@%.*s %.2f %.2f" "%c" "R%03uC%03u\n" "%.*s\n+\n%.*s\n",
			(int)sizeof(fstBCdata_p->name), fstBCdata_p->name,
			fstBCdata_p->newXY[0], fstBCdata_p->newXY[1],
			_US_CHR_,
			fstBCdata_p->fov_row, fstBCdata_p->fov_column,
			(int)sizeof(fstBCdata_p->seq), fstBCdata_p->seq,
			(int)sizeof(fstBCdata_p->qual), fstBCdata_p->qual
		);
		// clang-format on
		fflush(stdout);
	}
	// Update the workerArray's flag to indicate it is empty for filling.
	atomic_store(&worker->flag, 0);
}
