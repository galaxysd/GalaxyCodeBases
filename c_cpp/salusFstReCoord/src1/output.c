#include "common.h"

void output_worker(int_least16_t worker_id) {
	workerArray_t *worker = &Parameters.worksQuene[worker_id];
	regmatch_t matches[2];
	char **splitSets = worker->tokens;
#ifndef RELEASE
	char readName[MAXFQIDLEN + 1] = {0};
	char readSeq[BARCODELEN + 1] = {0};
	char readQual[BARCODELEN + 1] = {0};
	char readRowCol[ROWCOLSIZE + 1] = {0};
#endif
	for (uint64_t index = 0; index < JOBITEMSIZE; index++) {
		fstBCdata_t *fstBCdata_p = &worker->jobDatArray[index];
		if (fstBCdata_p->name[0] == 0) {
			break;
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
		ARRAYcpySTR(readRowCol, fstBCdata_p->RowCol);
		assert(readRowCol[sizeof(fstBCdata_p->RowCol)] == '\0');
		fprintf(stderr, "RoCo:[%s]\t", readRowCol);
		fprintf(stderr, "XY:[%.2f],[%.2f]\n\n", fstBCdata_p->newXY[0], fstBCdata_p->newXY[1]);
#endif
		fprintf(stdout, "@%.*s %.2f %.2f\t%.*s\n" "%.*s\n+\n%.*s\n",
			sizeof(fstBCdata_p->name), fstBCdata_p->name,
			fstBCdata_p->newXY[0], fstBCdata_p->newXY[1],
			sizeof(fstBCdata_p->RowCol), fstBCdata_p->RowCol,
			sizeof(fstBCdata_p->seq), fstBCdata_p->seq,
			sizeof(fstBCdata_p->qual), fstBCdata_p->qual
		);
		fflush(stdout);
	}
}
