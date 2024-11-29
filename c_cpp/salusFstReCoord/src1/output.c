#include "common.h"

void output_worker(int_least16_t worker_id) {
	workerArray_t *worker = &Parameters.worksQuene[worker_id];
	regmatch_t matches[2];
	char **splitSets = worker->tokens;
	char readName[MAXFQIDLEN + 1] = {0};
	char readSeq[BARCODELEN + 1] = {0};
}
