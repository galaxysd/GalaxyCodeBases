#include <stdio.h>   // printf, fopen
#include <stdlib.h>  // calloc

#include "common.h"

const char *argp_program_version = "fstBC transCorrd demo 0.1 @" __TIME__ "," __DATE__;
const char *argp_program_bug_address = "huxs@salus-bio.com";

/* Program documentation. */
static char doc[] =
    "fstBC transCorrd single threaded"
#if defined(DEBUG)
    " (Debug Version)"
#elif defined(NDEBUG)
    " (Release Version)"
#else
    " (with assert)"
#endif
    ;

/* Global Var for "common.h" */
parameters_t Parameters = {
    .fqSkipped = 0,
    .fqValid = 0,
};
// static_assert(VARTYPE(Parameters.jobDataState)==1, "It is not uint8_t");

int main(int argc, char *argv[]) {
	if (argc != 3) {
		fprintf(stderr, "====== %s ======\n[i]Usage: %s <fstBC.fq[.gz]> <spatial.txt>\n", doc, argv[0]);
		return 2;
	}
	printf("[i]%s %s %s\n", argv[0], argv[1], argv[2]);
	Parameters.inFastqFilename = argv[1];
	// 没printf拖时间就得加 barrier
	fqReader_init();
	Parameters.outSpatialFilename = argv[2];
	errno = 0; /* extern int, but do not declare errno manually */
	Parameters.outfp = fopen(Parameters.outSpatialFilename, "w");
	if (errno != 0) {
		fprintf(stderr, "[x]Error on opening output file [%s]: %s.\n", Parameters.outSpatialFilename, strerror(errno));
		exit(EXIT_FAILURE);
	}
	int_least16_t n_threads = 2;  // for multi-threads demo, we use threads 1 of [0,1].
	Parameters.workerArray = (workerArray_t *)calloc(n_threads, sizeof(workerArray_t));

	// defLoop_p = uv_default_loop();
	{
		fill_worker(1);
		worker(1);
	}
	fprintf(stderr, "done: %s.\n", strerror(errno));

#ifndef RELEASE
#if __has_builtin(__builtin_dump_struct)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-pedantic"
	__builtin_dump_struct(&Parameters, &printf);
#pragma GCC diagnostic pop
#endif
#endif

	free(Parameters.workerArray);
	fqReader_destroy();
	fclose(Parameters.outfp);
	fprintf(stderr, "[i]Run to the end.\n");
	return 0;
}
