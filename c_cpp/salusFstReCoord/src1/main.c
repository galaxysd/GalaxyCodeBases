#include <stdio.h>  // printf, fopen

#include "common.h"

#ifdef TARGET_OS_OSX
#include <malloc/malloc.h>
#endif

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
workerArray_t *workerArray;

int main(int argc, char *argv[]) {
	if (argc != 3) {
		fprintf(stderr, "====== %s ======\n[i]Usage: %s <fstBC.fq[.gz]> <spatial.txt>\n", doc, argv[0]);
		return 2;
	}
	printf("[i]%s %s %s\n", argv[0], argv[1], argv[2]);
	Parameters.inFastqFilename = argv[1];
	// 没printf拖时间就得加 barrier
	Parameters.outSpatialFilename = argv[2];
	errno = 0; /* extern int, but do not declare errno manually */
	Parameters.outfp = fopen(Parameters.outSpatialFilename, "w");
	if (errno != 0) {
		fprintf(stderr, "[x]Error on opening output file [%s]: %s.\n", Parameters.outSpatialFilename, strerror(errno));
		exit(EXIT_FAILURE);
	}
	// defLoop_p = uv_default_loop();

	int_least16_t n_threads = 1;
	workerArray_t *workerArray = (workerArray_t *)calloc(n_threads, sizeof(workerArray_t));

#ifndef RELEASE
#if __has_builtin(__builtin_dump_struct)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-pedantic"
	__builtin_dump_struct(&Parameters, &printf);
#pragma GCC diagnostic pop
#endif
#endif

	free(workerArray);
	fclose(Parameters.outfp);
	fprintf(stderr, "[i]Run to the end.\n");
	return 0;
}
