#include <stdio.h>   // printf, fopen
#include <stdlib.h>  // calloc

#include "common.h"

const char *argp_program_version = "fstBC transCorrd demo 0.1 @" __TIME__ "," __DATE__;
const char *argp_program_bug_address = "huxs@salus-bio.com";

/* Program documentation. */
static char doc[] =
    "fstBC transCorrd multi-threaded"
#if defined(DEBUG)
    " Debug Version"
#elif defined(NDEBUG)
    " Release Version"
#else
    " with assert"
#endif
#ifdef USE_ZLIBNG
    " (zlib-ng)"
#elif defined(USE_ZLIB)
    " (zlib)"
#else /* USE_LIBISAL */
    " (ISA-L)"
#endif
    ;

/* Global Var for "common.h" */
parameters_t Parameters = {
    .fqSkipped = 0,
    .fqValid = 0,
};
// static_assert(VARTYPE(Parameters.jobDataState)==1, "It is not uint8_t");

int main(int argc, char *argv[]) {
	uv_thread_t main_tid = uv_thread_self();
	errno = 0; /* extern int, but do not declare errno manually */
	if (argc < 2) {
		// https://github.com/facebook/zstd/issues/1155#issuecomment-400052833 Both rsync block-size and zfs recordsize are 128K max.
		fprintf(stderr, "====== %s ======\n[i]Usage: %s <fstBC.fq[.gz]> [unZoomRatio=1] 2>run.log | zstd -B131072 --rsyncable --adapt -T12 -fo spatialBC.fq.zstd \n", doc, argv[0]);
		return 2;
	} else {
		fprintf(stderr, "[!]Arguments[%d]:", argc);
		for (int i = 0; i < argc; i++) {
			fprintf(stderr, " [%s]", argv[i]);
		}
		fputs(".\n", stderr);
	}
	Parameters.inFastqFilename = argv[1];
	fqReader_init();
	if (argc >= 3) {
		Parameters.unZoomRatio = strtof(argv[2], NULL);
		if (errno != 0) {
			fprintf(stderr, "[x]Error on strtof(%s) to [%f]: %s.\n", argv[2], Parameters.unZoomRatio, strerror(errno));
			exit(EXIT_FAILURE);
		}
		if (Parameters.unZoomRatio < 0) {
			Parameters.unZoomRatio = -Parameters.unZoomRatio;
		} else if (Parameters.unZoomRatio == 0) {
			Parameters.unZoomRatio = 1.0f;
		} else if (Parameters.unZoomRatio > 1000) {
			Parameters.unZoomRatio = 1.0f;
		}
	} else {
		Parameters.unZoomRatio = 1.0f;
	}

#ifndef RELEASE
#if __has_builtin(__builtin_dump_struct)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-pedantic"
	__builtin_dump_struct(&Parameters, &fprintf, stderr);
#pragma GCC diagnostic pop
#endif
#endif

	int2hex(Parameters.buffer, &main_tid, sizeof(main_tid));
	fprintf(stderr, "[!]P:[%d]%sh, Available_Cores:[%u].\n", uv_os_getpid(), Parameters.buffer, uv_available_parallelism());
	// fprintf(stderr, "[!]P:[%d]0x%.*" PRIx64 ".\n", uv_os_getpid(), (int)2*sizeof(main_tid), (uint64_t) main_tid);
	// fprintf(stderr, "[!]P:[%d]%" PRIx64 "h, Available_Cores:[%u].\n", uv_os_getpid(), (uint64_t)main_tid, uv_available_parallelism());
	uv_loop_init(&Parameters.loop);

	uv_loop_close(&Parameters.loop);
	// defLoop_p = uv_default_loop();
	while (Parameters.ksflag > 0) {  // for multi-threads demo, we use workQueue 1 of [0,JOBQUEUESIZE-1].
		fill_worker(1);
		worker(1);
		output_worker(1);
	}
	fprintf(stderr, "done: %s.\n", strerror(errno));
	fqReader_destroy();
	fprintf(stderr, "[i]Run to the end.\n");
	return 0;
}

/*
seqkit fx2tab t.fq |awk -F'[\t \037]' -v OFS='\t' '{print $5,$6,$4,$2" "$3}'|awk -F'\t' '{print $1,$NF}'
*/
