#ifndef FSTRECOORD_COMMON_H_
#define FSTRECOORD_COMMON_H_
#ifdef __cplusplus
extern "C" {
#endif

/* ################ */

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
// define something for Windows (32-bit and 64-bit, this part is common)
#ifdef _WIN64
// define something for Windows (64-bit only)
#else
// define something for Windows (32-bit only)
#endif
#elif __APPLE__
#include <TargetConditionals.h>
#if TARGET_OS_OSX
// OS X devices
#include <malloc/malloc.h>
#define MALLOCSIZE malloc_size
#elif TARGET_OS_MACCATALYST
// Mac's Catalyst (ports iOS API into Mac, like UIKit).
#elif TARGET_OS_DRIVERKIT
// macOS, iOS, Apple TV OS, or Apple Watch OS
#elif TARGET_OS_IPHONE
// firmware, devices, or simulator
#else
#error "Unknown Apple platform"
#endif
#elif __ANDROID__
// Below __linux__ check should be enough to handle Android,
// but something may be unique to Android.
#elif __linux__
// linux
#include <malloc.h>
#define MALLOCSIZE malloc_usable_size
#elif __unix__  // all unices not caught above
// Unix
#elif defined(_POSIX_VERSION)
// POSIX
#else
#error "Unknown compiler"
#endif

// https://meghprkh.github.io/blog/posts/c++-force-inline/
#if defined(__clang__)
#define FORCE_INLINE __attribute__((__always_inline__, __gnu_inline__)) extern inline

#elif defined(__GNUC__)
#define FORCE_INLINE __attribute__((__always_inline__)) inline

#elif defined(_MSC_VER)
#pragma warning(error : 4714)
#define FORCE_INLINE __forceinline

#else
#warning Unsupported compiler, fall back to `static inline`.
#define FORCE_INLINE static inline
#endif

#include <assert.h>
#include <ctype.h>      // in "kseq.h"
#include <errno.h>      // extern int errno
#include <regex.h>      // regex_t
#include <stdatomic.h>  // atomic_int_least8_t
#include <stddef.h>     // in "stdatomic.h"
#include <stdint.h>     // in "stdatomic.h"
#include <stdio.h>      // FILE, fflush
#include <stdlib.h>     // in "kseq.h", strtof
#include <string.h>     // in "kseq.h", memccpy, and so on.

#ifndef KSEQ_INIT
// #define STR_HELPER(x) #x
// #define STR(x) STR_HELPER(x)
// #pragma message "ZLIB_ID: " STR(ZLIB_ID)
#ifdef USE_ZLIBNG
#include <zlib-ng.h>
#elif defined(USE_ZLIB)
#include <zlib.h>
#else /* USE_LIBISAL */
#include "izlib.h"
#endif
#include "kseq.h"
KSEQ_DECLARE(gzFile)
/*
zlib-ng:      KSEQ_INIT(gzFile, zng_gzread)
izlib(isa-l): KSEQ_INIT(gzFile, gzread)
*/
#endif /* KSEQ_INIT */

#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define likely(expr) __builtin_expect(!!(expr), 1)

#define UVERROR(msg, code)                                                                  \
	do {                                                                                    \
		fprintf(stderr, "%s: UV(%s),%s.\n", msg, uv_err_name((code)), uv_strerror((code))); \
		assert(0);                                                                          \
	} while (0)

#define VARTYPE(N) _Generic((N), uint8_t: 1, uint16_t: 2, uint32_t: 4, uint64_t: 8, int8_t: -1, int16_t: -2, int32_t: -4, int64_t: -8, default: 0)

#define ARRAYcpySTR(d, s)                                                    \
	do {                                                                     \
		size_t src_size = sizeof(s);                                         \
		assert((MALLOCSIZE(d) ? MALLOCSIZE(d) : sizeof(d)) >= 1 + src_size); \
		char* after_c = memccpy((d), (s), '\0', src_size);                   \
		if (after_c) {                                                       \
			*after_c = '\0';                                                 \
		} else {                                                             \
			(d)[src_size] = '\0';                                            \
		}                                                                    \
	} while (0)

#define STRcpyARRAY(d, s)                                  \
	do {                                                   \
		size_t dst_size = sizeof(d);                       \
		assert(dst_size >= strlen(s));                     \
		char* after_c = memccpy((d), (s), '\0', dst_size); \
		if (after_c) {                                     \
			*after_c = '\0';                               \
		}                                                  \
	} while (0)

#define STRcpySTR(d, s)                                    \
	do {                                                   \
		size_t dst_size = MALLOCSIZE(d);                   \
		if (dst_size == 0) {                               \
			dst_size = sizeof(d);                          \
		}                                                  \
		assert(dst_size >= 1 + strlen(s));                 \
		char* after_c = memccpy((d), (s), '\0', dst_size); \
		if (after_c) {                                     \
			*after_c = '\0';                               \
		}                                                  \
		assert(after_c == NULL);                           \
	} while (0)

#define FSTREPOS_VERSION "1.0.1"
#define FSTREPOS_VERNUM 0x010001F0L /* MMNNRRSM: major minor revision status modified */
#define FSTREPOS_VER_MAJOR 1
#define FSTREPOS_VER_MINOR 0
#define FSTREPOS_VER_REVISION 1
#define FSTREPOS_VER_STATUSH 0xF /* Hex values: 0=devel, 1-E=beta, F=Release */
#define FSTREPOS_VER_MODIFIED 0  /* non-zero if modified externally from fstRePos */

#define OLDDELIMCOUNT 12
// @202305161100_Pro018_B__0516_18B_BC_R001C001_6_-1_1_1452.870483_4.477099
#define NEWDELIMCOUNT 8
// @Pro019:S:123:2:000023:R001:C001:3390.2:0010.5 1:L:0
#define MAXDELIMITEMS 12  // Old items is still 12 for contiguous delimiter,

// um
#define StepLane 11800
#define FOVStepX 429.6
#define FOVStepY 709.8
// um
#define FOV_X_MIN 1024.00
#define FOV_X_MAX 3904.34
#define FOV_Y_MIN 208.35
#define FOV_Y_MAX 1951.65

/*
y=2159   y=0,x=0              w=0,h=0      h=2159
 ⎧‾‾‾‾‾‾‾‾‾‾⎫                    ⎧‾‾‾‾‾‾‾‾‾‾⎫
 ⎪          ⎪                    ⎪          ⎪
 ⎪          ⎪                    ⎪          ⎪
 ⎪          ⎪    w = 2159 - y    ⎪          ⎪
 ⎪  R08C66  ⎪   ⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯➛   ⎪  R08C66  ⎪
 ⎪          ⎪    h = x           ⎪          ⎪
 ⎪          ⎪                    ⎪          ⎪
 ⎪          ⎪                    ⎪          ⎪
 ⎩__________⎭                    ⎩__________⎭
           x=4095               h=4095
*/

// Px
#define FOVWidth 2160
// 2160*3.45/14:  532.2857 um X 0.807  = 429.55457 um ≈ FOVStepX; FOVStepX/3.45*14 ≈ 1743.3043px
#define FOVHight 4096
// 4096*3.45/14: 1009.3714 um X 0.7032 = 709.78999 um ≈ FOVStepY; FOVStepY/3.45*14 ≈ 2880.3478px
#define FOV_USED_WIDTH (FOVStepX * 14.0 / 3.45)
#define FOV_USED_HEIGHT (FOVStepY * 14.0 / 3.45)

/*
⎧‾‾‾‾‾⎫ . . . ⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫ . . . ⎧‾‾‾‾‾‾⎫
⎪R1C01⎪ . . . ⎪R1C63⎪⎪R1C64⎪⎪R1C65⎪⎪R1C66⎪⎪R1C66⎪⎪R1C68⎪ . . . ⎪R1C130⎪
⎩_____⎭ . . . ⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭ . . . ⎩______⎭
 . . .                              . . .                        . . .
 . . .                              . . .                        . . .
 . . .                              . . .                        . . .
⎧‾‾‾‾‾⎫ . . . ⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫ . . . ⎧‾‾‾‾‾‾⎫
⎪R7C01⎪ . . . ⎪R7C63⎪⎪R7C64⎪⎪R7C65⎪⎪R7C66⎪⎪R7C66⎪⎪R7C68⎪ . . . ⎪R7C130⎪
⎩_____⎭ . . . ⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭ . . . ⎩______⎭
⎧‾‾‾‾‾⎫ . . . ⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧<‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫ . . . ⎧‾‾‾‾‾‾⎫
⎪R8C01⎪ . . . ⎪R8C63⎪⎪R8C64⎪⎪R8C65⎪⎪R8C66⎪⎪R8C66⎪⎪R8C68⎪ . . . ⎪R8C130⎪
⎩_____⎭ . . . ⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭ . . . ⎩______⎭
 . . .                              . . .                        . . .
 . . .                              . . .                        . . .
 . . .                              . . .                        . . .
⎧‾‾‾‾‾⎫ . . . ⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫⎧‾‾‾‾‾⎫ . . . ⎧‾‾‾‾‾‾‾⎫
⎪R14C1⎪ . . . ⎪REC63⎪⎪REC64⎪⎪REC65⎪⎪REC66⎪⎪REC66⎪⎪REC68⎪ . . . ⎪R14C130⎪
⎩_____⎭ . . . ⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭⎩_____⎭ . . . ⎩_______⎭
*/
// #define LaneRow0 15
#define LaneCol 130
#define LaneRow 14
#define CenterFOV_ROW 8
#define CenterFOV_COL 66

#define MAXCPUCORES 1024
#define JOBQUEUESIZE ((MAXCPUCORES * 3) >> 1)
// 任务队列数设为CPU最大核心数的1.5倍，则肯定塞满 workers 后会有剩的。
#define JOBITEMSIZE 4

#define GZBUFSIZE (131072 * 2)
#define MAXFQIDLEN 80
#define BARCODELEN 30
#define MINBARCODELEN 3
#define MAXPOSTRLEN 10
/* -163473.62 FullShape=(40141, 224640). Size is same of 2 Lanes for both Solo and Dual */
#define MAXCOORDSTRLEN (9 + 1 + 10)
#define MAXSPATIALEN (BARCODELEN + 1 + MAXCOORDSTRLEN)
/* MAXSPATIALEN+1 == BARCODELEN+1+MAXCOORDSTRLEN+1 */
#define ROWCOLSIZE 8 /* R014C130 */
/*
https://stackoverflow.com/questions/3767284/using-printf-with-a-non-null-terminated-string
    printf("%.*s", stringLength, pointerToString);
    fwrite(your_string, sizeof(char), number_of_chars, stdout);
*/
#pragma pack(1)
// #pragma pack(push, 1)
struct fstBCdata_s {
	int8_t name[MAXFQIDLEN];
#ifndef RELEASE
	char* comment;
#endif
	int8_t seq[BARCODELEN];
	int8_t qual[BARCODELEN];
	int8_t RowCol[ROWCOLSIZE];
	double newXY[2];
};
typedef struct fstBCdata_s fstBCdata_t;
/*
union fstBCoutput_u {
    int8_t SpatiaStr[BARCODELEN + 1 + MAXCOORDSTRLEN];
    struct SpatiaDat_s {
        int8_t seq[BARCODELEN];
        int8_t delim;
        int8_t xy[MAXCOORDSTRLEN];
    } SpatiaDat;
};
typedef union fstBCoutput_u fstBCoutput_t;
*/
// #pragma pack(pop)
#pragma pack()

struct workerArray_s {
	int_least16_t workerID;
	uint64_t fqSliceID;
	atomic_int_least8_t flag;  // one day, it will be useful.
	fstBCdata_t jobDatArray[JOBITEMSIZE];
	char* tokens[MAXDELIMITEMS];
	// fstBCoutput_t output_array[JOBITEMSIZE];
};
typedef struct workerArray_s workerArray_t;

#define PARAMETERS_BUFFER_SIZE 120
struct parameters_s {
	// ArrayStateEnum_t jobDataState;
	// ArrayStateEnum_t outDataState;
	char* inFastqFilename;
	float unZoomRatio;  // 1 or 1.25, float is (1,8,23), thus enough for [- 2^{23} + 1, 2^{23} - 1]
	kseq_t* kseq;
	gzFile ksfp;
	regex_t regex;
	uint64_t fqRead;
	uint64_t fqValid;
	uint64_t fqSkipped;
	char buffer[PARAMETERS_BUFFER_SIZE];
	workerArray_t* worksQuene;
};
typedef struct parameters_s parameters_t;
extern parameters_t Parameters;
// extern struct parameters_s parameters;

void fqReader_init(void);
void fqReader_destroy(void);
void fill_worker(int_least16_t worker_id);
void worker(int_least16_t worker_id);
void output_worker(int_least16_t worker_id);

/* ################ */

#ifdef __cplusplus
}
#endif
#endif /* FSTRECOORD_COMMON_H_ */
