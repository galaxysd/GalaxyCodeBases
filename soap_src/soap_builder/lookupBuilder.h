#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#define Kibi *1024
#define Mibi *1024*1024
#define Gibi *1024ll*1024*1024
#ifndef MAX
#define MAX(a, b) ((a)<(b)?(b):(a))
#endif
#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif
#define TRY(a) (assert((a) == -1 ? (printf("Error: [%d] %s\n", errno, strerror(errno)), 0) : 1))
#define TRYEQ(a, b) (assert((a) != b ? (printf("Error: [%d] %s\n", errno, strerror(errno)), 0) : 1))
typedef unsigned long long ULL;
