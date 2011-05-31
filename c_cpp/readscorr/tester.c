#define _GNU_SOURCE
#include "gtypendef.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <err.h>
#include <argp.h>
#include "MurmurHash3.h"
#include "getch.h"
#include "2bitarray.h"
#include "gFileIO.h"
#include "sdleft.h"
#include "chrseq.h"

void test_types(void) {
    union {
        uint128_t u128;
        uint64_t u64[2];
    } test1;
    test1.u128 = 16u+(uint128_t)UINT64_MAX*65536u;
    printf("uint128_t [%lx %lx]\n",test1.u64[1],test1.u64[0]);
}

int main(void) {
    test_types();
    return 0;
}
