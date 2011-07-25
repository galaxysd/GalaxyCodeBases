#define _GNU_SOURCE
#include "gtypendef.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
//#include <errno.h>
#include <err.h>
#include <argp.h>
#include "MurmurHash3.h"
#include "getch.h"
#include "2bitarray.h"
#include "gFileIO.h"
#include "sdleft.h"
#include "chrseq.h"

// http://www.mersenneforum.org/showthread.php?t=14419
#include <endian.h>

typedef union
{
    __float128  val;
    float128    ourval;
    struct
    {
#if __BYTE_ORDER == __BIG_ENDIAN
        uint64_t sign    : 1;
        uint64_t exp     : 15;
        uint64_t frac1   : 48;
        uint64_t frac0   : 64;
#else
        uint64_t frac0   : 64;
        uint64_t frac1   : 48;
        uint64_t exp     : 15;
        uint64_t sign    : 1;
#endif
    } bits;
} qfloat;


void test_types(void) {
    union {
        uint128_t u128;
        uint64_t u64[2];
    } test1;
    test1.u128 = 16u+(uint128_t)UINT64_MAX*65536u;
    printf("uint128_t [%lx %lx]\n",test1.u64[1],test1.u64[0]);
// http://www.mersenneforum.org/showthread.php?t=14419
    qfloat      foo;
    qfloat      bar;
    foo.val = 65536.0;
    foo.val += 1.0;
    foo.val *= foo.val;
    bar.ourval = (__float128)1/(__float128)3;
    printf("%04X %012lX %016lX\n", (uint16_t)foo.bits.exp, (uint64_t)foo.bits.frac1, foo.bits.frac0);
    printf("%04X %012lX %016lX\n", (uint16_t)bar.bits.exp, (uint64_t)bar.bits.frac1, bar.bits.frac0);
}

void test_baseswitcher(void) {
    const unsigned char baseswitcher[]={
        4,0,4,1,4,   // 64-68
        4,4,2,4,4,   // 69-73
        4,4,4,4,4,   // 74-78
        4,4,4,4,4,   // 79-83
        3            // 84
    };
    char *seq="ATCGatcgNn";
    while (*seq) {
        printf("%c -> %u\n",*seq,baseswitcher[*seq-64]);
        ++seq;
    }
}

int main(void) {
    test_types();
    test_baseswitcher();
    return 0;
}
