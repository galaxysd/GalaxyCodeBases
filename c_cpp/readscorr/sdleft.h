// by Hu Xuesong
#ifndef _G_SDLEFT_H
#define _G_SDLEFT_H

#include <stdint.h>

typedef struct __DLeftArray_t {
    unsigned char CountBit, rBit, ArrayBit;
    unsigned char HashCnt;
    size_t ArraySize;
    //unsigned char ArrayCount;
    uint64_t ItemInsideAll, CellOverflowCount; // ItemInsideAll = ItemInsideArray + CellOverflowCount
    double FalsePositiveRatio;
    char *dlap, *extreep;
    uint64_t *outhash;
} SDLeftArray_t;

/*
#ifndef KMER_T
#define KMER_T
typedef struct __dbitseq_t {
    size_t l, m;
    char *s;
} dBitSeq_t;
#endif
*/

SDLeftArray_t *dleft_arrayinit(unsigned char CountBit, unsigned char rBit, size_t ArraySize);
int_fast8_t dleft_insert_read(char const *const inseq, size_t len, SDLeftArray_t *dleftobj);

void fprintSDLAnfo(FILE *stream, const SDLeftArray_t * dleftobj);
void dleft_arraydestroy(SDLeftArray_t * const dleftobj);

#endif /* sdleft.h */

