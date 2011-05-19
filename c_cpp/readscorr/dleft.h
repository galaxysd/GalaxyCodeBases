// by Hu Xuesong
#ifndef _G_DLEFT_H
#define _G_DLEFT_H

#include <stdint.h>

typedef struct __DLeftArray_t {
    unsigned char CountBit, rBit, ArrayBit;
    unsigned char HashCnt;
    size_t ArraySize;
    unsigned char ArrayCount;
    uint64_t ItemInsideAll, CellOverflowCount; // ItemInsideAll = ItemInsideArray + CellOverflowCount
    double FalsePositiveRatio;
    char *dlap, *extreep;
} DLeftArray_t;

/*
#ifndef KMER_T
#define KMER_T
typedef struct __dbitseq_t {
    size_t l, m;
    char *s;
} dBitSeq_t;
#endif
*/

DLeftArray_t *dleft_arrayinit(unsigned char CountBit, unsigned char rBit, size_t ArraySize, unsigned char ArrayCount);
int_fast8_t dleft_insert_read(char const *const inseq, size_t len, DLeftArray_t *dleftobj);

void fprintDLAnfo(FILE *stream, const DLeftArray_t * dleftobj);
void dleft_arraydestroy(DLeftArray_t * const dleftobj);

#endif /* dleft.h */

