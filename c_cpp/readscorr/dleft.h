// by Hu Xuesong
#ifndef _G_DLEFT_H
#define _G_DLEFT_H

#include <stdint.h>

typedef struct __DLeftArray_t {
    unsigned char CountBit;
    unsigned char rBit;
    size_t ArraySize;
    unsigned char ArrayCount;
    uint64_t ItemInsideAll, CellOverflowCount; // ItemInsideAll = ItemInsideArray + CellOverflowCount
    double FalsePositiveRatio;
    char *dlap, *extreep;
} DLeftArray_t;

#ifndef KMER_T
#define KMER_T
typedef struct __dbitseq_t {
    size_t l, m;
    char *s;
} dBitSeq_t;
#endif

DLeftArray_t *dleft_arrayinit(unsigned char CountBit, unsigned char rBit, size_t ArraySize, unsigned char ArrayCount);


#endif /* dleft.h */

