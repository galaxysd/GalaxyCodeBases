#include <stdlib.h>	//calloc
#include <stdint.h>	// uint_fast8_t
#include "dleft.h"
#include "MurmurHash3.h"

//#ifndef ROUNDUPMOD8
//#define ROUNDUPMOD8(x) ( ((x)+7) & ~(x) )
//#endif
#define	FORCE_INLINE static inline __attribute__((always_inline))

DLeftArray_t *dleft_arrayinit(unsigned char CountBit, unsigned char rBit, size_t ArraySize, unsigned char ArrayCount) {
    unsigned char itemByte = (CountBit+rBit+7) >> 3;
    DLeftArray_t *dleftobj = malloc(sizeof(DLeftArray_t));
    dleftobj->dlap = calloc(ArrayCount,itemByte*ArraySize);
    dleftobj->CountBit = CountBit;
    dleftobj->rBit = (itemByte<<3) - CountBit;
    dleftobj->ArraySize = ArraySize;
    dleftobj->ArrayCount = ArrayCount;
    return dleftobj;
}

FORCE_INLINE int_fast8_t dleft_insert_kmer(dBitSeq_t * str, DLeftArray_t * dleftobj) {
}

int_fast8_t dleft_insert_read(dBitSeq_t * str, DLeftArray_t * dleftobj) {
}
