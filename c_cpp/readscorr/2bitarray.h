// by Hu Xuesong

#ifndef _GA_TBITARRAY_H
#define _GA_TBITARRAY_H

#include "gtypendef.h"
#include <stdint.h>
#include <stdlib.h>

#include <limits.h>		/* for CHAR_BIT */

#define CHAR_TBIT (CHAR_BIT>>1)	// 8/2 = 4

#define TBITSHIFT(b) (((b) % CHAR_TBIT)<<1)

#define TBITMASK(b) /* Returns the bit mask of a byte */ \
    (3 << TBITSHIFT(b))

#define TBITSLOT(b) ((b) / CHAR_TBIT) // Just the index of char[]

#define TBITGETVALUE(a, b) /* Get Value */ \
    (((a)[TBITSLOT(b)] >> TBITSHIFT(b)) & 3)
//  (((a)[2BITSLOT(b)] & 2BITMASK(b)) >> 2BITSHIFT(b))

#define TBITNSLOTS(nb) ((2*(nb) + CHAR_BIT - 1) / CHAR_BIT) // caltulate the length of char[]. x*2 can overflow while x<<1 cannot.

/*
__inline void TBITSetValue ( char arr[], size_t index, uint_fast8_t value );
__inline uint_fast8_t TBITSaturatedADD ( char arr[], size_t index, int_fast8_t value );
__inline uint_fast8_t TBITSaturatedINC ( char arr[], size_t index );
__inline uint_fast8_t TBITSaturatedDEC ( char arr[], size_t index );
*/

// I find gcc will not inline functions in other c files, thus move them here.

FORCE_INLINE void TBITSetValue ( char arr[], size_t index, uint_fast8_t value ) {
    if (value>3) value=3;
    unsigned char theByte=(arr[TBITSLOT(index)]) & ~TBITMASK(index);
    theByte |= value << TBITSHIFT(index);
    arr[TBITSLOT(index)] = theByte;
}

FORCE_INLINE uint_fast8_t TBITSaturatedADD ( char arr[], size_t index, int_fast8_t value ) {
    unsigned char theByte = (arr[TBITSLOT(index)]);
    int_fast8_t theValue = (theByte >> TBITSHIFT(index)) & 3 ;
    if (value==0) {
        return (uint_fast8_t) theValue;
    } else if (value > 0) {
        if (theValue > 3-value) theValue=3;
          else theValue += value;
    } else {
        if (theValue < -value) theValue =0;
          else theValue += value;
    }
    theByte &= ~TBITMASK(index);
    theByte |= theValue << TBITSHIFT(index);
    arr[TBITSLOT(index)] = theByte;
    return (uint_fast8_t) theValue;
}

FORCE_INLINE uint_fast8_t TBITSaturatedINC ( char arr[], size_t index ) {
    unsigned char theByte = (arr[TBITSLOT(index)]);
    uint_fast8_t theValue = (theByte >> TBITSHIFT(index)) & 3 ;
    if (theValue==3) {
        return 3;
    } else {
        ++theValue;
        theByte &= ~TBITMASK(index);
        theByte |= theValue << TBITSHIFT(index);
        arr[TBITSLOT(index)] = theByte;
        return theValue;
    }
}

FORCE_INLINE uint_fast8_t TBITSaturatedDEC ( char arr[], size_t index ) {
    unsigned char theByte = (arr[TBITSLOT(index)]);
    uint_fast8_t theValue = (theByte >> TBITSHIFT(index)) & 3 ;
    if (theValue==0) {
        return 0;
    } else {
        --theValue;
        theByte &= ~TBITMASK(index);
        theByte |= theValue << TBITSHIFT(index);
        arr[TBITSLOT(index)] = theByte;
        return theValue;
    }
}



#endif /* 2bitarray.h */

// Original comments for bitarray: 
/*
(If you don't have <limits.h>, try using 8 for CHAR_BIT.)

Here are some usage examples. To declare an ``array'' of 47 bits:

	char bitarray[BITNSLOTS(47)];

To set the 23rd bit:

	BITSET(bitarray, 23);

To test the 35th bit:

	if(BITTEST(bitarray, 35)) ...

To compute the union of two bit arrays and place it in a third array 
(with all three arrays declared as above):

	for(i = 0; i < BITNSLOTS(47); i++)
		array3[i] = array1[i] | array2[i];

To compute the intersection, use & instead of |.

As a more realistic example, here is a quick implementation of the Sieve 
of Eratosthenes, for computing prime numbers:

#include <stdio.h>
#include <string.h>

#define MAX 10000

int main()
{
	char bitarray[BITNSLOTS(MAX)];
	int i, j;

	memset(bitarray, 0, BITNSLOTS(MAX));

	for(i = 2; i < MAX; i++) {
		if(!BITTEST(bitarray, i)) {
			printf("%d\n", i);
			for(j = i + i; j < MAX; j += i)
				BITSET(bitarray, j);
		}
	}
	return 0;
}
*/

