// From http://c-faq.com/misc/bitsets.html
// See also /usr/include/X11/extensions/xtrapbits.h
// Added BitToggle and BitCopy with ByteInArray from xtrapbits.h. Not tested yet.
// by Hu Xuesong

#include <limits.h>		/* for CHAR_BIT */

typedef unsigned char * UByteP;  /* Pointer to an unsigned byte array */

#define BITMASK(b) /* Returns the bit mask of a byte */ \
    (1 << ((b) % CHAR_BIT))

#define ByteInArray(array,bit) /* Returns the byte offset to get to a bit */ \
    (((UByteP)(array))[(bit) / CHAR_BIT)])

#define BITSLOT(b) ((b) / CHAR_BIT) // Just the index of char[]

#define BITSET(a, b) /* Set a specific bit to be True */ \
    ((a)[BITSLOT(b)] |= BITMASK(b))

#define BITCLEAR(a, b) /* Set a specific bit to be False */ \
    ((a)[BITSLOT(b)] &= ~BITMASK(b))

#define BitToggle(array,bit)    /* Toggle a specific bit */ \
    (ByteInArray(array,bit) ^= BITMASK(bit))

#define BITTEST(a, b) /* Test to see if a specific bit is True */ \
    ((a)[BITSLOT(b)] & BITMASK(b))

#define BitCopy(dest,src,bit)   /* Copy a specific bit */ \
    BITTEST((src),(bit)) ? BITSET((dest),(bit)) : BITCLEAR((dest),(bit))

#define BITNSLOTS(nb) ((nb + CHAR_BIT - 1) / CHAR_BIT) // caltulate the length of char[]

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

