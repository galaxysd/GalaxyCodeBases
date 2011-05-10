#ifndef _G_DBSEQ_INLINE_H
#define _G_DBSEQ_INLINE_H

#include <stdint.h>	// int_fast8_t
#include <stdlib.h> // malloc
#include "2bitseq.h"

#ifndef FORCE_INLINE
#define	FORCE_INLINE __attribute__((always_inline))
#endif

size_t Ncount;

FORCE_INLINE unsigned char singlebase2dbit(const char *const base, unsigned char *const hexBQchr) {
	switch (*base) {
	case 'a': case 'A':
		return 0;
		break;
	case 't': case 'T':
		return 3;
		break;
	case 'c': case 'C':
		return 1;
		break;
	case 'g': case 'G':
		return 2;
		break;
	default:
    	*hexBQchr = 1<<7;    // 128 for N
    	++Ncount;
    	// DONOT use `|=` since the memory is just malloced
		return 0;
		break;
	}
}

FORCE_INLINE int_fast8_t base2dbit(size_t seqlen,
 const char *const baseschr, const char *const qstr,
 uint64_t *diBseq, unsigned char *hexBQ) {
// First, do the bases, overwrite output arrays since we have just malloc without setting to zero. 
    Ncount=0;
    size_t i;
    uint_fast8_t j;
	const size_t seqlenDowntoDW = seqlen & ~((2<<6)-1);
	uint64_t tmpqdbase;
    for (i=0;i<seqlenDowntoDW;i+=32u) {
        tmpqdbase=0;
        for (j=0;j<64;j+=2) {
            tmpqdbase |= singlebase2dbit(baseschr+i+j,hexBQ+i+j)<<(62u-j);
        }
/*
        *diBseq++ = (singlebase2dbit(i,hexBQ+i)<<6)       |
                    ((singlebase2dbit(i+1,hexBQ+i+1)&3)<<4) |
                    ((singlebase2dbit(i+2,hexBQ+i+2)&3)<<2) |
                     (singlebase2dbit(i+3,hexBQ+i+3)&3);
*/
        *diBseq++ = tmpqdbase;
        // should be better than diBseq[i/32]
    }
    tmpqdbase = 0;
    for (j=0;j<=seqlen-i-1;j++) {   // seqlen starts from 0
        tmpqdbase |= singlebase2dbit(baseschr+i+j,hexBQ+i+j)<<(62u-j*2);
        //printf("{%lx,%s,%lx}",tmpqdbase,baseschr+i+j,singlebase2dbit(baseschr+i+j,hexBQ+i+j)<<((31u-j)<<1));
    }
	*diBseq++ = tmpqdbase;
    // Everything should done in for-cycle when seqlen%4 == 0.
    // No need to add '/0' since we have got seqlen, and there may be poly-A.
// Then, for the Q values ...
    if (qstr != NULL) {
        for (i=0;i<=seqlen;i++) {
            if ( qstr[i] >= '@' && qstr[i] <= 'h' )
                tmpqdbase = qstr[i] - '?';  // @=1, A=2, B=3
            else tmpqdbase = 0;
            hexBQ[i] |= tmpqdbase;  // lower 7bit already Zero-ed.
        }
    }
// Eamss-masking will be considered later.    
    return Ncount;
}

FORCE_INLINE uint64_t unitReverseComp(uint64_t seq32mer){
	seq32mer = ~seq32mer;
	seq32mer = ((seq32mer & 0x3333333333333333LLU)<< 2) | ((seq32mer & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
	seq32mer = ((seq32mer & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq32mer & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
	seq32mer = ((seq32mer & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq32mer & 0xFF00FF00FF00FF00LLU)>> 8);
	seq32mer = ((seq32mer & 0x0000FFFF0000FFFFLLU)<<16) | ((seq32mer & 0xFFFF0000FFFF0000LLU)>>16);
	seq32mer = ((seq32mer & 0x00000000FFFFFFFFLLU)<<32) | ((seq32mer & 0xFFFFFFFF00000000LLU)>>32);
	return seq32mer;
}

#endif  // 2bitseq.h

