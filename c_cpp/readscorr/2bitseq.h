#ifndef _G_DBSEQ_H
#define _G_DBSEQ_H

#include <stdint.h>	// int_fast8_t

#ifndef FORCE_INLINE
#define	FORCE_INLINE __attribute__((always_inline))
#endif

size_t Ncount;

FORCE_INLINE unsigned char singlebase2dbit(const char base, unsigned char *const hexBQchr) {
	switch (base) {
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

FORCE_INLINE int_fast8_t base2dbit(const size_t seqlen,
 const char *const baseschr, const char *const qstr,
 unsigned char *diBseq, unsigned char *hexBQ) {
// First, do the bases, overwrite output arrays since we have just malloc without setting to zero. 
    Ncount=0;
    size_t i;
    for (i=0;i<=seqlen;i+=4) {
        *diBseq++ = (singlebase2dbit(i,hexBQ+i)<<6)       |
                    ((singlebase2dbit(i+1,hexBQ+i+1)&3)<<4) |
                    ((singlebase2dbit(i+2,hexBQ+i+2)&3)<<2) |
                     (singlebase2dbit(i+3,hexBQ+i+3)&3);
        // should be better than diBseq[i/4]
    }
    // now, i=seqlen+1
    unsigned char tmpqdbase = 0;
    switch (seqlen % 4) {
    case 3:
        tmpqdbase |= (singlebase2dbit(i+2,hexBQ+i+2)&3)<<6;
    case 2:
        tmpqdbase |= (singlebase2dbit(i+1,hexBQ+i+1)&3)<<4;
    case 1:
        tmpqdbase |= (singlebase2dbit(i,hexBQ+i)&3)<<2;
        *diBseq++ = tmpqdbase;
    default:
        // Everything should done in for-cycle when seqlen%4 == 0.
        // No need to add '/0' since we have got seqlen, and there may be poly-A.
        break;
    }
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




#endif  // 2bitseq.h

