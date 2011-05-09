#ifndef _G_DBSEQ_H
#define _G_DBSEQ_H

#include <stdint.h>	// int_fast8_t

#ifndef FORCE_INLINE
#define	FORCE_INLINE __attribute__((always_inline))
#endif

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
    	*hexBQchr |= 1<<7;    // 128 for N
		return 0;
		break;
	}
}

FORCE_INLINE int_fast8_t base2dbit(const size_t seqlen,
 const char *const baseschr, const char *const qstr,
 unsigned char *diBseq, unsigned char *hexBQ) {
    size_t Ncount=0;
    size_t i;
    unsigned char tmpqdbase;
    for (i=0;i<=seqlen;i+=4) {
        tmpqdbase = (singlebase2dbit(i,hexBQ+i)<<6)       |
                    ((singlebase2dbit(i+1,hexBQ+i+1)&3)<<4) |
                    ((singlebase2dbit(i+2,hexBQ+i+2)&3)<<2) |
                     (singlebase2dbit(i+3,hexBQ+i+3)&3);
        diBseq[i/4]=tmpqdbase;
    }
    tmpqdbase = 0;
    switch (seqlen % 4) {
    case 3:
    case 2:
    case 1:
    default:
        break;
    }
    return Ncount;
}




#endif  // 2bitseq.h

