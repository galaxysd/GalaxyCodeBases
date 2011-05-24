#ifndef _G_CHRSEQ_H
#define _G_CHRSEQ_H

#include <stdint.h>	// int_fast8_t
#include <stdlib.h> // malloc

#ifndef FORCE_INLINE
#define	FORCE_INLINE static inline __attribute__((always_inline))
#endif

FORCE_INLINE void NormalizeChrSeq(char *seq){
    while(*seq){
        switch (*seq) {
	    case 'a': case 'A':
		    *seq = 'A';
		    break;
	    case 't': case 'T':
		    *seq = 'T';
		    break;
	    case 'c': case 'C':
		    *seq = 'C';
		    break;
	    case 'g': case 'G':
		    *seq = 'G';
		    break;
	    default:
            *seq = 'N';
		    break;
        }   // Whether a looking-up array[64] can be faster ?
        ++seq;
    }
}

// *seq must have been Normalized to /[ATCGN]*/
FORCE_INLINE char *ChrSeqRevComp(char const * seq, size_t len){
    char *revcomseq=malloc(len+1);
    const char *tmpseqin=seq+len-1; // if len is shorter than real length, trim it.
    char *tmpseqrc=revcomseq;
    while(tmpseqin>=seq){
        switch (*tmpseqin) {
	    case 'A':
		    *tmpseqrc++ = 'T';
		    break;
	    case 'T':
		    *tmpseqrc++ = 'A';
		    break;
	    case 'C':
		    *tmpseqrc++ = 'G';
		    break;
	    case 'G':
		    *tmpseqrc++ = 'C';
		    break;
	    default:
            *tmpseqrc++ = 'N';
		    break;
        }   // Whether a looking-up array[32] can be faster ?
        --tmpseqin;
    }
    *tmpseqrc++ = '\0';
    return revcomseq;
}

#endif  // chrseq.h

