#ifndef _G_DBSEQ_INLINE_H
#define _G_DBSEQ_INLINE_H

#include <stdint.h>	// int_fast8_t
#include <stdlib.h> // malloc
//#include "2bitseq.h"
#ifdef DEBUG
#include <stdio.h>
#endif

#ifndef FORCE_INLINE
#define	FORCE_INLINE static inline __attribute__((always_inline))
#endif

FORCE_INLINE uint64_t *dibmalloc(size_t len){
    size_t needtomallocQQW = (len+31u)>>5;
    uint64_t *outseq = malloc(needtomallocQQW*8);
    return outseq;
}

FORCE_INLINE uint64_t singlebase2dbitplus(const char *const base, unsigned char *const hexBQchr, size_t *const Ncountp) {
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
    	*hexBQchr = 1u<<7;    // 128 for N
    	++(*Ncountp);
    	// DONOT use `|=` since the memory is just malloced
		return 0;
		break;
	}   // Whether a looking-up array[64] can be faster ?
}

FORCE_INLINE size_t base2dbit(size_t seqlen,
 const char *const baseschr, const char *const qstr,
 uint64_t *diBseq, unsigned char *hexBQ) {
// First, do the bases, overwrite output arrays since we have just malloc without setting to zero. 
    size_t Ncount=0;
    size_t i,j;
	const size_t seqlenDowntoDW = seqlen & ~((2u<<4)-1);
#ifdef DEBUG
 printf("[a]%zu %zu\n",seqlen,seqlenDowntoDW);
#endif
	uint64_t tmpqdbase;
    for (i=0;i<seqlenDowntoDW;i+=32u) {
        tmpqdbase=0;
        for (j=0;j<32;j++) {
            //tmpqdbase |= singlebase2dbit(baseschr+i+j,hexBQ+i+j)<<(62u-j);
            tmpqdbase |= singlebase2dbitplus(baseschr+i+j,hexBQ+i+j,&Ncount)<<(j*2);
// printf(" A{%lx,%.1s,%lx}",tmpqdbase,baseschr+i+j,singlebase2dbit(baseschr+i+j,hexBQ+i+j)<<(j*2));
        }
        *diBseq++ = tmpqdbase;
    }
// printf("[b]%zu %zu\n",i,j);
    tmpqdbase = 0;
    for (j=0;j<seqlen-i;j++) {   // seqlen starts from 0.
// Cannot use 'j<=seqlen-i-1' here since 'size_t j' cannot be negative.
        tmpqdbase |= singlebase2dbitplus(baseschr+i+j,hexBQ+i+j,&Ncount)<<(j*2);
// printf(" B{%lx,%.1s,%lx}",tmpqdbase,baseschr+i+j,singlebase2dbit(baseschr+i+j,hexBQ+i+j)<<(j*2));
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
#if !defined(__GNUC__) || (__GNUC__ < 4) || (__GNUC__ == 4 && __GNUC_MINOR__ < 3)
// see also https://github.com/damm/ffi/raw/7d11ff675759ca47243bb901e4226ee3ec6de0c6/ext/ffi_c/AbstractMemory.c
	seq32mer = ((seq32mer & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq32mer & 0xFF00FF00FF00FF00LLU)>> 8);
	seq32mer = ((seq32mer & 0x0000FFFF0000FFFFLLU)<<16) | ((seq32mer & 0xFFFF0000FFFF0000LLU)>>16);
	seq32mer = ((seq32mer & 0x00000000FFFFFFFFLLU)<<32) | ((seq32mer & 0xFFFFFFFF00000000LLU)>>32);
#else
	seq32mer = __builtin_bswap64(seq32mer);
#endif
	return seq32mer;
}
/*
0xABCDEF0123456789 -> cgagtcgcccactagacaaattgtctattggg
revcomp ->            cccaatagacaatttgtctagtgggcgactcg
*/
//FORCE_INLINE uint64_t QQWkmerMovHigher(uint64_t *seq32mer, unsigned char bphigher, uint64_t ){}

// *base must have been Normalized to /[ATCGN]*/
FORCE_INLINE uint64_t singlebase2dbit(const char *const base, size_t *const Ncountp) {
	switch (*base) {
	case 'A':
		return 0;
		break;
	case 'T':
		return 3;
		break;
	case 'C':
		return 1;
		break;
	case 'G':
		return 2;
		break;
	default:
    	++(*Ncountp);
    	// DONOT use `|=` since the memory is just malloced
		return 0;
		break;
	}   // Whether a looking-up array[32] can be faster ?
}
FORCE_INLINE size_t ChrSeq2dib(const char *const baseschr, size_t seqlen, uint64_t **diBseqp, size_t *const puint64cnt){
    size_t needtomallocQQW = (seqlen+31u)>>5;  // in fact length in sizeof(uint64_t)
    if (needtomallocQQW > *puint64cnt) {
        *diBseqp = realloc(*diBseqp, needtomallocQQW*8);
        *puint64cnt = needtomallocQQW;
    }
    uint64_t *tmpdiBseqp = *diBseqp;
    size_t Ncount = 0;
    size_t i,j;
	const size_t seqlenDowntoDW = seqlen & ~((2u<<4)-1);
	uint64_t tmpqdbase;
    for (i=0;i<seqlenDowntoDW;i+=32u) {
        tmpqdbase=0;
        for (j=0;j<32;j++) {
            tmpqdbase |= singlebase2dbit(baseschr+i+j,&Ncount)<<(j*2);
        }
        *tmpdiBseqp++ = tmpqdbase;
    }
    tmpqdbase = 0;
    for (j=0;j<seqlen-i;j++) {   // seqlen starts from 0.
        tmpqdbase |= singlebase2dbit(baseschr+i+j,&Ncount)<<(j*2);
    }
	*tmpdiBseqp++ = tmpqdbase;
#ifdef EP	
	while (tmpdiBseqp < (*diBseqp)+*puint64cnt) {
	    *tmpdiBseqp++ = 0LLU;
	}
#endif  // not nessary since we have seqlen and the last uint64 is always with tailing 0s.
#ifdef DEBUG
 printf("[a]%zu %zu [%lx] [%s]\n",seqlen,seqlenDowntoDW,**diBseqp,baseschr);  // just show the 1st one.
#endif
    return Ncount;
}

#endif  // 2bitseqinline.h

