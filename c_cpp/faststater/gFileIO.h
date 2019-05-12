// by Hu Xuesong
#ifndef _G_FILEIO_H
#define _G_FILEIO_H

#include "gFileType.h"	// for file ID
#include <stdint.h>	// uint64_t
//#include "sdleft.h"

#define GFIOCHRBASE 1
#define GFIODIBBASE 2

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;	// l for sequence length in base-pairs, m for malloc in bytes.
	char *s;	// for 2bit ATCG, 1bp=2bit; for 6bit quality, 1bp=6bit.
} kstring_t;
#endif	// Yes, it is the same as that in kseq.h

//typedef ssize_t (*G_ssize_t_oneIN)(void * const);
//typedef void (*G_p_oneIN)(void * const);
struct __SeqFileObj;
typedef ssize_t (G_ssize_t_oneIN)(struct __SeqFileObj * const);
typedef void (G_p_oneIN)(struct __SeqFileObj * const);

typedef struct __SeqFileObj {
   size_t readlength,binNcount,binMallocedQQWord;
   const char *name, *comment, *seq, *qual;
	// $ cdecl explain "char * const* name"
	// declare name as pointer to const pointer to char
   uint64_t *diBseq;   // 2bit, {A,C,G,T}={0,1,2,3}
   unsigned char *hexBQ;   // 0~63 for Quality, 128 for N, 64 for Eamss-masked or smallcase-masked
   //int (*getNextSeq)(void *);   // void * fh
   G_ssize_t_oneIN *getNextSeq;
   G_p_oneIN *closefh; // remember to close file handle
   //= void (*closefh)(struct __SeqFileObj * const);
   void *fobj;   // Not just FILE *fp
   long datePos[2];   // [1] for item id, [0] for offset. For Text file, item id == 0.
   //long datePosOffset,datePosItem;
   uint_fast8_t type;   // 1->hasBaseChr, 2->hasQchar, 4->hasBase2bit,
   // 8->hashexBQ or just array of {0,128} for N,
   // 16->Eamss bit set, 32->Q value not raw, after Eamss.
} SeqFileObj;   // We may support both FA/FQ and binary formats. So, object is a good thing.

SeqFileObj * inSeqFinit(const char * const, unsigned char);
//ssize_t inSeqFreadNext(SeqFileObj * const);
int inSeqFseek(SeqFileObj * const, const fpos_t datePos[]);
void inSeqFdestroy(SeqFileObj * const);

#endif /* gFileIO.h */
