// by Hu Xuesong
#ifndef _G_FILEIO_H
#define _G_FILEIO_H

#include <stdint.h>

#define GFIOCHRBASE 1
#define GFIODIBBASE 2

#define G_HEADER_LENGTH 4
#define G_TYPE_FAQC "FAQC"
#define G_TYPE_GZ "0x1f0x8b0x08"


struct __SeqFileObj;
typedef ssize_t (G_ssize_t_oneIN)(struct __SeqFileObj * const);
typedef void (G_p_oneIN)(struct __SeqFileObj * const);

typedef struct __SeqFileObj {
   size_t readlength,binNcount,binMallocedQQWord;
   const char *name, *comment, *seq, *qual;
   uint64_t *diBseq;   // 2bit, {A,C,G,T}={0,1,2,3}
   unsigned char *hexBQ;   // 0~63 for Quality, 128 for N, 64 for Eamss-masked or smallcase-masked
   G_ssize_t_oneIN *getNextSeq;
   G_p_oneIN *closefh;
   void *fobj;
   long datePos[1];   // [1] for item id, [0] for offset. For Text file, item id == 0.
   uint_fast8_t type;
/* 1->hasBaseChr, 2->hasQchar, 4->hasBase2bit,
   8->hashexBQ or just array of {0,128} for N,
  16->Eamss bit set, 32->Q value not raw, after Eamss. */
} SeqFileObj;   // We may support both FA/FQ and binary formats. So, object is a good thing.

SeqFileObj * inSeqFinit(const char * const, unsigned char);

void inSeqFdestroy(SeqFileObj * const);

#endif /* libgFileIO.h */
