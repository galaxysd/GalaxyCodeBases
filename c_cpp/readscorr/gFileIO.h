// by Hu Xuesong
#ifndef _G_FILEIO_H
#define _G_FILEIO_H

#include <stdio.h>	// for fpos_t

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;	// l for sequence length in base-pairs, m for malloc in bytes.
	char *s;	// for 2bit ATCG, 1bp=2bit; for 6bit quality, 1bp=6bit.
} kstring_t;
#endif	// Yes, it is the same as that in kseq.h

typedef struct __SeqFileObj {
   void * fh;   // Not just FILE *fp
   int (*getNextSeq)(void *);   // void * fh
   fpos_t datePos[1];   // [1] for item id, [0] for offset. For Text file, item id == 0.
   kstring_t * name, comment;   // when SeqFileObj.ReadNext, they are being updated.
   kstring_t diBseq;   // 2bit, {A,C,G,T}={0,1,2,3}
   kstring_t hexBQ;   // 0~63 for Quality, 128 for N, 64 for Eamss-masked
} SeqFileObj;   // We may support both FA/FQ and binary formats. So, object is a good thing.


#endif /* gFileIO.h */

