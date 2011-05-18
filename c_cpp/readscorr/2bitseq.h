#ifndef _G_DBSEQ_H
#define _G_DBSEQ_H

#include <stdint.h>	// int_fast8_t
#include <stdlib.h> // malloc
#include "2bitseqinline.h"

#ifndef DBIT2BASES
#define DBIT2BASES
#define DBIT2BASE(dbit) ("acgt"[dbit]) // ("ACGT"[dbit])
#define DBIT2COMPBASE(dbit) ("TGCA"[dbit])
#endif

char *unit2basechr(uint64_t);
char *dib2basechr(uint64_t *, size_t);

uint64_t *dibrevcomp(uint64_t const *const, size_t);
//uint64_t *dibmalloc(size_t len);

#endif  // 2bitseq.h

