
#include <stdio.h>
#include <stdlib.h>
#include "MiscUtilities.h"
#include "MemManager.h"
#include "TextConverter.h"
#include "BWT.h"
#include "HSP.h"
#include "Types.h"

unsigned int REVBWTForwardSearch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *rev_bwt,unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight,unsigned int *rev_resultSaIndexLeft, unsigned int *rev_resultSaIndexRight);

unsigned int REVBWTContForwardSearch(const unsigned char *convertedKey, const unsigned int start, const unsigned int len,const BWT *rev_bwt,unsigned int *saL, unsigned int *saR,unsigned int *rev_saL, unsigned int *rev_saR);

unsigned int BWTContBackwardSearch(const unsigned char *convertedKey, const unsigned int start, const unsigned int len, const BWT *bwt, unsigned int *saL, unsigned int *saR);
unsigned int BWTBackward1Error(char *querypattern, int chain, BWT *bwt, unsigned int start, unsigned int len, unsigned int pl, unsigned int pr, unsigned int allele2, HitInfo *hits, unsigned int *numOfHits);
unsigned int REVBWTForward1Error(char *queryPattern,int chain, BWT *bwt, BWT * rev_bwt, unsigned int start,unsigned int len, unsigned int pl,unsigned int pr, unsigned int rev_pl,unsigned int rev_pr, unsigned int allele2, HitInfo *hits, unsigned int *numOfHits);
