#include <stdlib.h>	//calloc
#include <stdint.h>	// uint_fast8_t
#include <math.h>	//log2, ceil
#include <stdio.h>	//fprintf
#include <err.h>
#include <string.h> //strcmp
#include "sdleft.h"
#include "MurmurHash3.h"
#include "2bitseq.h"
#include "chrseq.h"

#define	FORCE_INLINE static inline __attribute__((always_inline))

SDLeftArray_t *dleft_arrayinit(unsigned char CountBit, unsigned char rBit, size_t ArraySize) {
    if (ArraySize<2 || CountBit<1 || rBit<1)
       err(EXIT_FAILURE, "[x]Wrong D Left Array Parameters:(%d+%d)x%zd ",rBit,CountBit,ArraySize);
    unsigned char itemByte = (CountBit+rBit+7u) >> 3;	// 2^3=8
    unsigned char ArrayBit = ceil(log2(ArraySize));
    SDLeftArray_t *dleftobj = malloc(sizeof(SDLeftArray_t));
    dleftobj->pDLA = calloc(SDLA_ITEMARRAY,itemByte*ArraySize);
    //unsigned char SDLArray[ArraySize][SDLA_ITEMARRAY][itemByte];
    //the GNU C Compiler allocates memory for VLAs on the stack, >_<
    dleftobj->CountBit = CountBit;
    dleftobj->rBit = (itemByte<<3) - CountBit;
    dleftobj->ArrayBit = ArrayBit;
    dleftobj->itemByte = itemByte;
    dleftobj->ArraySize = ArraySize;
    //dleftobj->ArrayCount = ArrayCount;
    dleftobj->HashCnt = (rBit+ArrayBit+(HASH_LENB-1))/HASH_LENB;
    dleftobj->outhash = malloc(dleftobj->HashCnt * HASH_LENB/8);
    dleftobj->FalsePositiveRatio = (SDLA_ITEMARRAY*3/4)*exp2(-rBit);
    dleftobj->ItemInsideAll=0;
    dleftobj->CellOverflowCount=0;
    return dleftobj;
}

void fprintSDLAnfo(FILE *stream, const SDLeftArray_t * dleftobj){
    fprintf(stream,"SDLA(%#zx) -> {\n\
 Size:[r:%uB+cnt:%uB]*Item:%zd(%.3f~%uB)*subArray:%u\n\
 Hash:%u*%uB   Designed Capacity:%lu\n\
 ItemCount:%lu, with Overflow:%lu\n\
 FP:%g, estimated FP item count:%.2f\n\
",
      (size_t)dleftobj,
      dleftobj->rBit,dleftobj->CountBit,dleftobj->ArraySize,log2(dleftobj->ArraySize),
        dleftobj->ArrayBit,SDLA_ITEMARRAY,
      dleftobj->HashCnt,HASH_LENB, dleftobj->ArraySize*(SDLA_ITEMARRAY*3/4),
      dleftobj->ItemInsideAll,dleftobj->CellOverflowCount,
      dleftobj->FalsePositiveRatio,dleftobj->ItemInsideAll*dleftobj->FalsePositiveRatio);
    fputs("}\n", stream);
}

FORCE_INLINE uint32_t rotl32 ( uint32_t x, int8_t r ){
  return (x << r) | (x >> (32 - r));
}

// this is POP, thus the input *pdat will be changed/shifted.
// Max bits is the same as size_t
// this is a static function, bits has already been <= 8*sizeof(size_t)
FORCE_INLINE size_t popLowestBits(unsigned char bits, uint64_t *pdat, uint_fast8_t *pdatLen64u){
    //unsigned char lastBits = bits % (8*sizeof(size_t));
    unsigned char nullBits = 8*sizeof(size_t)-bits;
    //char MoreUnit = (bits-lastBits)/sizeof(size_t);
    size_t bitMask = (1LLU<<bits)-1U;
    size_t outUnit = pdat[0] & bitMask;
printf("[b]%d,%d,%016zx,%016zx\n",bits,nullBits,bitMask,outUnit);
    uint_fast8_t i;
    for(i=1;i<*pdatLen64u;i++){
        size_t tmpUnit = *(pdat+i) & bitMask;
        *(pdat+i-1) = (*(pdat+i-1) >> bits) | (tmpUnit << nullBits);
    }
    *(pdat+i-1) = *(pdat+i-1) >> bits;
    *pdatLen64u -= bits/(8*sizeof(size_t));
//printf("[c][%016lx,%016lx]\n",pdat[0],pdat[1]);
    return outUnit;
}

FORCE_INLINE int_fast8_t dleft_insert_kmer(const char *const kmer, const size_t len, SDLeftArray_t *dleftobj) {
    char* revcomkmer = ChrSeqRevComp(kmer,len);
char xx = strcmp(kmer,revcomkmer);
    const char *const smallerkmer = (strcmp(kmer,revcomkmer)<=0)?kmer:revcomkmer;
printf("[%zd]->[%s,%s,%s] (%d)\n",len,kmer,revcomkmer,smallerkmer,xx);
    uint64_t *dibskmer=NULL;
    size_t uint64cnt = 0;
    size_t bytelen = (len+3u)/4;
    size_t Ncount = ChrSeq2dib(smallerkmer,len,&dibskmer,&uint64cnt);
printf("%zd:%zd:[%016lx][%016lx]->[%s] (%zd)\n",uint64cnt,bytelen,dibskmer[0],dibskmer[1],dib2basechr(dibskmer,len),Ncount);
    uint64_t *ptmpout = dleftobj->outhash;
    for(uint_fast8_t i=0;i<dleftobj->HashCnt;i++){
        uint32_t seed=rotl32(0x3ab2ae35-i,i);
        MurmurHash3_x64_128(dibskmer,bytelen,seed,ptmpout);
printf("[%016lx,%016lx]\n",ptmpout[0],ptmpout[1]);
        ptmpout += 2;
    }
    uint_fast8_t datLen64u = dleftobj->HashCnt * HASH_LENB/(8*sizeof(uint64_t));
printf("[x][%016lx,%016lx]\n",popLowestBits(60,dleftobj->outhash,&datLen64u),popLowestBits(56,dleftobj->outhash,&datLen64u));
    free(revcomkmer);
}

int_fast8_t dleft_insert_read(char const *const inseq, size_t len, SDLeftArray_t *dleftobj) {
    dleft_insert_kmer(inseq,len,dleftobj);
}

void dleft_arraydestroy(SDLeftArray_t * const dleftobj){
	free(dleftobj->pDLA);
	free(dleftobj->outhash);
	free(dleftobj);
}
