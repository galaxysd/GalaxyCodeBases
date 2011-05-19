#include <stdlib.h>	//calloc
#include <stdint.h>	// uint_fast8_t
#include <math.h>	//log2, ceil
#include <stdio.h>	//fprintf
#include <err.h>
#include <string.h> //strcmp
#include "dleft.h"
#include "MurmurHash3.h"
#include "2bitseq.h"
#include "chrseq.h"

#define HASH_LENB 128u
#define DLA_ITEMARRAY 8u

#define	FORCE_INLINE static inline __attribute__((always_inline))

DLeftArray_t *dleft_arrayinit(unsigned char CountBit, unsigned char rBit, size_t ArraySize, unsigned char ArrayCount) {
    if (ArraySize<2 || ArrayCount<1 || CountBit<1 || rBit<1)
       err(EXIT_FAILURE, "[x]Wrong D Left Array Parameters:(%d+%d)x%zdx%d",rBit,CountBit,ArraySize,ArrayCount);
    unsigned char itemByte = (CountBit+rBit+7u) >> 3;	// 2^3=8
    unsigned char ArrayBit = ceil(log2(ArraySize/ArrayCount));
    DLeftArray_t *dleftobj = malloc(sizeof(DLeftArray_t));
    dleftobj->dlap = calloc(ArrayCount,itemByte*ArraySize*DLA_ITEMARRAY);
    dleftobj->CountBit = CountBit;
    dleftobj->rBit = (itemByte<<3) - CountBit;
    dleftobj->ArrayBit = ArrayBit;
    dleftobj->ArraySize = ArraySize;
    dleftobj->ArrayCount = ArrayCount;
    dleftobj->HashCnt = (rBit+ArrayBit+HASH_LENB-1)/HASH_LENB;
    dleftobj->outhash = malloc(dleftobj->HashCnt * HASH_LENB/8);
    dleftobj->FalsePositiveRatio = ArrayCount*(DLA_ITEMARRAY*3/4)*exp2(-rBit);
    dleftobj->ItemInsideAll=0;
    dleftobj->CellOverflowCount=0;
    return dleftobj;
}

void fprintDLAnfo(FILE *stream, const DLeftArray_t * dleftobj){
    fprintf(stream,"%#zx -> {\n\
 Size:[r:%uB+cnt:%uB]*Item:%zd(%.3f~%uB)*Array:%u\n\
 Hash:%u*%uB   Designed Capacity:%lu\n\
 ItemCount:%lu, with Overflow:%lu\n\
 FP:%g, estimated FP item count:%.2f\n\
",
      (size_t)dleftobj,
      dleftobj->rBit,dleftobj->CountBit,dleftobj->ArraySize,log2(dleftobj->ArraySize/dleftobj->ArrayCount),
        dleftobj->ArrayBit,dleftobj->ArrayCount,
      dleftobj->HashCnt,HASH_LENB, dleftobj->ArrayCount*dleftobj->ArraySize*(DLA_ITEMARRAY*3/4),
      dleftobj->ItemInsideAll,dleftobj->CellOverflowCount,
      dleftobj->FalsePositiveRatio,dleftobj->ItemInsideAll*dleftobj->FalsePositiveRatio);
    fputs("}\n", stream);
}

FORCE_INLINE uint32_t rotl32 ( uint32_t x, int8_t r ){
  return (x << r) | (x >> (32 - r));
}

FORCE_INLINE int_fast8_t dleft_insert_kmer(const char *const kmer, const size_t len, DLeftArray_t *dleftobj) {
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
        uint32_t seed=0x3ab2ae35-i;
        MurmurHash3_x64_128(dibskmer,bytelen,seed,ptmpout);
printf("[%016lx,%016lx]\n",ptmpout[0],ptmpout[1]);
        ptmpout += 2;
    }
    free(revcomkmer);
}

int_fast8_t dleft_insert_read(char const *const inseq, size_t len, DLeftArray_t *dleftobj) {
    dleft_insert_kmer(inseq,len,dleftobj);
}

void dleft_arraydestroy(DLeftArray_t * const dleftobj){
	free(dleftobj->dlap);
	free(dleftobj->outhash);
	free(dleftobj);
}
