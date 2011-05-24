#include <stdlib.h>	//calloc
#include <stdint.h>	// uint_fast8_t
#include <math.h>	//log2, ceil
#include <stdio.h>	//fprintf
#include <err.h>
#include <string.h> //strcmp
//#include <asm/byteorder.h>  // __LITTLE_ENDIAN_BITFIELD or __BIG_ENDIAN_BITFIELD
#include "sdleft.h"
#include "MurmurHash3.h"
#include "2bitseq.h"
#include "chrseq.h"

#define	FORCE_INLINE static inline __attribute__((always_inline))

#define HASHSEED 0x3ab2ae35

SDLeftArray_t *dleft_arrayinit(unsigned char CountBit, unsigned char rBit, size_t ArraySize) {
    if (ArraySize<2 || CountBit<1 || rBit<1 || rBit>8*sizeof(uint64_t) || CountBit>8*sizeof(uint64_t) )
       err(EXIT_FAILURE, "[x]Wrong D Left Array Parameters:(%d+%d)x%zd ",rBit,CountBit,ArraySize);
    unsigned char itemByte = (CountBit+rBit+7u) >> 3;	// 2^3=8
    unsigned char ArrayBit = ceil(log2(ArraySize));
    SDLeftArray_t *dleftobj = calloc(1,sizeof(SDLeftArray_t));    // set other int to 0
    dleftobj->pDLA = calloc(SDLA_ITEMARRAY,itemByte*ArraySize);
    //unsigned char SDLArray[ArraySize][SDLA_ITEMARRAY][itemByte];
    //the GNU C Compiler allocates memory for VLAs on the stack, >_<
    dleftobj->CountBit = CountBit;
    //dleftobj->maxCount = (1LLU<<CountBit)-iu;
    rBit = (itemByte<<3) - CountBit;
    dleftobj->rBit = rBit;
    dleftobj->ArrayBit = ArrayBit;
    dleftobj->itemByte = itemByte;
    dleftobj->ArraySize = ArraySize;
    //dleftobj->ArrayCount = ArrayCount;
    //dleftobj->HashCnt = (rBit+ArrayBit+(HASH_LENB-1))/HASH_LENB;
    //dleftobj->outhash = malloc(dleftobj->HashCnt * HASH_LENB/8);
    dleftobj->FalsePositiveRatio = (SDLA_ITEMARRAY*3/4)*exp2(-rBit);
    dleftobj->ItemInsideAll=0;
    dleftobj->CellOverflowCount=0;
    dleftobj->Item_rBitMask=(uint128_t)((1LLU<<rBit)-1u) << CountBit;
    dleftobj->Item_CountBitMask=(1LLU<<CountBit)-1u;    // == maxCount
    dleftobj->Hash_ArrayBitMask=(1LLU<<ArrayBit)-1u;
    dleftobj->Hash_rBitMask=(1LLU<<rBit)-1u;
    return dleftobj;
}

void fprintSDLAnfo(FILE *stream, const SDLeftArray_t * dleftobj){
    fprintf(stream,"SDLA(%#zx) -> {\n\
 Size:[r:%uB+cnt:%uB]*Item:%zd(%.3f~%uB)*subArray:%u = %g MiB\n\
 Hash:%u*%uB   ItemByte:%u\n\
 Designed Capacity:%lu   ItemCount:%lu, with Overflow:%lu\n\
 FP:%g, estimated FP item count:%.2f\n\
",
      (size_t)dleftobj,
      dleftobj->rBit,dleftobj->CountBit,dleftobj->ArraySize,log2(dleftobj->ArraySize),
        dleftobj->ArrayBit,SDLA_ITEMARRAY,(double)SDLA_ITEMARRAY*dleftobj->itemByte*dleftobj->ArraySize/1048576,
      1,HASH_LENB,dleftobj->itemByte,
      dleftobj->ArraySize*(SDLA_ITEMARRAY*3/4),
      dleftobj->ItemInsideAll,dleftobj->CellOverflowCount,
      dleftobj->FalsePositiveRatio,dleftobj->ItemInsideAll*dleftobj->FalsePositiveRatio);
/*
    fprintf(stream," Item_rBitMask:[%016lX %016lX]\n", (uint64_t)(dleftobj->Item_rBitMask>>64), (uint64_t)dleftobj->Item_rBitMask);
    fprintf(stream," Item_CountBitMask:[%016lX]\n", dleftobj->Item_CountBitMask);
    uint128_t check = dleftobj->Item_rBitMask + dleftobj->Item_CountBitMask;
    fprintf(stream," Sum:[%016lX %016lX]\n", (uint64_t)(check>>64), (uint64_t)check);
    fprintf(stream," Hash_ArrayBitMask:[%016lX]\n", dleftobj->Hash_ArrayBitMask);
    fprintf(stream," Hash_rBitMask:[%016lX]\n", dleftobj->Hash_rBitMask);
*/
    fputs("}\n", stream);
}

FORCE_INLINE uint32_t rotl32 ( uint32_t x, int8_t r ){
  return (x << r) | (x >> (32 - r));
}

// this is POP, thus the input *pdat will be changed/shifted.
// Max bits is the same as size_t
// this is a static function, bits has already been <= 8*sizeof(size_t)
// static function, thus the input parameters should be clean
FORCE_INLINE uint64_t popLowestBits(unsigned char bits, uint64_t *pdat, uint_fast8_t *pdatLenu64t){
    //unsigned char lastBits = bits % (8*sizeof(size_t));
    unsigned char nullBits = 8*sizeof(uint64_t)-bits;
    //char MoreUnit = (bits-lastBits)/sizeof(size_t);
    uint64_t bitMask = (1LLU<<bits)-1U;
    uint64_t outUnit = pdat[0] & bitMask;
//printf("[b]%d,%d,%016zx,%016zx\n",bits,nullBits,bitMask,outUnit);
    uint_fast8_t i;
    for(i=1;i<*pdatLenu64t;i++){
        uint64_t tmpUnit = *(pdat+i) & bitMask;
        *(pdat+i-1) = (*(pdat+i-1) >> bits) | (tmpUnit << nullBits);
    }
    *(pdat+i-1) = *(pdat+i-1) >> bits;
    *pdatLenu64t -= bits/(8*sizeof(uint64_t));
//printf("[c][%016lx,%016lx]\n",pdat[0],pdat[1]);
    return outUnit;
}

// rBits is (0,64]
FORCE_INLINE void incSDLArray(size_t ArrayPos, uint64_t rBits, SDLeftArray_t *dleftobj){
    size_t relAddr = ArrayPos*SDLA_ITEMARRAY*dleftobj->itemByte;
    unsigned char* pChunk = (unsigned char*)dleftobj->pDLA + relAddr;
    unsigned char* pEndChunk = (unsigned char*)pChunk + SDLA_ITEMARRAY*dleftobj->itemByte;
    //unsigned char BoundaryByteRel = rBits % 8;    // Well, I will not touch the egg-head ...
    uint64_t Item_rBits, Item_CountBits;
    uint128_t theItem;
/*
    union u128t {              
        uint128_t all;         
        unsigned char byte[16];
    } theItem;                 
*/
    while (pChunk < pEndChunk) {
        theItem = 0;
        for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
            theItem |= ((uint128_t)*(pChunk+i)) << (i*8u);
            //theItem.byte[i] = *(pChunk+i);
        }
        Item_CountBits = theItem & dleftobj->Item_CountBitMask;
        if (Item_CountBits == 0) {  // reaching the pre-end
            Item_CountBits = 1;
            ++dleftobj->ItemInsideAll;
//printf("<%lu>",dleftobj->ItemInsideAll);
            break;
        }
        Item_rBits = (theItem & dleftobj->Item_rBitMask) >> dleftobj->CountBit;
        if (Item_rBits == rBits) {
            if (Item_CountBits < dleftobj->Item_CountBitMask) {
                ++Item_CountBits;
                break;
            } else {
                ++dleftobj->CountBitOverflow;
                return;
            }
        }
        pChunk += dleftobj->itemByte;
    }  // reaching the structure-end
    if (pChunk < pEndChunk) {
//printf("Old:%zu[%lx %lx]<-[%lx][%lx][%lx %lx] ",(size_t)((char*)pChunk - (char*)dleftobj->pDLA)-relAddr,(uint64_t)(theItem>>64),(uint64_t)theItem,rBits,Item_CountBits,(uint64_t)((((uint128_t)rBits)<<dleftobj->CountBit)>>64),(uint64_t)(((uint128_t)rBits)<<dleftobj->CountBit));
        theItem = (((uint128_t)rBits)<<dleftobj->CountBit) | Item_CountBits;
        for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
            uint128_t tmpMask = ((uint128_t)0xffLLU) << (i*8u);
            *pChunk++ = (theItem & tmpMask)>>(i*8u);
        }
/*printf("New:%zu[%lx %lx] ",(size_t)((char*)pChunk - (char*)dleftobj->pDLA)-relAddr,(uint64_t)(theItem>>64),(uint64_t)theItem);
printf("Mem:%zu[",(size_t)((char*)pChunk - (char*)dleftobj->pDLA)-relAddr);
for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
    printf("%x ",*(pChunk-dleftobj->itemByte+i));
}
puts("]");*/
    } else {  // reaching the structure-end
        ++dleftobj->CellOverflowCount;
        // deal(*pextree);
    }
}
// return 0 for not found
FORCE_INLINE uint64_t querySDLArray(size_t ArrayPos, uint64_t rBits, SDLeftArray_t *dleftobj){
    size_t relAddr = ArrayPos*SDLA_ITEMARRAY*dleftobj->itemByte;
    unsigned char* pChunk = (unsigned char*)dleftobj->pDLA + relAddr;
    unsigned char* pEndChunk = (unsigned char*)pChunk + SDLA_ITEMARRAY*dleftobj->itemByte;
    uint64_t Item_rBits;
    uint64_t Item_CountBits=0;  // set value in case SDLA_ITEMARRAY*dleftobj->itemByte == 0 (EP ?)
    uint128_t theItem;
    while (pChunk < pEndChunk) {
        theItem = 0;
        for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
            theItem |= ((uint128_t)*(pChunk+i)) << (i*8u);
        }
        Item_CountBits = theItem & dleftobj->Item_CountBitMask;
        if (Item_CountBits == 0) {  // reaching the pre-end, not found
            break;
        }
        Item_rBits = (theItem & dleftobj->Item_rBitMask) >> dleftobj->CountBit;
        if (Item_rBits == rBits) {  // found
            break;
        }
        pChunk += dleftobj->itemByte;
    }
    if (pChunk >= pEndChunk) {  // reaching the structure-end
        // deal(*pextree);
    }
    return Item_CountBits;
}

FORCE_INLINE int_fast8_t dleft_insert_kmer(const char *const kmer, const size_t len, SDLeftArray_t *dleftobj,
                                           uint64_t **dibskmer,size_t * const uint64cnt) {
    char* revcomkmer = ChrSeqRevComp(kmer,len);
//char xx = strcmp(kmer,revcomkmer);
    const char *const smallerkmer = (strcmp(kmer,revcomkmer)<=0)?kmer:revcomkmer;   // not strncmp since the first odd len bytes mush be different.
//printf("[%zd]->[%s,%s,%s] (%d)\n",len,kmer,revcomkmer,smallerkmer,xx);
    size_t bytelen = (len+3u)/4;
    size_t Ncount = ChrSeq2dib(smallerkmer,len,dibskmer,uint64cnt);
//printf("%zd:%zd:[%016lx][%016lx]->[%s] (%zd)\n",*uint64cnt,bytelen,*dibskmer[0],*dibskmer[1],dib2basechr(*dibskmer,len),Ncount);
    //uint64_t *ptmpout = dleftobj->outhash;
    //for(uint_fast8_t i=0;i<dleftobj->HashCnt;i++){
        //uint32_t seed=rotl32(0x3ab2ae35-i,i);
        //MurmurHash3_x64_128(*dibskmer,bytelen,seed,ptmpout);
    free(revcomkmer);
    if (! Ncount) {
        MurmurHash3_x64_128(*dibskmer,bytelen,HASHSEED,dleftobj->outhash);
//printf("[%016lx,%016lx]\n",dleftobj->outhash[0],dleftobj->outhash[1]);
            //ptmpout += 2;
        //}
        uint_fast8_t datLenu64t = HASH_LENB/(8*sizeof(uint64_t));
//printf("[x][%016lx,%016lx]\n",popLowestBits(60,dleftobj->outhash,&datLenu64t),popLowestBits(56,dleftobj->outhash,&datLenu64t));
        size_t ArrayPos = popLowestBits(dleftobj->ArrayBit,dleftobj->outhash,&datLenu64t) % dleftobj->ArraySize;
        uint64_t rBits = popLowestBits(dleftobj->rBit,dleftobj->outhash,&datLenu64t);
        incSDLArray(ArrayPos, rBits, dleftobj);
        return 1;
    }
    return 0;
}

size_t dleft_insert_read(unsigned int k, char const *const inseq, size_t len, SDLeftArray_t *dleftobj) {
    if (len<k) return 0;
    size_t insertedCount=0;
    uint64_t *dibskmer=NULL;
    size_t uint64cnt = 0;
    for (size_t i=0;i<=len-k;i++) {
        if (dleft_insert_kmer(inseq+i,k,dleftobj,&dibskmer,&uint64cnt)) {
            ++insertedCount;
        }
    }
    free(dibskmer);
    return insertedCount;
}

void dleft_arraydestroy(SDLeftArray_t * const dleftobj){
	free(dleftobj->pDLA);
	//free(dleftobj->outhash);
	free(dleftobj);
}
