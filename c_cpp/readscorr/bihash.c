#include <stdlib.h>	//calloc
#include <stdint.h>	// uint_fast8_t
#include <math.h>	//log2, ceil
#include <stdio.h>	//fprintf, fseek
#include <err.h>
#include <string.h> //strcmp, strncpy
#include <sys/mman.h>
#include <endian.h> //BYTE_ORDER, LITTLE_ENDIAN 1234

#ifdef PTHREAD
#include <pthread.h>
#endif

#include "bihash.h"
#include "MurmurHash3.h"
#include "2bitseq.h"
#include "chrseq.h"

#define HASHSEED (0x3ab2ae35)

unsigned char GETitemByte_PADrBit_trimSubItemCount(unsigned char CountBit, unsigned char *prBit, uint16_t *pSubItemCount){
    unsigned char itemByte = (CountBit+*prBit+7u) >> 3;	// 2^3=8
    *prBit = (itemByte<<3) - CountBit;
#ifndef TEST    /* SubItemCount is trimed on rBit only */
    double maxSubItemByR = floor(pow(2.0,(double)*prBit));
    if ( *pSubItemCount > maxSubItemByR ) *pSubItemCount = (uint16_t)maxSubItemByR; // safe since *pSubItemCount was bigger
#endif
    return itemByte;
}

// the smarter one
SDLeftArray_t *dleft_arraynew(const unsigned char CountBit, const SDLConfig * const psdlcfg, int kmer){
    unsigned char rBit;
    size_t ArraySize;
    uint16_t SubItemCount;
    ;
    return dleft_arrayinit(CountBit, ArraySize, SubItemCount, kmer);
}

// the native one
SDLeftArray_t *dleft_arrayinit(unsigned char CountBit, size_t ArraySize, uint16_t SubItemCount, int kmer) {
    if (ArraySize<2u || ArraySize>(1LLU<<63) || CountBit<3u || CountBit>8u*sizeof(uint64_t) || SubItemCount<1u ) {
       err(EXIT_FAILURE, "[x]Wrong D Left Array Parameters:(%d)[%u]x%zd ",CountBit,SubItemCount,ArraySize);
    }   // CountBit+rBit >= 9, makes uint16_t always OK
    unsigned char ArrayBit = ceil(log2(ArraySize));
    unsigned char rBit = ArrayBit + ENCODEDBIT_ENTROPY_PAD;
    int toEncode = kmer*2 - rBit;
    if (toEncode > ArrayBit) {
        toEncode = ArrayBit;
        rBit = kmer*2 - ArrayBit;
    }
    unsigned char EncodedBit=(toEncode>0)?toEncode:0;

#ifdef TEST    /* Test mode, keep rBit, pad CountBit */
    unsigned char itemByte = GETitemByte_PADrBit_trimSubItemCount(rBit,&CountBit,&SubItemCount);
#else   /* Normal, keep CountBit, pad rBit */
    unsigned char itemByte = GETitemByte_PADrBit_trimSubItemCount(CountBit,&rBit,&SubItemCount);
#endif
    SDLeftArray_t *dleftobj = calloc(1,sizeof(SDLeftArray_t));    // set other int to 0
    dleftobj->SDLAbyte = (SubItemCount*itemByte*ArraySize+127u)&(~(size_t)127u);    // We are reading in uint128_t now.
    dleftobj->pDLA = calloc(1,dleftobj->SDLAbyte);
    int mlock_r = mlock(dleftobj->pDLA,dleftobj->SDLAbyte);
    if (mlock_r) warn("[!]Cannot lock SDL array in memory. Performance maybe reduced.");
    dleftobj->CountBit = CountBit;
    dleftobj->rBit = rBit;
    dleftobj->ArrayBit = ArrayBit;
    dleftobj->EncodedBit = EncodedBit;
    dleftobj->itemByte = itemByte;
    dleftobj->ArraySize = ArraySize;
    dleftobj->SubItemCount = SubItemCount;
    dleftobj->SubItemByUnit = SubItemCount/SDL_SUBARRAY_UNIT;
    dleftobj->maxCountSeen = 1; //if SDLA is not empty, maxCountSeen>=1.
    dleftobj->FalsePositiveRatio = exp2(-rBit)/(double)ArraySize;
    //dleftobj->ItemInsideAll=0;
    //dleftobj->CellOverflowCount=0;
    dleftobj->Item_rBitMask=(uint128_t)((1LLU<<rBit)-1u) << CountBit;
    dleftobj->Item_CountBitMask=(1LLU<<CountBit)-1u;    // == maxCount
    dleftobj->Hash_ArrayBitMask=(1LLU<<ArrayBit)-1u;
    dleftobj->Hash_rBitMask=(1LLU<<rBit)-1u;
    return dleftobj;
}

void fprintSDLAnfo(FILE *stream, const SDLeftArray_t * dleftobj){
    fprintf(stream,"SDLA(%#zx) -> {\n\
 Size:[r:%uB+cnt:%uB]*Item:%zd(%.3f~%uB[enc:%u])*subArray:%u = %g MiB\n\
 Hash:%u*%uB   ItemByte:%u   MaxCountSeen:%lu%s\n\
 Designed Capacity:%lu   ItemCount:%lu, with Overflow:%lu\n\
 FP:%g, estimated FP item count:%.10g\n\
 Mem:%zu bytes\n\
",
      (size_t)dleftobj,
      dleftobj->rBit,dleftobj->CountBit,dleftobj->ArraySize,log2(dleftobj->ArraySize),
        dleftobj->ArrayBit,dleftobj->EncodedBit,dleftobj->SubItemCount,(double)dleftobj->SubItemCount*dleftobj->itemByte*dleftobj->ArraySize/1048576,
      1,HASH_LENB,dleftobj->itemByte,dleftobj->maxCountSeen,(dleftobj->ItemInsideAll)?"":"(=0, as SDLA is empty)",
      dleftobj->ArraySize*(dleftobj->SubItemCount*3/4),
      dleftobj->ItemInsideAll,dleftobj->CellOverflowCount,
      dleftobj->FalsePositiveRatio,dleftobj->ItemInsideAll*dleftobj->FalsePositiveRatio,
      dleftobj->SDLAbyte);
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


