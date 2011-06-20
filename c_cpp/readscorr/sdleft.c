#include <stdlib.h>	//calloc
#include <stdint.h>	// uint_fast8_t
#include <math.h>	//log2, ceil
#include <stdio.h>	//fprintf, fseek
#include <err.h>
#include <string.h> //strcmp, strncpy
#include <sys/mman.h>
//#include <asm/byteorder.h>  // __LITTLE_ENDIAN_BITFIELD or __BIG_ENDIAN_BITFIELD
#include "sdleft.h"
//#include "sdleftTF.h"
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
SDLeftArray_t *dleft_arraynew(const unsigned char CountBit, const SDLConfig * const psdlcfg){
    unsigned char rBit;
    size_t ArraySize;
    uint16_t SubItemCount;
    ;
    return dleft_arrayinit(CountBit, rBit, ArraySize, SubItemCount);
}

// the native one
SDLeftArray_t *dleft_arrayinit(unsigned char CountBit, unsigned char rBit, size_t ArraySize, uint16_t SubItemCount) {
    if (ArraySize<2u || CountBit<1u || rBit<1u || rBit>8u*sizeof(uint64_t) || CountBit>8u*sizeof(uint64_t) || SubItemCount<1u ) {
       err(EXIT_FAILURE, "[x]Wrong D Left Array Parameters:(%d+%d)[%u]x%zd ",rBit,CountBit,SubItemCount,ArraySize);
    }
#ifdef TEST    /* Test mode, keep rBit, pad CountBit */
    unsigned char itemByte = GETitemByte_PADrBit_trimSubItemCount(rBit,&CountBit,&SubItemCount);
#else   /* Normal, keep CountBit, pad rBit */
    unsigned char itemByte = GETitemByte_PADrBit_trimSubItemCount(CountBit,&rBit,&SubItemCount);
#endif
    unsigned char ArrayBit = ceil(log2(ArraySize));
    SDLeftArray_t *dleftobj = calloc(1,sizeof(SDLeftArray_t));    // set other int to 0
    dleftobj->SDLAbyte = SubItemCount*itemByte*ArraySize;
    dleftobj->pDLA = calloc(1,dleftobj->SDLAbyte);
    int mlock_r = mlock(dleftobj->pDLA,dleftobj->SDLAbyte);
    if (mlock_r) warn("[!]Cannot lock SDL array in memory. Performance maybe reduced.");
    //unsigned char SDLArray[ArraySize][SubItemCount][itemByte];
    //the GNU C Compiler allocates memory for VLAs on the stack, >_<
    dleftobj->CountBit = CountBit;
    dleftobj->rBit = rBit;
    dleftobj->ArrayBit = ArrayBit;
    dleftobj->itemByte = itemByte;
    dleftobj->ArraySize = ArraySize;
    dleftobj->SubItemCount = SubItemCount;
    dleftobj->maxCountSeen = 1; //if SDLA is not empty, maxCountSeen>=1.
    /* only one 128bit hash needed.
    dleftobj->ArrayCount = ArrayCount;
    dleftobj->HashCnt = (rBit+ArrayBit+(HASH_LENB-1))/HASH_LENB;
    dleftobj->outhash = malloc(dleftobj->HashCnt * HASH_LENB/8);*/
    dleftobj->FalsePositiveRatio = exp2(-rBit)/(double)ArraySize;
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
 Hash:%u*%uB   ItemByte:%u   MaxCountSeen:%lu%s\n\
 Designed Capacity:%lu   ItemCount:%lu, with Overflow:%lu\n\
 FP:%g, estimated FP item count:%.10g\n\
 Mem:%zu bytes\n\
",
      (size_t)dleftobj,
      dleftobj->rBit,dleftobj->CountBit,dleftobj->ArraySize,log2(dleftobj->ArraySize),
        dleftobj->ArrayBit,dleftobj->SubItemCount,(double)dleftobj->SubItemCount*dleftobj->itemByte*dleftobj->ArraySize/1048576,
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
FORCE_INLINE void incSDLArray(size_t ArrayBits, uint64_t rBits, SDLeftArray_t *dleftobj){
    size_t ArrayPos = ArrayBits % dleftobj->ArraySize;
    size_t relAddr = ArrayPos*dleftobj->SubItemCount*dleftobj->itemByte;
    unsigned char* pChunk = (unsigned char*)dleftobj->pDLA + relAddr;
    unsigned char* pEndChunk = (unsigned char*)pChunk + dleftobj->SubItemCount*dleftobj->itemByte;
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
                if (Item_CountBits > dleftobj->maxCountSeen) {
                    dleftobj->maxCountSeen = Item_CountBits;
#ifdef DEBUG
fprintf(stderr,"[sdlm][%zu][%lu]:[%lx],[%lu] [%016lx %016lx]\n",
    ArrayPos,dleftobj->SubItemCount-(pEndChunk-pChunk)/dleftobj->itemByte,rBits,Item_CountBits,(uint64_t)(theItem>>64u),(uint64_t)theItem);
#endif
                }   // if SDLA is not empty, maxCountSeen>=1, no need to check when Item_CountBits == 0.
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
FORCE_INLINE uint64_t querySDLArray(size_t ArrayBits, uint64_t rBits, SDLeftArray_t *dleftobj){
    size_t ArrayPos = ArrayBits % dleftobj->ArraySize;
    size_t relAddr = ArrayPos*dleftobj->SubItemCount*dleftobj->itemByte;
    unsigned char* pChunk = (unsigned char*)dleftobj->pDLA + relAddr;
    unsigned char* pEndChunk = (unsigned char*)pChunk + dleftobj->SubItemCount*dleftobj->itemByte;
    uint64_t Item_rBits;
    uint64_t Item_CountBits=0;  // set value in case SubItemCount*dleftobj->itemByte == 0 (EP ?)
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
        uint64_t rBits = popLowestBits(dleftobj->rBit,dleftobj->outhash,&datLenu64t);
        size_t ArrayBits = popLowestBits(dleftobj->ArrayBit,dleftobj->outhash,&datLenu64t);
        incSDLArray(ArrayBits, rBits, dleftobj);
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

void dleft_dump(const SDLeftArray_t * const dleftobj, SDLdumpHead * const pSDLeftStat, FILE *stream){
    rewind(stream);  // for Binary file, position is important.
    strncpy(pSDLeftStat->FileID,"GDSD",4);
    pSDLeftStat->FileVersion[0]=0u;
    pSDLeftStat->FileVersion[1]=1u;
    //pSDLeftStat->extreebyte = 0;
    pSDLeftStat->SubItemCount=dleftobj->SubItemCount;
    pSDLeftStat->CountBit=dleftobj->CountBit;
    pSDLeftStat->rBit=dleftobj->rBit;
    pSDLeftStat->ArraySize=dleftobj->ArraySize;
    pSDLeftStat->SDLAbyte=dleftobj->SDLAbyte;
    pSDLeftStat->ItemInsideAll=dleftobj->ItemInsideAll;
    pSDLeftStat->CellOverflowCount=dleftobj->CellOverflowCount;
    pSDLeftStat->CountBitOverflow=dleftobj->CountBitOverflow;
    pSDLeftStat->maxCountSeen=dleftobj->maxCountSeen;
    pSDLeftStat->crc32c=0xffffffff; //later
    size_t unitwritten;
    unitwritten=fwrite(pSDLeftStat,sizeof(SDLdumpHead),1u,stream);
    if (unitwritten != 1)
        err(EXIT_FAILURE, "Fail to write dat file ! [%zd]",unitwritten);
    unitwritten=fwrite(dleftobj->pDLA,pSDLeftStat->SDLAbyte,1u,stream);
    if (unitwritten != 1)
        err(EXIT_FAILURE, "Cannot write to dat file ! [%zd]",unitwritten);
printf("Cb %x rB %x AS %lx Size %lx HMC %lx HMH %lx\n",pSDLeftStat->CountBit,pSDLeftStat->rBit,pSDLeftStat->ArraySize,pSDLeftStat->SDLAbyte,pSDLeftStat->HistMaxCntVal,pSDLeftStat->HistMaxHistVal);
}

void dleft_arraydestroy(SDLeftArray_t * const dleftobj){
	free(dleftobj->pDLA);
	//free(dleftobj->outhash);
	free(dleftobj);
}
