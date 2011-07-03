#include <stdlib.h>	//calloc
#include <stdint.h>	// uint_fast8_t
#include <math.h>	//log2, ceil
#include <stdio.h>	//fprintf, fseek
#include <err.h>
#include <string.h> //strcmp, strncpy
#include <sys/mman.h>
#include <endian.h> //BYTE_ORDER, LITTLE_ENDIAN 1234
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
    if (ArraySize<2u || CountBit<3u || rBit<6u || rBit>8u*sizeof(uint64_t) || CountBit>8u*sizeof(uint64_t) || SubItemCount<1u ) {
       err(EXIT_FAILURE, "[x]Wrong D Left Array Parameters:(%d+%d)[%u]x%zd ",rBit,CountBit,SubItemCount,ArraySize);
    }   // CountBit+rBit >= 9, makes uint16_t always OK
#ifdef TEST    /* Test mode, keep rBit, pad CountBit */
    unsigned char itemByte = GETitemByte_PADrBit_trimSubItemCount(rBit,&CountBit,&SubItemCount);
#else   /* Normal, keep CountBit, pad rBit */
    unsigned char itemByte = GETitemByte_PADrBit_trimSubItemCount(CountBit,&rBit,&SubItemCount);
#endif
    unsigned char ArrayBit = ceil(log2(ArraySize));
    SDLeftArray_t *dleftobj = calloc(1,sizeof(SDLeftArray_t));    // set other int to 0
    dleftobj->SDLAbyte = (SubItemCount*itemByte*ArraySize+127u)&(~(size_t)127u);    // We are reading in uint128_t now.
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
    //size_t ArrayPos = ArrayBits % dleftobj->ArraySize;
    size_t ArrayPos = ((uint128_t)ArrayBits*(uint128_t)dleftobj->ArraySize) >> dleftobj->ArrayBit;
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
#ifdef OLD
        theItem = 0;
        for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
            theItem |= ((uint128_t)*(pChunk+i)) << (i*8u);
            //theItem.byte[i] = *(pChunk+i);
        }
#elif defined NEW  /* faster for one less register shift operation for memory uint8_t */
        theItem = 0;
        for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
            theItem = theItem << 8u;
            theItem |= *(pChunk+i);
            //theItem.byte[i] = *(pChunk+i);
        }
#elif BYTE_ORDER == LITTLE_ENDIAN
        theItem = *(uint128_t*)pChunk;
#else
    #error Faster version is Little Endian, choose OLD or NEW to define !
#endif
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
/*
        theItem &= ~dleftobj->CountBit;
        theItem |= Item_CountBits;
        *(uint128_t*)pChunk = theItem;
*/
        theItem = (((uint128_t)rBits)<<dleftobj->CountBit) | Item_CountBits;
#ifdef NEW
//printf("Old: %lx, %lx\n",(uint64_t)(theItem>>64u),(uint64_t)theItem);
        pChunk += dleftobj->itemByte-1;
        for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
            *pChunk-- = (uint8_t)(theItem>>(i*8u));
        }
/*
        ++pChunk;
        theItem = 0;
        for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
            theItem = theItem << 8u;
            theItem |= (uint128_t)*(pChunk+i);
            //theItem.byte[i] = *(pChunk+i);
        }
printf("New: %lx, %lx\n",(uint64_t)(theItem>>64u),(uint64_t)theItem);
*/
#elif (BYTE_ORDER == LITTLE_ENDIAN) || (defined OLD)
        for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
            uint128_t tmpMask = ((uint128_t)0xffLLU) << (i*8u);
            *pChunk++ = (theItem & tmpMask)>>(i*8u);
        }
#else
    #error Faster version is Little Endian, choose OLD or NEW to define !
#endif
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
    size_t ArrayPos = ((uint128_t)ArrayBits*(uint128_t)dleftobj->ArraySize) >> dleftobj->ArrayBit;
    size_t relAddr = ArrayPos*dleftobj->SubItemCount*dleftobj->itemByte;
    unsigned char* pChunk = (unsigned char*)dleftobj->pDLA + relAddr;
    unsigned char* pEndChunk = (unsigned char*)pChunk + dleftobj->SubItemCount*dleftobj->itemByte;
    uint64_t Item_rBits;
    uint64_t Item_CountBits=0;  // set value in case SubItemCount*dleftobj->itemByte == 0 (EP ?)
    uint128_t theItem;
    while (pChunk < pEndChunk) {
#ifdef OLD
        theItem = 0;
        for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
            theItem |= ((uint128_t)*(pChunk+i)) << (i*8u);
        }
#elif defined NEW  /* faster for one less register shift operation for memory uint8_t */
        theItem = 0;
        for (uint_fast8_t i=0;i<dleftobj->itemByte;i++) {
            theItem = theItem << 8u;
            theItem |= *(pChunk+i);
            //theItem.byte[i] = *(pChunk+i);
        }
#elif BYTE_ORDER == LITTLE_ENDIAN
        theItem = *(uint128_t*)pChunk;
#else
    #error Faster version is Little Endian, choose OLD or NEW to define !
#endif
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
    const char *const smallerkmer = (strcmp(kmer,revcomkmer)<=0)?kmer:revcomkmer;   // not strncmp since the first odd len bytes mush be different.
#ifdef DEBUGMORE
 printf("[%zd]->[%s,%s,%s]\n",len,kmer,revcomkmer,smallerkmer);
#endif
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

#ifdef TEST
SDLeftStat_t * dleft_stat(SDLeftArray_t * const dleftobj, FILE *stream, FILE *fpdat) {
#else
SDLeftStat_t * dleft_stat(SDLeftArray_t * const dleftobj, FILE *stream) {
#endif
    SDLeftStat_t *pSDLeftStat = malloc(sizeof(SDLeftStat_t));
    uint64_t * const pCountHistArray = calloc(sizeof(uint64_t),1+dleftobj->maxCountSeen);
    size_t totalDLAsize = dleftobj->SubItemCount * dleftobj->itemByte * dleftobj->ArraySize;
    const unsigned char * const pDLA = dleftobj->pDLA;
    uint64_t Item_CountBits=0;  // set value in case SDLA_ITEMARRAY*dleftobj->itemByte == 0 (EP ?)
    uint128_t theItem;
#ifdef TEST
    uint16_t SubItemCount = dleftobj->SubItemCount;
    uint64_t * const pCountSubArray = calloc(sizeof(uint64_t),1+SubItemCount);
    uint64_t * const pCountSubNoOneArray = calloc(sizeof(uint64_t),1+SubItemCount);
    uint16_t SubItemUsed=0;
    uint16_t SubItemUsedCountIsOne=0;
    uint64_t ArraySize = dleftobj->ArraySize;
    size_t mainArrayID=0;
    rewind(fpdat);
    size_t unitwritten;
    unitwritten=fwrite(&ArraySize,sizeof(uint64_t),1u,fpdat);
    if (unitwritten != 1)
        err(EXIT_FAILURE, "Fail to write dat file ! [%zd]",unitwritten);
#endif
    for (size_t i=0;i<totalDLAsize;i+=dleftobj->itemByte) {
        //const unsigned char * pChunk = pDLA + i;
#ifdef OLD
        theItem = 0;
        for (uint_fast8_t j=0;j<dleftobj->itemByte;j++) {
            theItem |= ((uint128_t)*(pDLA + i + j)) << (j*8u);
        }
#elif defined NEW
        theItem = 0;
        for (uint_fast8_t j=0;j<dleftobj->itemByte;j++) {
            theItem = theItem << 8u;
            theItem |= *(pDLA + i + j);
        }
#elif BYTE_ORDER == LITTLE_ENDIAN
        theItem = *(uint128_t*)(pDLA+i);
#else
    #error Faster version is Little Endian, choose OLD or NEW to define !
#endif
        Item_CountBits = (uint64_t)theItem & dleftobj->Item_CountBitMask;   // Item_CountBitMask is uint64_t.
        ++pCountHistArray[Item_CountBits];
#ifdef TEST
        ++mainArrayID;
        if (Item_CountBits) {
            ++SubItemUsed;
            if (Item_CountBits==1) ++SubItemUsedCountIsOne;
        }
        if (!(mainArrayID % SubItemCount)) {
            ++pCountSubArray[SubItemUsed];
            ++pCountSubNoOneArray[SubItemUsed-SubItemUsedCountIsOne];
            fwrite(&SubItemUsed,sizeof(uint8_t),1u,fpdat);  // Well, just write 1 byte each. (LE)
//printf("-%zd %u\n",mainArrayID,SubItemUsed);
            SubItemUsed=0;
            SubItemUsedCountIsOne=0;
        } else {    // mainArrayID mod SubItemCount != 0. mainArrayID == (i+1)/dleftobj->itemByte .
//printf("---%zd %% %d = %zd %d\n",mainArrayID,SubItemCount,mainArrayID % SubItemCount,SubItemUsed);
        }
#endif
    }
    //THETYPE HistSum=0;    // HistSum == dleftobj->ItemInsideAll
    float128 HistSumSquared=0.0;
    double SStd;    // We need to return in a fixed type for printf
    for (size_t p=1;p<=dleftobj->maxCountSeen;p++) {
        //HistSum += pCountHistArray[p];
        HistSumSquared += pCountHistArray[p] * pCountHistArray[p];
    }
    //http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    SStd = sqrtl( ( HistSumSquared-((long double)dleftobj->ItemInsideAll*(long double)dleftobj->ItemInsideAll/(long double)dleftobj->maxCountSeen) ) / (long double)(dleftobj->maxCountSeen -1) );
    pSDLeftStat->HistSStd = SStd;
    pSDLeftStat->HistMean = (double)dleftobj->ItemInsideAll / (double)dleftobj->maxCountSeen;
    pSDLeftStat->HistMaxCntVal = 1; //later
    pSDLeftStat->HistMaxHistVal = 1; //later
    fprintf(stream,"#Kmer_real_count: %ld\n#Kmer_count_hist: %ld\n"
        "#Kmer_depth_mean: %f\n#Kmer_depth_sStd: %f\n"
        "#CountBit_overflow: %lu\n"
        "\n#Kmer_frequence\tHist_value\tKmer_count\tHist_ratio\n",
        dleftobj->ItemInsideAll,dleftobj->maxCountSeen,pSDLeftStat->HistMean,SStd,dleftobj->CountBitOverflow);
    for (size_t p=1;p<=dleftobj->maxCountSeen;p++) {
        fprintf(stream,"%zu\t%lu\t%lu\t%g\n",p,(uint64_t)pCountHistArray[p],
            (uint64_t)pCountHistArray[p]*(uint64_t)p,(double)pCountHistArray[p]/(double)dleftobj->ItemInsideAll);
    }
    free(pCountHistArray);
    // deal(*pextree);
#ifdef TEST
    fprintf(stream,"#-----------------------------------------------\n\n#Array_size: %lu\n\n"
        "SubArrayFilled\tCount\tCount_without_1\n",ArraySize);
    for (uint16_t i=0;i<=SubItemCount;i++) {
        fprintf(stream,"%u\t%lu, %g\t%lu, %g\n",i,pCountSubArray[i],(double)pCountSubArray[i]/(double)ArraySize,pCountSubNoOneArray[i],(double)pCountSubNoOneArray[i]/(double)ArraySize);
    }
    free(pCountSubArray);
#endif
    return pSDLeftStat;
}

/*
Speed test with OLD 2bitseqinline.h:
./readscorr Saccharomyces_cerevisiae.cfg -o ttmp 2> ttmp.log

==> tnew.log <==
   User: 170.445088(s), System: 0.400939(s). Real: 174.093140(s).
   Sleep: 3.247113(s). Block I/O times: 0/0. MaxRSS: 0 kiB.
   Wait(s): 62(nvcsw) + 54770(nivcsw). Page Fault(s): 41294(minflt) + 0(majflt).

==> told.log <==
   User: 204.444919(s), System: 2.459626(s). Real: 215.708506(s).
   Sleep: 8.803961(s). Block I/O times: 0/0. MaxRSS: 0 kiB.
   Wait(s): 143(nvcsw) + 58522(nivcsw). Page Fault(s): 41294(minflt) + 0(majflt).

==> ttmp.log <==
   User: 133.809657(s), System: 1.160823(s). Real: 137.159722(s).
   Sleep: 2.189242(s). Block I/O times: 0/0. MaxRSS: 0 kiB.
   Wait(s): 64(nvcsw) + 82900(nivcsw). Page Fault(s): 41301(minflt) + 0(majflt).
*/

/*
zz <- file("oSCa8ms64r25k17.dat", "rb")
l=readBin(zz,"integer",size=8,signed=F)
a=readBin(zz,"integer",l,size=1,signed=F)
dim(a)<-c(l/8192,8192)
b=colMeans(a)
png("oSCa8ms64r25k17.png",2048,600)
plot(x=1:length(b),y=b)
dev.off();
close(zz);
*/
