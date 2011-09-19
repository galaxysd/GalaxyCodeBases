// by Hu Xuesong
#ifndef _G_BIHASH_H
#define _G_BIHASH_H

#include "gtypendef.h"
#include <stddef.h> //size_t
#include <stdint.h> //uint64_t
#include <stdio.h>  //FILE
#include "gFileIO.h"
#include "cfgparser.h"

#define MAX_KMER_LEN (127)
#define HASH_LENB (128u)
#define SDL_SUBARRAY_UNIT (8u)
//#define SDLA_ITEMARRAY 32u
#define SUBARRAY_SIZE (32u*1024u)
//gcc -### -march=native -E /usr/include/stdlib.h 2>&1 | grep l1-cache-size
#define ENCODEDBIT_ENTROPY_PAD (1u);
/* 
http://math.stackexchange.com/questions/53353/is-there-some-reversible-mapping-that-as-uniform-as-a-hash 
Let $H$ be a hash function taking an arbitrary string as input and producing an $n$-bit string as output. 
Given an input string $s$, let $s_0$ denote the first $n$ bits of $s$ and $s_1$ the rest 
(i.e. $s = s_0\;||\;s_1$, where $s_0$ is $n$ bits long and $||$ denotes concatenation). 
Then define $$f(s) = (s_0 \oplus H(s_1))\;||\;s_1,$$ where $\oplus$ means bitwise XOR. 
It's easy to see that $f$ is its own inverse, i.e.$f(f(s)) = s$, 
so that the input string $s$ can be recovered from $f(s)$ just by running it through $f$ again. 

f(s)=(s0 XOR H(s1))||s1, length(s0)=EncodedBit, length(s1)=rBit 
 
EncodedBit < ArrayBit < rBit
*/ 
typedef struct __SDLeftArray_t {
    unsigned char CountBit, rBit, ArrayBit, EncodedBit;
    unsigned char itemByte; //, HashCnt;
    uint16_t SubItemCount,SubItemByUnit;  //max(SubItemCount) should be 32768
    size_t ArraySize,SDLAbyte;
    //uint64_t maxCount; == Item_CountBitMask
    //unsigned char ArrayCount;
    uint64_t ItemInsideAll, CellOverflowCount, CountBitOverflow; // ItemInsideAll = ItemInsideArray + CellOverflowCount
    double FalsePositiveRatio;
    void *pDLA, *pextree;
    uint64_t maxCountSeen;
    //uint64_t *outhash;
    uint64_t outhash[2];    // both ArrayBit and rBit is (0,64], so HashCnt==1 for MurmurHash3_x64_128
    uint128_t Item_rBitMask;
    uint64_t Hash_ArrayBitMask, Hash_rBitMask, Item_CountBitMask;
} SDLeftArray_t;

typedef struct __SDLdumpHead_t {
    char FileID[4]; //"GDSD"
    unsigned char FileVersion[2];    //0,1
    uint16_t kmersize;
    unsigned char CountBit, rBit;
    uint16_t SubItemCount;
    uint64_t ArraySize, SDLAbyte, extreebyte;
    uint64_t ItemInsideAll, CellOverflowCount, CountBitOverflow;
    uint64_t maxCountSeen;
    uint64_t HistMaxCntVal;
    uint64_t HistMaxHistVal;
    double HistMean;
    double HistSStd;
    uint32_t crc32c;
} __attribute__ ((packed)) SDLdumpHead;

SDLeftArray_t *dleft_arraynew(unsigned char CountBit, const SDLConfig * const psdlcfg, int kmer);
SDLeftArray_t *dleft_arrayinit(unsigned char CountBit, size_t ArraySize, uint16_t SubItemCount, int kmer);
size_t dleft_insert_read(unsigned int k, char const *const inseq, size_t len, SDLeftArray_t *dleftobj);

unsigned char GETitemByte_PADrBit_trimSubItemCount(const unsigned char CountBit, unsigned char *prBit, uint16_t *pSubItemCount);

void fprintSDLAnfo(FILE *stream, const SDLeftArray_t * dleftobj);
void dleft_arraydestroy(SDLeftArray_t * const dleftobj);
void dleft_dump(const SDLeftArray_t * const, SDLdumpHead * const, FILE *);

//#include "sdleftTF.h"
typedef struct __SDLeftStat_t {
    uint64_t HistMaxCntVal;
    uint64_t HistMaxHistVal;
    double HistMean;
    double HistSStd;
} SDLeftStat_t;

typedef SDLeftStat_t *(G_SDLeftArray_IN)(SDLeftArray_t * const, FILE *);    //*G_SDLeftArray_IN() is OK,too .

#ifdef TEST
SDLeftStat_t * dleft_stat(SDLeftArray_t * const dleftobj, FILE *stream, FILE *fpdat);
#else
SDLeftStat_t * dleft_stat(SDLeftArray_t * const dleftobj, FILE *stream);
#endif

#endif /* bihash.h */

/*  reprobe(i)=i(i+1)/2
    pos(m,i)=(hash(m)+reprobe(i)) mod M
*/
