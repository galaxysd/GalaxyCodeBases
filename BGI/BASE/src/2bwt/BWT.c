/*

   BWT.c    BWT-Index

   This module contains an implementation of BWT-index for alphabet size = 4.
   The functions provided include:
    Load functions for loading BWT to memory;
    Core functions for accessing core Inverse Psi values;
    Search functions for searching patterns from text;
    Text retrieval functions for retrieving text from BWT.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.L

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include "BWT.h"
#include "MiscUtilities.h"
#include "DNACount.h"
#include "TextConverter.h"
#include "MemManager.h"
#include "Socket.h"
#include "r250.h"
#include "HSP.h"
#include "HSPstatistic.h"

// static functions
static INLINE unsigned long long BWTOccValueExplicit(const BWT *bwt, const unsigned long long occIndexExplicit, const unsigned int character);
static INLINE void BWTAllOccValueExplicit(const BWT *bwt, const unsigned long long occIndexExplicit, unsigned long long* __restrict occValueExplicit);
static INLINE unsigned long long BWTSaIndexToChar(const BWT *bwt, const unsigned long long saIndex);
static INLINE unsigned long long BWTGetWordPackedText(const unsigned int *packedText, const unsigned long long index, const unsigned long long shift, const unsigned long long numOfBit);

static INLINE void BWTPrefetchOccValueExplicit(const BWT *bwt, const unsigned long long occIndexExplicit);
static INLINE void BWTPrefetchBWT(const BWT *bwt, const unsigned long long index);


int SaIndexGroupDPHitOrder1(const void *saIndexGroup, const long long index1, const long long index2);
int SaIndexGroupDPHitOrder2(const void *saIndexGroup, const long long index1, const long long index2);


static INLINE unsigned long long BWTSaIndexToChar(const BWT *bwt, const unsigned long long saIndex) {

    return (saIndex > bwt->cumulativeFreq[1]) + (saIndex > bwt->cumulativeFreq[2])
                                           + (saIndex > bwt->cumulativeFreq[3]);

}

BWT *BWTCreate(MMPool *mmPool, const unsigned long long textLength, unsigned long long  *decodeTable) {

    BWT *bwt;

    bwt = (BWT*)MMPoolDispatch(mmPool, sizeof(BWT));

    bwt->textLength = 0;
    bwt->inverseSa = 0;

    bwt->cumulativeFreq = (unsigned long long*) MMPoolDispatch(mmPool, (ALPHABET_SIZE + 1) * sizeof(unsigned long long));
    initializeVALLONG(bwt->cumulativeFreq, ALPHABET_SIZE + 1, 0);

    bwt->bwtSizeInWord = 0;
    bwt->saValueOnBoundary = NULL;

    // Generate decode tables
    if (decodeTable == NULL) {
        bwt->decodeTable = (unsigned long long*) MMPoolDispatch(mmPool, DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned long long));
        GenerateDNAOccCountTable(bwt->decodeTable);
    } else {
        bwt->decodeTable = decodeTable;
    }

    bwt->occMajorSizeInWord = BWTOccValueMajorSizeInWord(textLength);
    bwt->occValueMajor = (unsigned long long*) MMPoolDispatch(mmPool, bwt->occMajorSizeInWord * sizeof(unsigned long long));

    bwt->occSizeInWord = 0;
    bwt->occValue = NULL;

    bwt->saInterval = ALL_ONE_MASK_64;
    bwt->saValueSizeInWord = 0;
    bwt->saValue = NULL;

    bwt->inverseSaInterval = ALL_ONE_MASK_64;
    bwt->inverseSaSizeInWord = 0;
    bwt->inverseSa = NULL;

    return bwt;

}

BWT *BWTLoadLite(MMPool *mmPool, const char *bwtCodeFileName, const char *occValueFileName, unsigned long long *decodeTable) {
    unsigned long long i;
    FILE *bwtCodeFile, *occValueFile;
    BWT *bwt;
    unsigned long long bwtCodeLengthInFile;
    unsigned long long numOfCachedSaIndex;

    bwtCodeFile = (FILE*)fopen64(bwtCodeFileName, "rb");
    if (bwtCodeFile == NULL) {
        fprintf(stderr, "BWTLoad() : cannot open bwtCodeFile!\n");
        exit(1);
    }

    occValueFile = (FILE*)fopen64(occValueFileName, "rb");
    if (occValueFile == NULL) {
        fprintf(stderr, "BWTLoad() : cannot open occValueFile!\n");
        exit(1);
    }

    bwt = (BWT*) MMPoolDispatch(mmPool, sizeof(BWT));

    // fread(&bwt->inverseSa0, sizeof(unsigned long long), 1, bwtCodeFile);
    bwt->inverseSa0 = 0x3FFFFFFFFFFFFFFFLL; // hardcode as inf

    bwt->cumulativeFreq = (unsigned long long*) MMPoolDispatch(mmPool, (ALPHABET_SIZE + 1) * sizeof(unsigned long long));
    bwt->cumulativeFreq[0] = 0;
    fread(bwt->cumulativeFreq + 1, sizeof(unsigned long long), ALPHABET_SIZE, bwtCodeFile);
    bwt->textLength = bwt->cumulativeFreq[ALPHABET_SIZE];

    bwt->bwtSizeInWord = BWTResidentSizeInWord(bwt->textLength) + WORD_BETWEEN_OCC / 2;
    bwtCodeLengthInFile = BWTFileSizeInWord(bwt->textLength);
    bwt->bwtCode = (unsigned int*) MMUnitAllocate(bwt->bwtSizeInWord * sizeof(unsigned int));
    fread(bwt->bwtCode, sizeof(unsigned int), bwtCodeLengthInFile, bwtCodeFile);
    fclose(bwtCodeFile);
    BWTClearTrailingBwtCode(bwt);

    bwt->occSizeInWord = BWTOccValueMinorSizeInWord(bwt->textLength) ;
    bwt->occMajorSizeInWord = BWTOccValueMajorSizeInWord(bwt->textLength);
    bwt->occValue = (unsigned int*) MMUnitAllocate(bwt->occSizeInWord * sizeof(unsigned int));
    fread(bwt->occValue, sizeof(unsigned int), bwt->occSizeInWord, occValueFile);
    bwt->occValueMajor = (unsigned long long*) MMUnitAllocate(bwt->occMajorSizeInWord * sizeof(unsigned long long));
    fread(bwt->occValueMajor, sizeof(unsigned long long), bwt->occMajorSizeInWord, occValueFile);
    fclose(occValueFile);

    if (decodeTable == NULL) {
        bwt->decodeTable = (unsigned long long*) MMUnitAllocate(DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned long long));
        GenerateDNAOccCountTable(bwt->decodeTable);
        bwt->decodeTableGenerated = TRUE;
    } else {
        bwt->decodeTable = decodeTable;
        bwt->decodeTableGenerated = FALSE;
    }

    bwt->saValueOnBoundary = NULL;
    bwt->saInterval = ALL_ONE_MASK_64;
    bwt->saValueSizeInWord = 0;
    bwt->saValue = NULL;

    bwt->inverseSaInterval = ALL_ONE_MASK_64;
    bwt->inverseSaSizeInWord = 0;
    bwt->inverseSa = NULL;

    bwt->cachedSaIndex = (unsigned long long*) MMUnitAllocate((ALPHABET_SIZE + 1) * sizeof(unsigned long long));
    for (i=0;i<ALPHABET_SIZE;i++) {
        bwt->cachedSaIndex[i] = bwt->cumulativeFreq[i] + 1;
    }
    bwt->cachedSaIndex[ALPHABET_SIZE] = bwt->textLength + 1;    // To handle boundary case
    bwt->cachedSaIndexNumOfChar = 1;
    bwt->cachedSaIndexSizeInWord = (ALPHABET_SIZE + 1);
}

BWT *BWTLoad(MMPool *mmPool, const char *bwtCodeFileName, const char *occValueFileName, 
             const char *saValueFileName, const char *inverseSaFileName, const char *cachedSaIndexFileName,
             unsigned long long *decodeTable) {

    unsigned long long i;
    FILE *bwtCodeFile, *occValueFile, *saValueFile = NULL, *inverseSaFile = NULL, *cachedSaIndexFile = NULL;
    BWT *bwt;
    unsigned long long tmp;
    unsigned long long bwtCodeLengthInFile;
    unsigned long long numOfCachedSaIndex;

    bwtCodeFile = (FILE*)fopen64(bwtCodeFileName, "rb");
    if (bwtCodeFile == NULL) {
        fprintf(stderr, "BWTLoad() : cannot open bwtCodeFile:%s!\n", bwtCodeFileName);
        exit(1);
    }

    occValueFile = (FILE*)fopen64(occValueFileName, "rb");
    if (occValueFile == NULL) {
        fprintf(stderr, "BWTLoad() : cannot open occValueFile:%s!\n", occValueFileName);
        exit(1);
    }

    if (saValueFileName != NULL && saValueFileName[0] != '\0' && saValueFileName[0] != '-') {
        saValueFile = (FILE*)fopen64(saValueFileName, "rb");
        if (saValueFile == NULL) {
            fprintf(stderr, "BWTLoad() : cannot open saValueFile!\n");
            exit(1);
        }
    }

    if (inverseSaFileName != NULL && inverseSaFileName[0] != '\0' && inverseSaFileName[0] != '-') {
        inverseSaFile = (FILE*)fopen64(inverseSaFileName, "rb");
        if (inverseSaFile == NULL) {
            fprintf(stderr, "BWTLoad() : cannot open inverseSaFile!\n");
            exit(1);
        }
    }

    if (cachedSaIndexFileName != NULL && cachedSaIndexFileName[0] != '\0' && cachedSaIndexFileName[0] != '-') {
        cachedSaIndexFile = (FILE*)fopen64(cachedSaIndexFileName, "rb");
        if (cachedSaIndexFile == NULL) {
            fprintf(stderr, "BWTLoad() : cannot open cachedSaIndexFile!\n");
            exit(1);
        }
    }

    bwt = (BWT*) MMPoolDispatch(mmPool, sizeof(BWT));

    fread(&bwt->inverseSa0, sizeof(unsigned long long), 1, bwtCodeFile);

    bwt->cumulativeFreq = (unsigned long long*) MMPoolDispatch(mmPool, (ALPHABET_SIZE + 1) * sizeof(unsigned long long));
    bwt->cumulativeFreq[0] = 0;
    fread(bwt->cumulativeFreq + 1, sizeof(unsigned long long), ALPHABET_SIZE, bwtCodeFile);
    bwt->textLength = bwt->cumulativeFreq[ALPHABET_SIZE];

    fread(&tmp, sizeof(unsigned long long), 1, occValueFile);
    if (tmp != bwt->inverseSa0) {
        fprintf(stderr, "BWTLoad(): OccValue inverseSa0 not match!\n");
        exit(1);
    }
    for (i=1; i<=ALPHABET_SIZE; i++) {
        fread(&tmp, sizeof(unsigned long long), 1, occValueFile);
        if (tmp != bwt->cumulativeFreq[i]) {
            fprintf(stderr, "BWTLoad(): OccValue cumulativeFreq not match!\n");
            exit(1);
        }
    }

    bwt->bwtSizeInWord = BWTResidentSizeInWord(bwt->textLength) + WORD_BETWEEN_OCC / 2;    // + 8 words so that the 128 bits before and after an explicit occ are in the same aligned 64 byte
    bwtCodeLengthInFile = BWTFileSizeInWord(bwt->textLength);
    bwt->bwtCode = (unsigned int*) MMUnitAllocate(bwt->bwtSizeInWord * sizeof(unsigned int));
    fread(bwt->bwtCode, sizeof(unsigned int), bwtCodeLengthInFile, bwtCodeFile);
    fclose(bwtCodeFile);
    BWTClearTrailingBwtCode(bwt);

    bwt->occSizeInWord = BWTOccValueMinorSizeInWord(bwt->textLength) ;
    bwt->occMajorSizeInWord = BWTOccValueMajorSizeInWord(bwt->textLength);
    bwt->occValue = (unsigned int*) MMUnitAllocate(bwt->occSizeInWord * sizeof(unsigned int));
    fread(bwt->occValue, sizeof(unsigned int), bwt->occSizeInWord, occValueFile);
    bwt->occValueMajor = (unsigned long long*) MMUnitAllocate(bwt->occMajorSizeInWord * sizeof(unsigned long long));
    fread(bwt->occValueMajor, sizeof(unsigned long long), bwt->occMajorSizeInWord, occValueFile);
    fclose(occValueFile);

    if (decodeTable == NULL) {
        bwt->decodeTable = (unsigned long long*) MMUnitAllocate(DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned long long));
        GenerateDNAOccCountTable(bwt->decodeTable);
        bwt->decodeTableGenerated = TRUE;
    } else {
        bwt->decodeTable = decodeTable;
        bwt->decodeTableGenerated = FALSE;
    }

    bwt->saValueOnBoundary = NULL;
    if (saValueFile == NULL) {
        bwt->saInterval = ALL_ONE_MASK_64;
        bwt->saValueSizeInWord = 0;
        bwt->saValue = NULL;
    } else {
        fread(&tmp, sizeof(unsigned long long), 1, saValueFile);
        if (tmp != bwt->inverseSa0) {
            fprintf(stderr, "BWTLoad(): SaValue inverseSa0 not match!\n");
            exit(1);
        }
        for (i=1; i<=ALPHABET_SIZE; i++) {
            fread(&tmp, sizeof(unsigned long long), 1, saValueFile);
            if (tmp != bwt->cumulativeFreq[i]) {
                fprintf(stderr, "BWTLoad(): SaValue cumulativeFreq not match!\n");
                exit(1);
            }
        }
        fread(&bwt->saInterval, sizeof(unsigned long long), 1, saValueFile);
        bwt->saValueSizeInWord = (bwt->textLength + bwt->saInterval) / bwt->saInterval;
        bwt->saValue = (unsigned long long*) MMUnitAllocate(bwt->saValueSizeInWord * sizeof(unsigned long long));
        fread(bwt->saValue, sizeof(unsigned long long), bwt->saValueSizeInWord, saValueFile);
        bwt->saValue[0] = (unsigned long long)-1;    // Special handling for bwt
        fclose(saValueFile);

        BWTGenerateSaValueOnBoundary(mmPool, bwt);
    }

    if (inverseSaFile == NULL) {
        bwt->inverseSaInterval = ALL_ONE_MASK_64;
        bwt->inverseSaSizeInWord = 0;
        bwt->inverseSa = NULL;
    } else {
        fread(&tmp, sizeof(unsigned long long), 1, inverseSaFile);
        if (tmp != bwt->inverseSa0) {
            fprintf(stderr, "BWTLoad(): InverseSaValue inverseSa0 not match!\n");
            exit(1);
        }
        for (i=1; i<=ALPHABET_SIZE; i++) {
            fread(&tmp, sizeof(unsigned long long), 1, inverseSaFile);
            if (tmp != bwt->cumulativeFreq[i]) {
                fprintf(stderr, "BWTLoad(): InverseSaValue cumulativeFreq not match!\n");
                exit(1);
            }
        }
        fread(&bwt->inverseSaInterval, sizeof(unsigned long long), 1, inverseSaFile);
        bwt->inverseSaSizeInWord = (bwt->textLength + bwt->inverseSaInterval) / bwt->inverseSaInterval;
        bwt->inverseSa = (unsigned long long*) MMUnitAllocate(bwt->inverseSaSizeInWord * sizeof(unsigned long long));
        fread(bwt->inverseSa, sizeof(unsigned long long), bwt->inverseSaSizeInWord, inverseSaFile);
        fclose(inverseSaFile);
    }

    // Load Sa index range
    if (cachedSaIndexFile == NULL) {
        // Create a range from cumulative freq
        bwt->cachedSaIndex = (unsigned long long*) MMUnitAllocate((ALPHABET_SIZE + 1) * sizeof(unsigned long long));
        for (i=0;i<ALPHABET_SIZE;i++) {
            bwt->cachedSaIndex[i] = bwt->cumulativeFreq[i] + 1;
        }
        bwt->cachedSaIndex[ALPHABET_SIZE] = bwt->textLength + 1;    // To handle boundary case
        bwt->cachedSaIndexNumOfChar = 1;
        bwt->cachedSaIndexSizeInWord = (ALPHABET_SIZE + 1);
    } else {
        fread(&tmp, sizeof(unsigned long long), 1, cachedSaIndexFile);
        if (tmp != bwt->inverseSa0) {
            fprintf(stderr, "BWTLoad(): SaIndex inverseSa0 not match!\n");
            exit(1);
        }
        for (i=1; i<=ALPHABET_SIZE; i++) {
            fread(&tmp, sizeof(unsigned long long), 1, cachedSaIndexFile);
            if (tmp != bwt->cumulativeFreq[i]) {
                fprintf(stderr, "BWTLoad(): SaIndex cumulativeFreq not match!\n");
                exit(1);
            }
        }
        fread(&bwt->cachedSaIndexNumOfChar, sizeof(unsigned long long), 1, cachedSaIndexFile);
        numOfCachedSaIndex = 1 << (bwt->cachedSaIndexNumOfChar * 2);    // 4^cachedSaIndexNumOfChar
        bwt->cachedSaIndex = (unsigned long long*) MMUnitAllocate((numOfCachedSaIndex + 1) * sizeof(unsigned long long));
        fread(bwt->cachedSaIndex, sizeof(unsigned long long), numOfCachedSaIndex, cachedSaIndexFile);
        bwt->cachedSaIndex[numOfCachedSaIndex] = bwt->textLength + 1;    // To handle boundary case
        bwt->cachedSaIndexSizeInWord = (numOfCachedSaIndex + 1);
        fclose(cachedSaIndexFile);
    }

    return bwt;

}


void BWTFree(MMPool *mmPool, BWT *bwt) {

    MMPoolReturn(mmPool, bwt->cumulativeFreq, ALPHABET_SIZE * sizeof(unsigned long long));
    MMUnitFree(bwt->bwtCode, bwt->bwtSizeInWord * sizeof(unsigned int));

    if (bwt->occValue != NULL) {
        MMUnitFree(bwt->occValue, bwt->occSizeInWord * sizeof(unsigned int));
    }
    if (bwt->occValueMajor != NULL) {
        MMUnitFree(bwt->occValueMajor, bwt->occMajorSizeInWord * sizeof(unsigned long long));
    }

    if (bwt->saValue != NULL) {
        MMUnitFree(bwt->saValue, bwt->saValueSizeInWord * sizeof(unsigned long long));
    }
    if (bwt->inverseSa != NULL) {
        MMUnitFree(bwt->inverseSa, bwt->inverseSaSizeInWord * sizeof(unsigned long long));
    }

    if (bwt->decodeTableGenerated == TRUE) {
        MMUnitFree(bwt->decodeTable, DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned long long));
    }

    if (bwt->cachedSaIndex != NULL) {
        MMUnitFree(bwt->cachedSaIndex, bwt->cachedSaIndexSizeInWord * sizeof(unsigned long long));
    }

    if (bwt->saValueOnBoundary != NULL) {
        MMPoolReturn(mmPool, bwt->saValueOnBoundary, sizeof(unsigned long long) * 2 * ALPHABET_SIZE);
    }

    MMPoolReturn(mmPool, bwt, sizeof(BWT));

}

void BWTPrintMemoryUsage(const BWT *bwt, FILE *output, const unsigned long long packedDNASize) {

    unsigned long long totalMemorySize;

    fprintf(output, "BWT code size    : %llu\n", bwt->bwtSizeInWord * sizeof(unsigned int));
    fprintf(output, "Occ value size   : %llu\n", (bwt->occSizeInWord + bwt->occMajorSizeInWord) * sizeof(unsigned long long));
    if (bwt->saValueSizeInWord > 0) {
        fprintf(output, "SA value size    : %llu\n", bwt->saValueSizeInWord * sizeof(unsigned long long));
    }
    if (bwt->inverseSaSizeInWord > 0) {
        fprintf(output, "Inverse SA size  : %llu\n", bwt->inverseSaSizeInWord * sizeof(unsigned long long));
    }
    if (bwt->cachedSaIndex > 0) {
        fprintf(output, "SA index rangee  : %llu\n", bwt->cachedSaIndexSizeInWord * sizeof(unsigned long long));
    }
    if (packedDNASize > 0) {
        fprintf(output, "Packed DNA size  : %llu\n", packedDNASize);
    }
    
    totalMemorySize = (bwt->bwtSizeInWord) * sizeof(unsigned int) + (bwt->occSizeInWord + bwt->occMajorSizeInWord + bwt->saValueSizeInWord + bwt->inverseSaSizeInWord + bwt->cachedSaIndexSizeInWord) * sizeof(unsigned long long)
                       + packedDNASize;
    fprintf(output, "Total memory     : %llu\n", totalMemorySize);
    fprintf(output, "Bit per char     : %.2f\n", 
            (float)totalMemorySize / ((float)bwt->textLength / BITS_IN_BYTE));

}

void BWTGenerateSaValueOnBoundary(MMPool *mmPool, BWT *bwt) {

    unsigned long long i;

    if (bwt->saValueOnBoundary == NULL) {
        bwt->saValueOnBoundary = (unsigned long long*) MMPoolDispatch(mmPool, sizeof(unsigned long long) * 2 * ALPHABET_SIZE);
    }

    for (i=0; i<ALPHABET_SIZE; i++) {
        bwt->saValueOnBoundary[i * 2 + 1] = BWTSaValue(bwt, bwt->cumulativeFreq[i + 1]);
        if (bwt->cumulativeFreq[i] < bwt->textLength) {
            bwt->saValueOnBoundary[i * 2] = BWTSaValue(bwt, bwt->cumulativeFreq[i] + 1);
        } else {
            bwt->saValueOnBoundary[i * 2] = bwt->saValueOnBoundary[i * 2 + 1];
        }
    }

}

// Ordering of index1 and index2 is not important; this module will handle the ordering
// index1 and index2 can be on the same aligned 128 bit region or can be on adjacant aligned 128 bit region
// If index1 and index2 are in the same aligned 128 bit region, one of them must be on the boundary
// These requirements are to reduce the no. of branches in the program flow

unsigned long long BWTDecode(const BWT *bwt, const unsigned long long index1, const unsigned long long index2, const unsigned int character) {
    unsigned long long numChar1, numChar2, numChar3, numChar4, minIndex, maxIndex, minIndex128, maxIndex128;
    unsigned long long r,i,j,k,count=0;

    const static unsigned int ALIGN_16 partitionOne1[4] = { 23, 15, 7, 0 };
    const static unsigned int ALIGN_16 partitionOne4[4] = { 0, 7, 15, 23 };
    
    const static unsigned int ALIGN_16 partitionZero1[4]  = { 31, 23, 15, 7 };
    const static unsigned int ALIGN_16 partitionZero4[4]  = { 7, 15, 23, 31 };

    // SSE registers
    __m128i r1e, r2e, r3e, r4e;
    __m128i mcl;
    __m128i m0, m1, m2, m3;
    __m128i r1a, r1b, r1c, r1d;
    __m128i r2a, r2b, r2c, r2d;
    __m128i r3a, r3b, r3c, r3d;
    __m128i r4a, r4b, r4c, r4d;

    // to mask the highest bit of a character
    __m128i mhb;

    // Sort index1 and index2
    r = (index1 - index2) & -(index1 < index2);
    minIndex = index2 + r;
    maxIndex = index1 - r;

    // Locate 128 bit boundary
    //#define lastAlignedBoundary(offset, alignBoundary) ((offset) & (-alignBoundary))
    minIndex128 = lastAlignedBoundary(minIndex, CHAR_PER_128);
    maxIndex128 = lastAlignedBoundary(maxIndex - (maxIndex - minIndex > CHAR_PER_128), CHAR_PER_128);

    // Determine no.of characters to count
    numChar1 = maxIndex128 - minIndex;
    numChar2 = CHAR_PER_128 * (numChar1 > CHAR_PER_128);
    numChar3 = CHAR_PER_128 * (numChar1 > CHAR_PER_128*2);
    numChar4 = maxIndex - maxIndex128;
    numChar1 = numChar1 - numChar2 - numChar3;
    
//    printf("%u(%u) - %u(%u) [Difference=%u]\n",minIndex128,minIndex,maxIndex128,maxIndex,maxIndex-minIndex);
//    printf("%u-%u-%u-%u\n",numChar1,numChar2,numChar3,numChar4);
    
    //Load all necessary words
    r1e = _mm_load_si128((__m128i *)(bwt->bwtCode + minIndex128 / CHAR_PER_WORD));
    r2e = _mm_load_si128((__m128i *)(bwt->bwtCode + minIndex128 / CHAR_PER_WORD +4));
    r3e = _mm_load_si128((__m128i *)(bwt->bwtCode + minIndex128 / CHAR_PER_WORD +8));
    r4e = _mm_load_si128((__m128i *)(bwt->bwtCode + maxIndex128 / CHAR_PER_WORD));

    // to mask the highest bit of a character
    mhb = _mm_set1_epi32(0x77777777);
    r1e = _mm_and_si128(r1e, mhb);
    r2e = _mm_and_si128(r2e, mhb);
    r3e = _mm_and_si128(r3e, mhb);
    r4e = _mm_and_si128(r4e, mhb);
    
    // Set character extraction masks 
    m0 = _mm_set1_epi32(0xFFFFFFFF + (character & 1));    // Character selection mask for 0th bits
    m1 = _mm_set1_epi32(0xFFFFFFFF + ((character >> 1)&1));    // Character selection mask for 1st bits
    m2 = _mm_set1_epi32(0xFFFFFFFF + ((character >> 2)&1));    // Character selection mask for 2th bits
    m3 = _mm_set1_epi32(0xFFFFFFFF + ((character >> 3)&1));    // Character selection mask for 3st bits
    mcl = _mm_set1_epi32(0x11111111);                    // Set bit-clearing mask to 0x1111 1111....(alternate 1-bit)
 
    // This version of counting 1-4 x 128 bits whenever necessary
    // Count this pattern = numChar1[numChar2][numChar3]numChar4
    
    r1a = _mm_set1_epi32(numChar1);                   // Load number of characters into register
    r2a = _mm_set1_epi32(numChar2);                   // Load number of characters into register
    r3a = _mm_set1_epi32(numChar3);                   // Load number of characters into register
    r4a = _mm_set1_epi32(numChar4);                   // Load number of characters into register
    
    r1b = _mm_load_si128((__m128i*)partitionOne1);    // Load partition into register
    r1c = _mm_load_si128((__m128i*)partitionZero1);    // Load partition into register
    r4b = _mm_load_si128((__m128i*)partitionOne4);    // Load partition into register
    r4c = _mm_load_si128((__m128i*)partitionZero4);    // Load partition into register
    r2b = _mm_set1_epi32(CHAR_PER_128);                    // Load partition into register
    
    r1b = _mm_cmpgt_epi32(r1a, r1b);                // Compare to generate 4x32 bit mask; the word with counting boundary is all ones
    r1c = _mm_cmpgt_epi32(r1a, r1c);                // Compare to generate 4x32 bit mask; the word with counting boundary is all zeros
    r4b = _mm_cmpgt_epi32(r4a, r4b);                // Compare to generate 4x32 bit mask; the word with counting boundary is all ones
    r4c = _mm_cmpgt_epi32(r4a, r4c);                // Compare to generate 4x32 bit mask; the word with counting boundary is all zeros
    r2c = _mm_cmpeq_epi32(r2a, r2b);                // Compare to generate 4x32 bit mask; the word with counting boundary is all ones
    r3c = _mm_cmpeq_epi32(r3a, r2b);                // Compare to generate 4x32 bit mask; the word with counting boundary is all zeros
    
    r1b = _mm_srli_epi32(r1b, (8 - numChar1 % 8) * 4);    // Shift bits so that all word comform to the requirement of counting the word with counting boundary 
    r4b = _mm_slli_epi32(r4b, (8 - numChar4 % 8) * 4);    // Shift bits so that all word comform to the requirement of counting the word with counting boundary
    
    r1c = _mm_or_si128(r1b, r1c);    // Combine two masks
    r4c = _mm_or_si128(r4b, r4c);    // Combine two masks
    
    r1c = _mm_and_si128(r1c, mcl);    // Combine with bit-clearing mask (now = 0x11111111....)
    r2c = _mm_and_si128(r2c, mcl);    // Combine with bit-clearing mask (now = 0x11111111....)
    r3c = _mm_and_si128(r3c, mcl);    // Combine with bit-clearing mask (now = 0x11111111....)
    r4c = _mm_and_si128(r4c, mcl);    // Combine with bit-clearing mask (now = 0x11111111....)
    //0 d
    r1d = _mm_xor_si128(r1e, m0);    // Check all bits with mask
    r2d = _mm_xor_si128(r2e, m0);    // Check all bits with mask
    r3d = _mm_xor_si128(r3e, m0);    // Check all bits with mask
    r4d = _mm_xor_si128(r4e, m0);    // Check all bits with mask
    //1 b
    r1b = _mm_srli_epi32(r1e, 1);    // Shift encoding to right by 1 bit
    r2b = _mm_srli_epi32(r2e, 1);    // Shift encoding to right by 1 bit
    r3b = _mm_srli_epi32(r3e, 1);    // Shift encoding to right by 1 bit
    r4b = _mm_srli_epi32(r4e, 1);    // Shift encoding to right by 1 bit
    
    r1a = _mm_xor_si128(r1b, m1);    // Check all bits with mask
    r2a = _mm_xor_si128(r2b, m1);    // Check all bits with mask
    r3a = _mm_xor_si128(r3b, m1);    // Check all bits with mask
    r4a = _mm_xor_si128(r4b, m1);    // Check all bits with mask
    
    r1b = _mm_and_si128(r1a, r1d);    // Combine all bits
    r2b = _mm_and_si128(r2a, r2d);    // Combine all bits
    r3b = _mm_and_si128(r3a, r3d);    // Combine all bits
    r4b = _mm_and_si128(r4a, r4d);    // Combine all bits
    //2 a
    r1a = _mm_srli_epi32(r1e, 2);    // Shift encoding to right by 2 bit
    r2a = _mm_srli_epi32(r2e, 2);    // Shift encoding to right by 2 bit
    r3a = _mm_srli_epi32(r3e, 2);    // Shift encoding to right by 2 bit
    r4a = _mm_srli_epi32(r4e, 2);    // Shift encoding to right by 2 bit
    
    r1d = _mm_xor_si128(r1a, m2);    // Check all bits with mask
    r2d = _mm_xor_si128(r2a, m2);    // Check all bits with mask
    r3d = _mm_xor_si128(r3a, m2);    // Check all bits with mask
    r4d = _mm_xor_si128(r4a, m2);    // Check all bits with mask
    
    r1a = _mm_and_si128(r1b, r1d);    // Combine all bits
    r2a = _mm_and_si128(r2b, r2d);    // Combine all bits
    r3a = _mm_and_si128(r3b, r3d);    // Combine all bits
    r4a = _mm_and_si128(r4b, r4d);    // Combine all bits
    //3 d
    r1d = _mm_srli_epi32(r1e, 3);    // Shift encoding to right by 3 bit
    r2d = _mm_srli_epi32(r2e, 3);    // Shift encoding to right by 3 bit
    r3d = _mm_srli_epi32(r3e, 3);    // Shift encoding to right by 3 bit
    r4d = _mm_srli_epi32(r4e, 3);    // Shift encoding to right by 3 bit
    
    r1b = _mm_xor_si128(r1d, m3);    // Check all bits with mask
    r2b = _mm_xor_si128(r2d, m3);    // Check all bits with mask
    r3b = _mm_xor_si128(r3d, m3);    // Check all bits with mask
    r4b = _mm_xor_si128(r4d, m3);    // Check all bits with mask
    
    r1d = _mm_and_si128(r1a, r1b);    // Combine all bits
    r2d = _mm_and_si128(r2a, r2b);    // Combine all bits
    r3d = _mm_and_si128(r3a, r3b);    // Combine all bits
    r4d = _mm_and_si128(r4a, r4b);    // Combine all bits
    
    r1a = _mm_and_si128(r1d, r1c);    // Combine with counting mask, which has been combined with bit-clearing mask of 0x11111111.... 
    r2a = _mm_and_si128(r2d, r2c);    // Combine with counting mask, which has been combined with bit-clearing mask of 0x11111111.... 
    r3a = _mm_and_si128(r3d, r3c);    // Combine with counting mask, which has been combined with bit-clearing mask of 0x11111111.... 
    r4a = _mm_and_si128(r4d, r4c);    // Combine with counting mask, which has been combined with bit-clearing mask of 0x11111111.... 


    r1a = _mm_add_epi32(r1a, r2a);        // Combine 2 x 128 bits by adding them together
    r1a = _mm_add_epi32(r1a, r3a);        // Combine 2 x 128 bits by adding them together
    r1a = _mm_add_epi32(r1a, r4a);        // Combine 2 x 128 bits by adding them together
    
    // Combine 2 x 128 bits and continue counting
    mcl = _mm_set1_epi32(0x0F0F0F0F);    // Set bit-clearing mask to 0x0F0F0F0F....(alternate 4-bits)

    r1b = _mm_srli_epi32(r1a, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(r1a, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together

    m0 = _mm_setzero_si128();            // Set an all-zero mask
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit

    //Cross Debugging
            //unsigned int ALIGN_16 test[ALPHABET_SIZE];
            //BWTDecodeAll(bwt,index1,index2,test);
            //if (_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4)!=test[character])
            //    printf("#>%u %u\n",  _mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4),test[character]);

    return _mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
}

// Ordering of index1 and index2 is not important; this module will handle the ordering
// index1 and index2 can be on the same aligned 128 bit region or can be on adjacant aligned 128 bit region
// If index1 and index2 are in the same aligned 128 bit region, one of them must be on the boundary
// These requirements are to reduce the no. of branches in the program flow

void BWTDecodeAll(const BWT *bwt, const unsigned long long index1, const unsigned long long index2, unsigned long long* __restrict occValue) {
    unsigned long long numChar1, numChar2, numChar3, numChar4, minIndex, maxIndex, minIndex128, maxIndex128;
    unsigned long long r,i,j,k,count=0;

    const static unsigned int ALIGN_16 partitionOne1[4] = { 23, 15, 7, 0 };
    const static unsigned int ALIGN_16 partitionOne4[4] = { 0, 7, 15, 23 };
    
    const static unsigned int ALIGN_16 partitionZero1[4]  = { 31, 23, 15, 7 };
    const static unsigned int ALIGN_16 partitionZero4[4]  = { 7, 15, 23, 31 };

    // SSE registers
    __m128i r1e, r2e, r3e, r4e;
    __m128i mhb; // to mask the highest bit of a char (which stores the quality of a base)
    __m128i mcl;
    __m128i m0, m1;
    __m128i s0, s1, s2, s3;
    __m128i r0a, r0b, r0c, r0e;
    __m128i r1a, r1b, r1c, r1d;
    __m128i r2a, r2b, r2c, r2d;
    __m128i r3a, r3b, r3c, r3d;
    __m128i r4a, r4b, r4c, r4d;

    __m128i counting0,counting1,counting2,counting3,   
            counting4,counting5,counting6,counting7,
            counting8,counting9,counting10,counting11,
            counting12,counting13,counting14,counting15;

    counting0 = _mm_setzero_si128();            // Set an all-zero mask
    counting1 = _mm_setzero_si128();            // Set an all-zero mask
    counting2 = _mm_setzero_si128();            // Set an all-zero mask
    counting3 = _mm_setzero_si128();            // Set an all-zero mask
    counting4 = _mm_setzero_si128();            // Set an all-zero mask
    counting5 = _mm_setzero_si128();            // Set an all-zero mask
    counting6 = _mm_setzero_si128();            // Set an all-zero mask
    counting7 = _mm_setzero_si128();            // Set an all-zero mask
    counting8 = _mm_setzero_si128();            // Set an all-zero mask
    counting9 = _mm_setzero_si128();            // Set an all-zero mask
    counting10 = _mm_setzero_si128();            // Set an all-zero mask
    counting11 = _mm_setzero_si128();            // Set an all-zero mask
    counting12 = _mm_setzero_si128();            // Set an all-zero mask
    counting13 = _mm_setzero_si128();            // Set an all-zero mask
    counting14 = _mm_setzero_si128();            // Set an all-zero mask
    counting15 = _mm_setzero_si128();            // Set an all-zero mask


    // Sort index1 and index2
    r = (index1 - index2) & -(index1 < index2);
    minIndex = index2 + r;
    maxIndex = index1 - r;

    // Locate 128 bit boundary
    minIndex128 = lastAlignedBoundary(minIndex, CHAR_PER_128);
    maxIndex128 = lastAlignedBoundary(maxIndex - (maxIndex - minIndex > CHAR_PER_128), CHAR_PER_128);

    // Determine no.of characters to count
    numChar1 = maxIndex128 - minIndex;
    numChar2 = CHAR_PER_128 * (numChar1 > CHAR_PER_128);
    numChar3 = CHAR_PER_128 * (numChar1 > CHAR_PER_128*2);
    numChar4 = maxIndex - maxIndex128;
    numChar1 = numChar1 - numChar2 - numChar3;
    
//    printf("%u(%u) - %u(%u) [Difference=%u]\n",minIndex128,minIndex,maxIndex128,maxIndex,maxIndex-minIndex);
//    printf("%u-%u-%u-%u\n",numChar1,numChar2,numChar3,numChar4);
    
    //Load all necessary words
    r1e = _mm_load_si128((__m128i *)(bwt->bwtCode + minIndex128 / CHAR_PER_WORD));
    r2e = _mm_load_si128((__m128i *)(bwt->bwtCode + minIndex128 / CHAR_PER_WORD +4));
    r3e = _mm_load_si128((__m128i *)(bwt->bwtCode + minIndex128 / CHAR_PER_WORD +8));
    r4e = _mm_load_si128((__m128i *)(bwt->bwtCode + maxIndex128 / CHAR_PER_WORD));

    // to mask the highest bit of a character
    mhb = _mm_set1_epi32(0x77777777);
    r1e = _mm_and_si128(r1e, mhb);
    r2e = _mm_and_si128(r2e, mhb);
    r3e = _mm_and_si128(r3e, mhb);
    r4e = _mm_and_si128(r4e, mhb);
    
    // Set character extraction masks 
    m0 = _mm_setzero_si128();
    m1 = _mm_set1_epi32(0xFFFFFFFF);
    mcl = _mm_set1_epi32(0x11111111);                    // Set bit-clearing mask to 0x1111 1111....(alternate 1-bit)
 
    // This version of counting 1-4 x 128 bits whenever necessary
    // Count this pattern = numChar1[numChar2][numChar3]numChar4
    
    r1a = _mm_set1_epi32(numChar1);                   // Load number of characters into register
    r2a = _mm_set1_epi32(numChar2);                   // Load number of characters into register
    r3a = _mm_set1_epi32(numChar3);                   // Load number of characters into register
    r4a = _mm_set1_epi32(numChar4);                   // Load number of characters into register
    
    r1b = _mm_load_si128((__m128i*)partitionOne1);    // Load partition into register
    r1c = _mm_load_si128((__m128i*)partitionZero1);    // Load partition into register
    r4b = _mm_load_si128((__m128i*)partitionOne4);    // Load partition into register
    r4c = _mm_load_si128((__m128i*)partitionZero4);    // Load partition into register
    r2b = _mm_set1_epi32(CHAR_PER_128);                    // Load partition into register
    
    r1b = _mm_cmpgt_epi32(r1a, r1b);                // Compare to generate 4x32 bit mask; the word with counting boundary is all ones
    r1c = _mm_cmpgt_epi32(r1a, r1c);                // Compare to generate 4x32 bit mask; the word with counting boundary is all zeros
    r4b = _mm_cmpgt_epi32(r4a, r4b);                // Compare to generate 4x32 bit mask; the word with counting boundary is all ones
    r4c = _mm_cmpgt_epi32(r4a, r4c);                // Compare to generate 4x32 bit mask; the word with counting boundary is all zeros
    r2c = _mm_cmpeq_epi32(r2a, r2b);                // Compare to generate 4x32 bit mask; the word with counting boundary is all ones
    r3c = _mm_cmpeq_epi32(r3a, r2b);                // Compare to generate 4x32 bit mask; the word with counting boundary is all zeros
    
    r1b = _mm_srli_epi32(r1b, (8 - numChar1 % 8) * 4);    // Shift bits so that all word comform to the requirement of counting the word with counting boundary 
    r4b = _mm_slli_epi32(r4b, (8 - numChar4 % 8) * 4);    // Shift bits so that all word comform to the requirement of counting the word with counting boundary
    
    r1c = _mm_or_si128(r1b, r1c);    // Combine two masks
    r4c = _mm_or_si128(r4b, r4c);    // Combine two masks
    
    r1d = _mm_and_si128(r1c, mcl);    // Combine with bit-clearing mask (now = 0x11111111....)
    r2d = _mm_and_si128(r2c, mcl);    // Combine with bit-clearing mask (now = 0x11111111....)
    r3d = _mm_and_si128(r3c, mcl);    // Combine with bit-clearing mask (now = 0x11111111....)
    r4d = _mm_and_si128(r4c, mcl);    // Combine with bit-clearing mask (now = 0x11111111....)
    
        r1a = _mm_srli_epi32(r1e, 1);    // Shift encoding to right by 1 bit
        r1b = _mm_srli_epi32(r1e, 2);    // Shift encoding to right by 2 bit
        r1c = _mm_srli_epi32(r1e, 3);    // Shift encoding to right by 3 bit
        
        r0e = _mm_xor_si128(r1e, m1);    // Check if bit = ---0
        r0a = _mm_xor_si128(r1a, m1);    // Check if bit = --0-
        r0b = _mm_xor_si128(r1b, m1);    // Check if bit = -0--
        r0c = _mm_xor_si128(r1c, m1);    // Check if bit = 0---
        
        //Character 0
        s0 = _mm_and_si128(r0e, r0a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r1d);
        counting0 = _mm_add_epi32(counting0, s3);        
        //Character 8
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r1d);
        counting8 = _mm_add_epi32(counting8, s3);        
        //Character 4
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r1d);
        counting4 = _mm_add_epi32(counting4, s3);        
        //Character 12
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r1d);
        counting12 = _mm_add_epi32(counting12, s3);
        
        //Character 1
        s0 = _mm_and_si128(r1e, r0a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r1d);
        counting1 = _mm_add_epi32(counting1, s3);
        //Character 9
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r1d);
        counting9 = _mm_add_epi32(counting9, s3);
        //Character 5
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r1d);
        counting5 = _mm_add_epi32(counting5, s3);
        //Character 13
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r1d);
        counting13 = _mm_add_epi32(counting13, s3);
        
        //Character 2
        s0 = _mm_and_si128(r0e, r1a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r1d);
        counting2 = _mm_add_epi32(counting2, s3);
        //Character 10
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r1d);
        counting10 = _mm_add_epi32(counting10, s3);
        //Character 6
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r1d);
        counting6 = _mm_add_epi32(counting6, s3);
        //Character 14
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r1d);
        counting14 = _mm_add_epi32(counting14, s3);
        
        //Character 3
        s0 = _mm_and_si128(r1e, r1a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r1d);
        counting3 = _mm_add_epi32(counting3, s3);
        //Character 11
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r1d);
        counting11 = _mm_add_epi32(counting11, s3);
        //Character 7
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r1d);
        counting7 = _mm_add_epi32(counting7, s3);
        //Character 15
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r1d);
        counting15 = _mm_add_epi32(counting15, s3);
        

        r1a = _mm_srli_epi32(r2e, 1);    // Shift encoding to right by 1 bit
        r1b = _mm_srli_epi32(r2e, 2);    // Shift encoding to right by 2 bit
        r1c = _mm_srli_epi32(r2e, 3);    // Shift encoding to right by 3 bit
        
        r0e = _mm_xor_si128(r2e, m1);    // Check if bit = ---0
        r0a = _mm_xor_si128(r1a, m1);    // Check if bit = --0-
        r0b = _mm_xor_si128(r1b, m1);    // Check if bit = -0--
        r0c = _mm_xor_si128(r1c, m1);    // Check if bit = 0---
        
        //Character 0
        s0 = _mm_and_si128(r0e, r0a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r2d);
        counting0 = _mm_add_epi32(counting0, s3);        
        //Character 8
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r2d);
        counting8 = _mm_add_epi32(counting8, s3);        
        //Character 4
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r2d);
        counting4 = _mm_add_epi32(counting4, s3);        
        //Character 12
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r2d);
        counting12 = _mm_add_epi32(counting12, s3);
        
        //Character 1
        s0 = _mm_and_si128(r2e, r0a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r2d);
        counting1 = _mm_add_epi32(counting1, s3);
        //Character 9
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r2d);
        counting9 = _mm_add_epi32(counting9, s3);
        //Character 5
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r2d);
        counting5 = _mm_add_epi32(counting5, s3);
        //Character 13
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r2d);
        counting13 = _mm_add_epi32(counting13, s3);
        
        //Character 2
        s0 = _mm_and_si128(r0e, r1a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r2d);
        counting2 = _mm_add_epi32(counting2, s3);
        //Character 10
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r2d);
        counting10 = _mm_add_epi32(counting10, s3);
        //Character 6
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r2d);
        counting6 = _mm_add_epi32(counting6, s3);
        //Character 14
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r2d);
        counting14 = _mm_add_epi32(counting14, s3);
        
        //Character 3
        s0 = _mm_and_si128(r2e, r1a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r2d);
        counting3 = _mm_add_epi32(counting3, s3);
        //Character 11
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r2d);
        counting11 = _mm_add_epi32(counting11, s3);
        //Character 7
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r2d);
        counting7 = _mm_add_epi32(counting7, s3);
        //Character 15
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r2d);
        counting15 = _mm_add_epi32(counting15, s3);


        r1a = _mm_srli_epi32(r3e, 1);    // Shift encoding to right by 1 bit
        r1b = _mm_srli_epi32(r3e, 2);    // Shift encoding to right by 2 bit
        r1c = _mm_srli_epi32(r3e, 3);    // Shift encoding to right by 3 bit
        
        r0e = _mm_xor_si128(r3e, m1);    // Check if bit = ---0
        r0a = _mm_xor_si128(r1a, m1);    // Check if bit = --0-
        r0b = _mm_xor_si128(r1b, m1);    // Check if bit = -0--
        r0c = _mm_xor_si128(r1c, m1);    // Check if bit = 0---
        
        //Character 0
        s0 = _mm_and_si128(r0e, r0a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r3d);
        counting0 = _mm_add_epi32(counting0, s3);        
        //Character 8
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r3d);
        counting8 = _mm_add_epi32(counting8, s3);        
        //Character 4
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r3d);
        counting4 = _mm_add_epi32(counting4, s3);        
        //Character 12
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r3d);
        counting12 = _mm_add_epi32(counting12, s3);
        
        //Character 1
        s0 = _mm_and_si128(r3e, r0a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r3d);
        counting1 = _mm_add_epi32(counting1, s3);
        //Character 9
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r3d);
        counting9 = _mm_add_epi32(counting9, s3);
        //Character 5
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r3d);
        counting5 = _mm_add_epi32(counting5, s3);
        //Character 13
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r3d);
        counting13 = _mm_add_epi32(counting13, s3);
        
        //Character 2
        s0 = _mm_and_si128(r0e, r1a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r3d);
        counting2 = _mm_add_epi32(counting2, s3);
        //Character 10
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r3d);
        counting10 = _mm_add_epi32(counting10, s3);
        //Character 6
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r3d);
        counting6 = _mm_add_epi32(counting6, s3);
        //Character 14
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r3d);
        counting14 = _mm_add_epi32(counting14, s3);
        
        //Character 3
        s0 = _mm_and_si128(r3e, r1a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r3d);
        counting3 = _mm_add_epi32(counting3, s3);
        //Character 11
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r3d);
        counting11 = _mm_add_epi32(counting11, s3);
        //Character 7
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r3d);
        counting7 = _mm_add_epi32(counting7, s3);
        //Character 15
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r3d);
        counting15 = _mm_add_epi32(counting15, s3);
        

        r1a = _mm_srli_epi32(r4e, 1);    // Shift encoding to right by 1 bit
        r1b = _mm_srli_epi32(r4e, 2);    // Shift encoding to right by 2 bit
        r1c = _mm_srli_epi32(r4e, 3);    // Shift encoding to right by 3 bit
        
        r0e = _mm_xor_si128(r4e, m1);    // Check if bit = ---0
        r0a = _mm_xor_si128(r1a, m1);    // Check if bit = --0-
        r0b = _mm_xor_si128(r1b, m1);    // Check if bit = -0--
        r0c = _mm_xor_si128(r1c, m1);    // Check if bit = 0---
        
        //Character 0
        s0 = _mm_and_si128(r0e, r0a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r4d);
        counting0 = _mm_add_epi32(counting0, s3);        
        //Character 8
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r4d);
        counting8 = _mm_add_epi32(counting8, s3);        
        //Character 4
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r4d);
        counting4 = _mm_add_epi32(counting4, s3);        
        //Character 12
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r4d);
        counting12 = _mm_add_epi32(counting12, s3);
        
        //Character 1
        s0 = _mm_and_si128(r4e, r0a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r4d);
        counting1 = _mm_add_epi32(counting1, s3);
        //Character 9
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r4d);
        counting9 = _mm_add_epi32(counting9, s3);
        //Character 5
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r4d);
        counting5 = _mm_add_epi32(counting5, s3);
        //Character 13
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r4d);
        counting13 = _mm_add_epi32(counting13, s3);
        
        //Character 2
        s0 = _mm_and_si128(r0e, r1a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r4d);
        counting2 = _mm_add_epi32(counting2, s3);
        //Character 10
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r4d);
        counting10 = _mm_add_epi32(counting10, s3);
        //Character 6
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r4d);
        counting6 = _mm_add_epi32(counting6, s3);
        //Character 14
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r4d);
        counting14 = _mm_add_epi32(counting14, s3);
        
        //Character 3
        s0 = _mm_and_si128(r4e, r1a);
        s1 = _mm_and_si128(s0, r0b);
        s3 = _mm_and_si128(s1, r0c);
        s3 = _mm_and_si128(s3, r4d);
        counting3 = _mm_add_epi32(counting3, s3);
        //Character 11
        s3 = _mm_and_si128(s1, r1c);
        s3 = _mm_and_si128(s3, r4d);
        counting11 = _mm_add_epi32(counting11, s3);
        //Character 7
        s2 = _mm_and_si128(s0, r1b);
        s3 = _mm_and_si128(s2, r0c);
        s3 = _mm_and_si128(s3, r4d);
        counting7 = _mm_add_epi32(counting7, s3);
        //Character 15
        s3 = _mm_and_si128(s2, r1c);
        s3 = _mm_and_si128(s3, r4d);
        counting15 = _mm_add_epi32(counting15, s3);
        
    mcl = _mm_set1_epi32(0x0F0F0F0F);    // Set bit-clearing mask to 0xFFFFFFFF....(alternate 4-bits)

    r1b = _mm_srli_epi32(counting0, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting0, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[0]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
//    printf("0 %u %u",occValue[0],BWTDecode(bwt,index1,index2,0));
//    if (occValue[0]!=BWTDecode(bwt,index1,index2,0)) printf("\tX");
//    printf("\n");
    
    r1b = _mm_srli_epi32(counting1, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting1, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[1]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
//    printf("1 %u %u\n",occValue[1],BWTDecode(bwt,index1,index2,1));
//    if (occValue[1]!=BWTDecode(bwt,index1,index2,1)) printf("\tX");
//    printf("\n");
//    
    r1b = _mm_srli_epi32(counting2, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting2, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[2]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
//    printf("2 %u %u\n",occValue[2],BWTDecode(bwt,index1,index2,2));
//    if (occValue[2]!=BWTDecode(bwt,index1,index2,2)) printf("\tX");
//    printf("\n");
    
    r1b = _mm_srli_epi32(counting3, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting3, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[3]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
//    printf("3 %u %u\n",occValue[3],BWTDecode(bwt,index1,index2,3));
//    if (occValue[3]!=BWTDecode(bwt,index1,index2,3)) printf("\tX");
//    printf("\n");
    
    r1b = _mm_srli_epi32(counting4, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting4, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[4]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    
    r1b = _mm_srli_epi32(counting5, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting5, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[5]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    
    r1b = _mm_srli_epi32(counting6, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting6, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[6]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    
    r1b = _mm_srli_epi32(counting7, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting7, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[7]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register

    /*r1b = _mm_srli_epi32(counting8, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting8, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[8]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    
    r1b = _mm_srli_epi32(counting9, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting9, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[9]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    
    r1b = _mm_srli_epi32(counting10, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting10, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[10]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    
    r1b = _mm_srli_epi32(counting11, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting11, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[11]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    
    r1b = _mm_srli_epi32(counting12, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting12, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[12]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    
    r1b = _mm_srli_epi32(counting13, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting13, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[13]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    
    r1b = _mm_srli_epi32(counting14, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting14, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[14]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    
    r1b = _mm_srli_epi32(counting15, 4);        // Shift intermediate result to right by 4 bit
    r1a = _mm_and_si128(counting15, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    occValue[15]=_mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract and return result from register
    //*/
}


unsigned long long BWTOccValue(const BWT *bwt, unsigned long long index, const unsigned int character) {

    unsigned long long occValue, decodeValue;
    unsigned long long occExplicitIndex, occIndex;
    unsigned long long r;

    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is subtracted by 1 for adjustment
    index -= (index > bwt->inverseSa0);

#ifdef DEBUG
    if (index > bwt->textLength) {
        fprintf(stderr, "BWTOccValue() : index > textLength!\n");
        exit(1);
    }
#endif

    occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;    // Bidirectional encoding
    occIndex = occExplicitIndex * OCC_INTERVAL;


    occValue = BWTOccValueExplicit(bwt, occExplicitIndex, character);
#ifdef DEBUG
    if (occValue > occIndex) {
        fprintf(stderr, "BWTOccValue() : occValueExplicit > occIndex!\n");
        exit(1);
    }
#endif

    if (occIndex != index) {
        decodeValue = BWTDecode(bwt, occIndex, index, character);
        r = -(occIndex > index);
        return occValue + (decodeValue & ~r) - (decodeValue & r);
    } else {
        return occValue;
    }

}

void BWTOccValueTwoIndex(const BWT *bwt, unsigned long long index1, unsigned long long index2, const unsigned int character, unsigned long long* __restrict occValue) {

    unsigned long long decodeValue, tempExplicit1, tempExplicit2, tempOccValue1, tempOccValue2;
    unsigned long long occExplicitIndex1, occIndex1;
    unsigned long long occExplicitIndex2, occIndex2;
    unsigned long long r;

    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is subtracted by 1 for adjustment
    index1 -= (index1 > bwt->inverseSa0);
    index2 -= (index2 > bwt->inverseSa0);

#ifdef DEBUG
    if (index1 > bwt->textLength) {
        fprintf(stderr, "BWTOccValueTwoIndex() : index1 > textLength!\n");
        exit(1);
    }
    if (index2 > bwt->textLength) {
        fprintf(stderr, "BWTOccValueTwoIndex() : index2 > textLength!\n");
        exit(1);
    }
#endif

    // Pre-fetch memory to be accessed
    BWTPrefetchBWT(bwt, index1);
    BWTPrefetchBWT(bwt, index2);

    occExplicitIndex1 = (index1 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;    // Bidirectional encoding
    occIndex1 = occExplicitIndex1 * OCC_INTERVAL;
    occExplicitIndex2 = (index2 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;    // Bidirectional encoding
    occIndex2 = occExplicitIndex2 * OCC_INTERVAL;

    // Pre-fetch memory to be accessed
    BWTPrefetchOccValueExplicit(bwt, occExplicitIndex1);
    BWTPrefetchOccValueExplicit(bwt, occExplicitIndex2);


    if (occIndex1 != index1) {
        decodeValue = BWTDecode(bwt, occIndex1, index1, character);
        r = -(occIndex1 > index1);
        tempOccValue1 = (decodeValue & ~r) - (decodeValue & r);
    } else {
        tempOccValue1 = 0;
    }

    if (occIndex2 != index2) {
        decodeValue = BWTDecode(bwt, occIndex2, index2, character);
        r = -(occIndex2 > index2);
        tempOccValue2 = (decodeValue & ~r) - (decodeValue & r);
    } else {
        tempOccValue2 = 0;
    }

    tempExplicit1 = BWTOccValueExplicit(bwt, occExplicitIndex1, character);
    tempExplicit2 = BWTOccValueExplicit(bwt, occExplicitIndex2, character);
#ifdef DEBUG
    if (tempExplicit1 > occIndex1) {
        fprintf(stderr, "BWTOccValueTwoIndex() : occValueExplicit1 > occIndex1!\n");
        exit(1);
    }
    if (tempExplicit2 > occIndex2) {
        fprintf(stderr, "BWTOccValueTwoIndex() : occValueExplicit2 > occIndex2!\n");
        exit(1);
    }
#endif

    occValue[0] = tempOccValue1 + tempExplicit1;
    occValue[1] = tempOccValue2 + tempExplicit2;

}


void BWTAllOccValue(const BWT *bwt, unsigned long long index, unsigned long long* __restrict occValue) {

    unsigned long long occExplicitIndex, occIndex;
    unsigned long long ALIGN_16 tempOccValue[ALPHABET_SIZE];
    unsigned long long r;

    // SSE registers
    __m128i rtov0 ,rtov1, rtov2, rtov3 , rov0, rov1, rov2, rov3 , rc, t1, t2;

    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is subtracted by 1 for adjustment
    index -= (index > bwt->inverseSa0);

#ifdef DEBUG
    if (index > bwt->textLength) {
        fprintf(stderr, "BWTOccValue() : index > textLength!\n");
        exit(1);
    }
#endif

    occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;    // Bidirectional encoding
    occIndex = occExplicitIndex * OCC_INTERVAL;

    BWTAllOccValueExplicit(bwt, occExplicitIndex, occValue);

    if (occIndex != index) {

        BWTDecodeAll(bwt, occIndex, index, tempOccValue);

        // The following code add tempOccvalue to occValue if index > occIndex and subtract tempOccValue from occValue if occIndex > index
        r = -(occIndex > index);
        rc = _mm_set1_epi64((__m64)r);                // Set rc = r r r r
        
        rtov0 = _mm_load_si128((__m128i*)tempOccValue);
        rtov1 = _mm_load_si128((__m128i*)tempOccValue+1);
        rtov2 = _mm_load_si128((__m128i*)tempOccValue+2);
        rtov3 = _mm_load_si128((__m128i*)tempOccValue+3);
        
        rov0 = _mm_load_si128((__m128i*)occValue);
        rov1 = _mm_load_si128((__m128i*)occValue+1);
        rov2 = _mm_load_si128((__m128i*)occValue+2);
        rov3 = _mm_load_si128((__m128i*)occValue+3);
                
        t1 = _mm_andnot_si128(rc, rtov0);
        t2 = _mm_and_si128(rc, rtov0);
        rov0 = _mm_add_epi64(rov0, t1);
        rov0 = _mm_sub_epi64(rov0, t2);
        
        t1 = _mm_andnot_si128(rc, rtov1);
        t2 = _mm_and_si128(rc, rtov1);
        rov1 = _mm_add_epi64(rov1, t1);
        rov1 = _mm_sub_epi64(rov1, t2);
        
        t1 = _mm_andnot_si128(rc, rtov2);
        t2 = _mm_and_si128(rc, rtov2);
        rov2 = _mm_add_epi64(rov2, t1);
        rov2 = _mm_sub_epi64(rov2, t2);
        
        t1 = _mm_andnot_si128(rc, rtov3);
        t2 = _mm_and_si128(rc, rtov3);
        rov3 = _mm_add_epi64(rov3, t1);
        rov3 = _mm_sub_epi64(rov3, t2);
        
        _mm_store_si128((__m128i*)occValue, rov0);
        _mm_store_si128((__m128i*)occValue+1, rov1);
        _mm_store_si128((__m128i*)occValue+2, rov2);
        _mm_store_si128((__m128i*)occValue+3, rov3);

    } else {
        return;
    }

}

void BWTAllOccValueTwoIndex(const BWT *bwt, unsigned long long index1, unsigned long long index2, unsigned long long* __restrict occValue1, unsigned long long* __restrict occValue2) {

    unsigned long long occExplicitIndex1, occIndex1;
    unsigned long long occExplicitIndex2, occIndex2;
    unsigned long long ALIGN_16 tempOccValue1[ALPHABET_SIZE];
    unsigned long long ALIGN_16 tempOccValue2[ALPHABET_SIZE];
    unsigned long long r;

    // SSE registers
    __m128i rtov, rc, t1, t2, o1, o2;

    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is subtracted by 1 for adjustment
    index1 -= (index1 > bwt->inverseSa0);
    index2 -= (index2 > bwt->inverseSa0);

#ifdef DEBUG
    if (index1 > index2) {
        fprintf(stderr, "BWTAllOccValueTwoIndex() : index1 > index2!\n");
        exit(1);
    }
    if (index2 > bwt->textLength) {
        fprintf(stderr, "BWTAllOccValueTwoIndex() : index2 > textLength!\n");
        exit(1);
    }
#endif

    // Pre-fetch memory to be accessed
    BWTPrefetchBWT(bwt, index1);
    BWTPrefetchBWT(bwt, index2);

    occExplicitIndex1 = (index1 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;    // Bidirectional encoding
    occIndex1 = occExplicitIndex1 * OCC_INTERVAL;
    occExplicitIndex2 = (index2 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;    // Bidirectional encoding
    occIndex2 = occExplicitIndex2 * OCC_INTERVAL;

    // Pre-fetch memory to be accessed
    BWTPrefetchOccValueExplicit(bwt, occExplicitIndex1);
    BWTPrefetchOccValueExplicit(bwt, occExplicitIndex2);

    if (occIndex1 != index1) {

        BWTDecodeAll(bwt, occIndex1, index1, tempOccValue1);

        // The following code add tempOccvalue to occValue if index > occIndex and subtract tempOccValue from occValue if occIndex > index
        r = -(occIndex1 > index1);
        rtov = _mm_load_si128((__m128i*)tempOccValue1);
        rc = _mm_set1_epi32(r);                // Set rc = r r r r
        t1 = _mm_andnot_si128(rc, rtov);
        t2 = _mm_and_si128(rc, rtov);
        o1 = _mm_sub_epi32(t1, t2);
    } else {
        o1 = _mm_setzero_si128();
    }

    if (occIndex2 != index2) {

        BWTDecodeAll(bwt, occIndex2, index2, tempOccValue2);

        // The following code add tempOccvalue to occValue if index > occIndex and subtract tempOccValue from occValue if occIndex > index
        r = -(occIndex2 > index2);
        rc = _mm_set1_epi32(r);                // Set rc = r r r r
        rtov = _mm_load_si128((__m128i*)tempOccValue2);
        t1 = _mm_andnot_si128(rc, rtov);
        t2 = _mm_and_si128(rc, rtov);
        o2 = _mm_sub_epi32(t1, t2);

    } else {
        o2 = _mm_setzero_si128();
    }

    BWTAllOccValueExplicit(bwt, occExplicitIndex1, occValue1);
    BWTAllOccValueExplicit(bwt, occExplicitIndex2, occValue2);

    t1 = _mm_load_si128((__m128i*)occValue1);
    t2 = _mm_load_si128((__m128i*)occValue2);

    t1 = _mm_add_epi32(t1, o1);
    t2 = _mm_add_epi32(t2, o2);

    _mm_store_si128((__m128i*)occValue1, t1);
    _mm_store_si128((__m128i*)occValue2, t2);

}

unsigned long long BWTOccValueOnSpot(const BWT *bwt, unsigned long long index, unsigned long long* __restrict character) {

    unsigned long long occExplicitIndex, occIndex;
    unsigned long long occValue, decodeValue;
    unsigned long long r;

    // The bwt character before index will be returned and the count will be up to that bwt character
    #ifdef DEBUG
    if (index == bwt->inverseSa0 + 1) {
        fprintf(stderr, "BWTOccValueOnSpot(): index = inverseSa0 + 1!\n");
        exit(1);
    }
    if (index > bwt->textLength + 1) {
        fprintf(stderr, "BWTOccValueOnSpot() : index > textLength!\n");
        exit(1);
    }
    if (index == 0) {
        fprintf(stderr, "BWTOccValueOnSpot() : index = 0!\n");
        exit(1);
    }
    #endif

    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is incremented for adjustment
    index -= (index > bwt->inverseSa0);

    // Bidirectional encoding
    occExplicitIndex = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;
    occIndex = occExplicitIndex * OCC_INTERVAL;

    *character = (bwt->bwtCode[(index - 1) / CHAR_PER_WORD] << (((index - 1) % CHAR_PER_WORD) * BIT_PER_CHAR) >> (BITS_IN_WORD - BIT_PER_CHAR)) & 0x7777777777777777ULL;
    occValue = BWTOccValueExplicit(bwt, occExplicitIndex, *character);

    if (occIndex != index) {
        decodeValue = BWTDecode(bwt, occIndex, index, *character);
        r = -(occIndex > index);
        return occValue + (decodeValue & ~r) - (decodeValue & r);
    } else {
        return occValue;
    }

}

unsigned long long BWTSearchOccValue(const BWT *bwt, const unsigned int character, const unsigned long long searchOccValue) {

    unsigned long long occValue;
    unsigned long long i,j;
    unsigned int c , e;
    unsigned long long bwtPos;
    unsigned long long occExplicitIndexLeft, occExplicitIndexRight, occExplicitIndexMiddle;

    #ifdef DEBUG
    if (searchOccValue == 0 || searchOccValue > bwt->textLength) {
        fprintf(stderr, "BWTSearchOccValue() : searchOccValue out of bound!\n");
        exit(1);
    }
    #endif

    // Search Occurrence value

    occExplicitIndexLeft = 0;
    occExplicitIndexRight = (bwt->textLength + OCC_INTERVAL - 1) / OCC_INTERVAL;

    while (occExplicitIndexLeft + 1 < occExplicitIndexRight) {
        occExplicitIndexMiddle = average(occExplicitIndexLeft, occExplicitIndexRight);
        if (searchOccValue > BWTOccValueExplicit(bwt, occExplicitIndexMiddle, character)) {
            occExplicitIndexLeft = occExplicitIndexMiddle;
        } else {
            occExplicitIndexRight = occExplicitIndexMiddle;
        }
    }

    // Not tuned for DNA
    occValue = BWTOccValueExplicit(bwt, occExplicitIndexLeft, character);
    bwtPos = occExplicitIndexLeft * OCC_INTERVAL / CHAR_PER_WORD;

    for (i=0; i < OCC_INTERVAL / CHAR_PER_WORD; i++) {
        c = bwt->bwtCode[bwtPos + i] & 0x77777777;
        for (j=0; j < CHAR_PER_WORD && occValue < searchOccValue; j++) {
            if (c >> (BITS_IN_WORD - BIT_PER_CHAR) == character) {
                occValue++;
                if (occValue >= searchOccValue) {
                    return occExplicitIndexLeft * OCC_INTERVAL + i * CHAR_PER_WORD + j;
                }
            }
            c <<= BIT_PER_CHAR;
        }
    }

    fprintf(stderr, "BWTSearchOccValue() : unexpected error!\n");
    exit(1);

}

static INLINE unsigned long long BWTOccValueExplicit(const BWT *bwt, const unsigned long long occIndexExplicit, const unsigned int character) {

    unsigned long long occIndexMajor;
    unsigned long long compareMask, shift, mask;

    occIndexMajor = occIndexExplicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR;

    compareMask = (-(occIndexExplicit % OCC_VALUE_PER_WORD == 0));
    shift = 16 & compareMask;
    mask = 0x0000FFFF | compareMask;
    
    unsigned long long stackIdx = (occIndexExplicit / OCC_VALUE_PER_LONG) * ALPHABET_SIZE * OCC_VALUE_PER_LONG;
    stackIdx += character * OCC_VALUE_PER_LONG + occIndexExplicit % OCC_VALUE_PER_LONG;
    stackIdx /= OCC_VALUE_PER_WORD;
    
    return bwt->occValueMajor[occIndexMajor * ALPHABET_SIZE + character] +
            ((bwt->occValue[stackIdx] >> shift) & mask);

}

static INLINE void BWTAllOccValueExplicit(const BWT *bwt, const unsigned long long occIndexExplicit, unsigned long long* __restrict occValueExplicit) {

    // Ed Wu 26/05/2012
    // The BWT OccValue in original BWT-SW was a two-level sampling arrays.
    // The first level of sampling array stores the number of each 
    // character from 0 position to every 65536 (OCC_INTERVAL_MAJOR) position. The second
    // level of sampling array stores the number of each character from LAST 65536 position
    // to every 256 (OCC_INTERVAL) positions.
    //
    // <---------------------------------------BWT------------------------------------------>
    //      ^             ^             ^             ^             ^             ^            First Level
    //  > > | < < < > > > | < < < > > > | < < < > > > | < < < > > > | < < < > > > | < < <      Second Level
    //
    // Every first level count^ was a 32-bit word and every second level count (<, >) was a 16-bit integer.
    //
    // In BWT alphabet8_64 version this must be updated to
    // Every first level count^ was a 64-bit word and every second level count (<, >) was still a 16-bit integer.
    //

    unsigned long long occIndexMajor;
    unsigned long long compareMask,shift;
    __m64 mask;
    __m128i v1, v2, m;
    
    unsigned long long stackIdx = (occIndexExplicit / OCC_VALUE_PER_LONG) * ALPHABET_SIZE * OCC_VALUE_PER_LONG / OCC_VALUE_PER_WORD;

    occIndexMajor = occIndexExplicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR;

    compareMask = (-(occIndexExplicit % OCC_VALUE_PER_WORD == 0));
    shift = 32 * ((occIndexExplicit % OCC_VALUE_PER_LONG) / 2) + (16 & compareMask);
    mask =  (__m64) 0x000000000000FFFFLL;
    m = _mm_set1_epi64(mask);
    
    //=======
    v2 = _mm_load_si128((__m128i *)(bwt->occValue + stackIdx));
    v1 = _mm_load_si128((__m128i *)(bwt->occValueMajor + occIndexMajor * ALPHABET_SIZE));

    v2 = _mm_srli_epi64(v2, shift);
    v2 = _mm_and_si128(v2, m);

    v1 = _mm_add_epi64(v1, v2);

    _mm_store_si128((__m128i *)(occValueExplicit), v1);
    //=======
    v2 = _mm_load_si128((__m128i *)(bwt->occValue + stackIdx + 4));
    v1 = _mm_load_si128((__m128i *)(bwt->occValueMajor + occIndexMajor * ALPHABET_SIZE + 2));

    v2 = _mm_srli_epi64(v2, shift);
    v2 = _mm_and_si128(v2, m);

    v1 = _mm_add_epi64(v1, v2);

    _mm_store_si128((__m128i*)(occValueExplicit+2), v1);
    //=======
    v2 = _mm_load_si128((__m128i *)(bwt->occValue + stackIdx + 8));
    v1 = _mm_load_si128((__m128i *)(bwt->occValueMajor + occIndexMajor * ALPHABET_SIZE + 4));

    v2 = _mm_srli_epi64(v2, shift);
    v2 = _mm_and_si128(v2, m);

    v1 = _mm_add_epi64(v1, v2);

    _mm_store_si128((__m128i*)(occValueExplicit+4), v1);
    //=======
    v2 = _mm_load_si128((__m128i *)(bwt->occValue + stackIdx + 12));
    v1 = _mm_load_si128((__m128i *)(bwt->occValueMajor + occIndexMajor * ALPHABET_SIZE + 6));

    v2 = _mm_srli_epi64(v2, shift);
    v2 = _mm_and_si128(v2, m);

    v1 = _mm_add_epi64(v1, v2);

    _mm_store_si128((__m128i*)(occValueExplicit+6), v1);
}

static INLINE void BWTPrefetchOccValueExplicit(const BWT *bwt, const unsigned long long occIndexExplicit) {

    unsigned long long occIndexMajor;

    occIndexMajor = occIndexExplicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR;

    _mm_prefetch((char*)(bwt->occValue + occIndexExplicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE), _MM_HINT_T0);
    _mm_prefetch((char*)(bwt->occValueMajor + occIndexMajor * ALPHABET_SIZE), _MM_HINT_T0);

}

static INLINE void BWTPrefetchBWT(const BWT *bwt, const unsigned long long index) {

    _mm_prefetch((char*)(bwt->bwtCode + index / CHAR_PER_WORD), _MM_HINT_NTA);

}


unsigned long long BWTResidentSizeInWord(const unsigned long long numChar) {

    unsigned long long numCharRoundUpToOccInterval;

    // The $ in BWT at the position of inverseSa0 is not encoded
    numCharRoundUpToOccInterval = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL;

    return (numCharRoundUpToOccInterval + CHAR_PER_WORD - 1) / CHAR_PER_WORD;

}

unsigned long long BWTFileSizeInWord(const unsigned long long numChar) {

    // The $ in BWT at the position of inverseSa0 is not encoded
    return (numChar + CHAR_PER_WORD - 1) / CHAR_PER_WORD;

}

unsigned long long BWTOccValueMinorSizeInWord(const unsigned long long numChar) {

    unsigned long long numOfOccValue;

    numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;        // Value at both end for bi-directional encoding
    numOfOccValue = (numOfOccValue + OCC_VALUE_PER_LONG - 1) / OCC_VALUE_PER_LONG * OCC_VALUE_PER_LONG; // Align to OCC_VALUE_PER_LONG
    return (numOfOccValue + OCC_VALUE_PER_WORD - 1) / OCC_VALUE_PER_WORD * ALPHABET_SIZE;

}

unsigned long long BWTOccValueMajorSizeInWord(const unsigned long long numChar) {

    unsigned long long numOfOccValue;
    unsigned long long numOfOccIntervalPerMajor;

    numOfOccValue = (numChar + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;                // Value at both end for bi-directional encoding
    numOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;

    return (numOfOccValue + numOfOccIntervalPerMajor - 1) / numOfOccIntervalPerMajor * ALPHABET_SIZE;

}

void BWTClearTrailingBwtCode(BWT *bwt) {

    unsigned long long bwtResidentSizeInWord;
    unsigned long long wordIndex, offset;
    unsigned long long i;

    bwtResidentSizeInWord = BWTResidentSizeInWord(bwt->textLength);

    wordIndex = bwt->textLength / CHAR_PER_WORD;
    offset = (bwt->textLength - wordIndex * CHAR_PER_WORD) * BIT_PER_CHAR;
    if (offset > 0) {
        bwt->bwtCode[wordIndex] = truncateRight(bwt->bwtCode[wordIndex], BITS_IN_WORD - offset) & 0x77777777;
    } else {
        if (wordIndex < bwtResidentSizeInWord) {
            bwt->bwtCode[wordIndex] = 0;
        }
    }

    for (i=wordIndex+1; i<bwtResidentSizeInWord; i++) {
        bwt->bwtCode[i] = 0;
    }

}

unsigned long long BWTPsiMinusValue(const BWT *bwt, const unsigned long long index) {

    unsigned long long c;
    unsigned long long occValue;

    #ifdef DEBUG
    if (index > bwt->textLength) {
        fprintf(stderr, "BWTPsiMinusValue() : index out of range!\n");
        exit(1);
    }
    #endif

    if (index != bwt->inverseSa0) {

        occValue = BWTOccValueOnSpot(bwt, index + 1, &c);
        occValue += bwt->cumulativeFreq[c];

        return occValue;

    } else {
        return 0;
    }

}

unsigned long long BWTPsiPlusValue(const BWT *bwt, const unsigned long long index) {

    unsigned int c;
    unsigned long long psiPlusValue;

    #ifdef DEBUG
    if (index > bwt->textLength) {
        fprintf(stderr, "BWTPsiPlusValue() : index out of range!\n");
        exit(1);
    }
    #endif

    if (index == 0) {
        return bwt->inverseSa0;
    }

    // Find the BWT of PSI+
    c = (index > bwt->cumulativeFreq[1]) + (index > bwt->cumulativeFreq[2])
                                         + (index > bwt->cumulativeFreq[3]);

    psiPlusValue = BWTSearchOccValue(bwt, c, index - bwt->cumulativeFreq[c]);
    if (psiPlusValue >= bwt->inverseSa0) {
        psiPlusValue++;
    }
    return psiPlusValue;

}

unsigned long long BWTSaValue(const BWT *bwt, unsigned long long saIndex) {

    unsigned long long saValueSkipped = 0;

    #ifdef DEBUG
    if (saIndex > bwt->textLength) {
        fprintf(stderr, "BWTSaValue() : Index out of range!\n");
        exit(1);
    }
    if (bwt->saValue == NULL) {
        fprintf(stderr, "BWTSaValue() : Explicit SA value is not loaded!\n");
        exit(1);
    }
    #endif

    while (saIndex % bwt->saInterval != 0) {
        saValueSkipped++;
        saIndex = BWTPsiMinusValue(bwt, saIndex);
    }
    
    #ifdef DEBUG
    if (bwt->saValue[saIndex/bwt->saInterval] + saValueSkipped > bwt->textLength) {
        fprintf(stderr, "BWTSaValue() : saValue out of range!\n");
        exit(1);
    }
    #endif
    // SA[0] stores -1 although it should be textLength
    // PsiMinusValue returns 0 on inverseSa0
    return bwt->saValue[saIndex/bwt->saInterval] + saValueSkipped;

}

unsigned long long BWTInverseSa(const BWT *bwt, unsigned long long saValue) {

    unsigned long long i;
    unsigned long long saIndex;
    unsigned long long inverseSaExplicitIndex;
    unsigned long long saValueToSkip;

    #ifdef DEBUG
    if (saValue > bwt->textLength) {
        fprintf(stderr, "BWTInverseSa() : Index out of range!\n");
        exit(1);
    }
    if (bwt->inverseSa == NULL) {
        fprintf(stderr, "BWTInverseSa() : Explicit inverse SA is not loaded!\n");
        exit(1);
    }
    #endif

    inverseSaExplicitIndex = (saValue + bwt->inverseSaInterval - 1) / bwt->inverseSaInterval;
    if (inverseSaExplicitIndex * bwt->inverseSaInterval > bwt->textLength) {
        saIndex = 0;
        saValueToSkip = bwt->textLength - saValue;
    } else {
        saIndex = bwt->inverseSa[inverseSaExplicitIndex];
        saValueToSkip = inverseSaExplicitIndex * bwt->inverseSaInterval - saValue;
    }

    for (i=0; i<saValueToSkip; i++) {
        saIndex = BWTPsiMinusValue(bwt, saIndex);
    }

    return saIndex;

}

static INLINE unsigned long long BWTGetWordPackedText(const unsigned int *packedText, const unsigned long long index, const unsigned long long shift, const unsigned long long numOfBit) {

    unsigned long long text;
    const static unsigned int mask[32] = { 0x00000000, 0x80000000, 0xC0000000, 0xE0000000,
                                  0xF0000000, 0xF8000000, 0xFC000000, 0xFE000000,
                                  0xFF000000, 0xFF800000, 0xFFC00000, 0xFFE00000,
                                  0xFFF00000, 0xFFF80000, 0xFFFC0000, 0xFFFE0000,
                                  0xFFFF0000, 0xFFFF8000, 0xFFFFC000, 0xFFFFE000,
                                  0xFFFFF000, 0xFFFFF800, 0xFFFFFC00, 0xFFFFFE00,
                                  0xFFFFFF00, 0xFFFFFF80, 0xFFFFFFC0, 0xFFFFFFE0,
                                  0xFFFFFFF0, 0xFFFFFFF8, 0xFFFFFFFC, 0xFFFFFFFE };

    if (shift > 0) {
        // packedText should be allocated with at least 1 Word buffer initialized to zero
        text = (packedText[index] << shift) | (packedText[index + 1] >> (BITS_IN_WORD - shift));
    } else {
        text = packedText[index];
    }

    if (numOfBit < BITS_IN_WORD) {
        // Fill unused bit with zero
        text &= mask[numOfBit];
    }

    return text;
}

int BWTForwardSearch(const unsigned int *packedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText) {

    unsigned long long startSaIndex, endSaIndex, saIndexMiddle;
    unsigned long long saExplicitIndexLeft, saExplicitIndexRight, saExplicitIndexMiddle;
    unsigned long long saValue;

    unsigned long long firstChar;
    unsigned long long index, shift;
    unsigned long long packedKeyLength, keyLengthInBit;
    unsigned long long llcp, rlcp, mlcp, maxlcp;
    unsigned long long p = 0;    // to avoid compiler warning only

    if (keyLength % CHAR_PER_WORD == 0) {
        packedKeyLength = keyLength / CHAR_PER_WORD;
        keyLengthInBit = packedKeyLength * BITS_IN_WORD;
    } else {
        packedKeyLength = keyLength / CHAR_PER_WORD + 1;
        keyLengthInBit = (keyLength / CHAR_PER_WORD) * BITS_IN_WORD + 
                         (keyLength % CHAR_PER_WORD) * BIT_PER_CHAR;
    }

    // Get the SA index initial range by retrieving cumulative frequency
    firstChar = packedKey[0] >> (BITS_IN_WORD - BIT_PER_CHAR);

    startSaIndex = bwt->cumulativeFreq[firstChar] + 1;
    endSaIndex = bwt->cumulativeFreq[firstChar + 1];

    if (startSaIndex > endSaIndex) {
        // The first character of search pattern does not exists in text
        return 0;
    }

    // Find lcp for left boundary
    saValue = bwt->saValueOnBoundary[firstChar * 2];        // Pre-calculated

    // restriction for positions near the end of text
    maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);

    shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
    index = saValue / CHAR_PER_WORD;

    llcp = 0;
    while (llcp < maxlcp && packedKey[llcp] == 
                    BWTGetWordPackedText(packedText, index + llcp, shift, keyLengthInBit - llcp * BITS_IN_WORD)) {
        llcp++;
    }
    if ((saValue + keyLength > bwt->textLength) && llcp == maxlcp) {
        llcp--;
    }
    if (llcp == packedKeyLength) {
        return 1;
    }

    // Find lcp for right boundary
    saValue = bwt->saValueOnBoundary[firstChar * 2 + 1];    // Pre-calculated

    // restriction for positions near the end of text
    maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);

    shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
    index = saValue / CHAR_PER_WORD;

    rlcp = 0;
    while (rlcp < maxlcp && packedKey[rlcp] == 
                    BWTGetWordPackedText(packedText, index + rlcp, shift, keyLengthInBit - rlcp * BITS_IN_WORD)) {
        rlcp++;
    }
    if ((saValue + keyLength > bwt->textLength) && rlcp == maxlcp) {
        rlcp--;
    }
    if (rlcp == packedKeyLength) {
        return 1;
    }

    // Locate in SA index explicitly stored
    saExplicitIndexLeft = startSaIndex / bwt->saInterval;
    saExplicitIndexRight = (endSaIndex - 1) / bwt->saInterval + 1;

    // loop until two adjacent SA explicit index is found
    while (saExplicitIndexLeft + 1 < saExplicitIndexRight) {

        saExplicitIndexMiddle = average(saExplicitIndexLeft, saExplicitIndexRight);

        saValue = bwt->saValue[saExplicitIndexMiddle];
        shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
        index = saValue / CHAR_PER_WORD;

        // Try to increase mlcp
        mlcp = min(llcp, rlcp);        // mlcp = the characters (in unit of 16 for DNA) matched so far
        // restriction for positions near the end of text
        maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);

        while (mlcp < maxlcp) {
            p = BWTGetWordPackedText(packedText, index + mlcp, shift, keyLengthInBit - mlcp * BITS_IN_WORD);
            if (packedKey[mlcp] != p) {
                break;
            }
            mlcp++;
        }
        if ((saValue + keyLength <= bwt->textLength) || mlcp != maxlcp) {
            if (mlcp == packedKeyLength) {
                return 1;
            }
            if (packedKey[mlcp] > p) {
                llcp = mlcp;
                saExplicitIndexLeft = saExplicitIndexMiddle;
            } else {
                rlcp = mlcp;
                saExplicitIndexRight = saExplicitIndexMiddle;
            }
        } else {
            if (packedKey[mlcp-1] >= p) {
                llcp = mlcp - 1;
                saExplicitIndexLeft = saExplicitIndexMiddle;
            } else {
                rlcp = mlcp - 1;
                saExplicitIndexRight = saExplicitIndexMiddle;
            }
            
        }

    }

    // Two adjacent SA explicit index is found, convert back to SA index
    if (saExplicitIndexLeft == startSaIndex / bwt->saInterval) {
        startSaIndex = bwt->cumulativeFreq[firstChar] + 1;
    } else {
        startSaIndex = saExplicitIndexLeft * bwt->saInterval;
    }
    if (saExplicitIndexRight == (endSaIndex - 1) / bwt->saInterval + 1) {
        endSaIndex = bwt->cumulativeFreq[firstChar + 1];
    } else {
        endSaIndex = saExplicitIndexRight * bwt->saInterval;
    }

    // binary search by decoding bwt

    while (startSaIndex < endSaIndex) {

        saIndexMiddle = average(startSaIndex, endSaIndex);

        saValue = BWTSaValue(bwt, saIndexMiddle);
        shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
        index = saValue / CHAR_PER_WORD;

        // Try to increase mlcp
        mlcp = min(llcp, rlcp);        // mlcp = the characters (in unit of 16 for DNA) matched so far
        // restriction for positions near the end of text
        maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);

        while (mlcp < maxlcp) {
            p = BWTGetWordPackedText(packedText, index + mlcp, shift, keyLengthInBit - mlcp * BITS_IN_WORD);
            if (packedKey[mlcp] != p) {
                break;
            }
            mlcp++;
        }
        if ((saValue + keyLength <= bwt->textLength) || mlcp != maxlcp) {
            if (mlcp == packedKeyLength) {
                return 1;
            }
            if (packedKey[mlcp] > p) {
                llcp = mlcp;
                startSaIndex = saIndexMiddle + 1;
            } else {
                rlcp = mlcp;
                endSaIndex = saIndexMiddle;
            }
        } else {
            if (packedKey[mlcp-1] >= p) {
                llcp = mlcp - 1;
                startSaIndex = saIndexMiddle + 1;
            } else {
                rlcp = mlcp - 1;
                endSaIndex = saIndexMiddle;
            }
            
        }

    }

    // no match found
    return 0;

}

int BWTForwardSearchSaIndex(const unsigned int *packedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText, 
                                unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight) {

    unsigned long long startSaIndex, endSaIndex, saIndexMiddle;
    unsigned long long saExplicitIndexLeft, saExplicitIndexRight, saExplicitIndexMiddle;
    unsigned long long tempResultSaIndexLeft;
    unsigned long long tempSaExplicitIndexLeft, tempSaExplicitIndexRight, tempSaIndexLeft, tempSaIndexRight;
    unsigned long long saValue;

    unsigned long long firstChar;
    unsigned long long index, shift;
    unsigned long long packedKeyLength, keyLengthInBit;
    unsigned long long llcp, rlcp, mlcp, maxlcp;
    unsigned long long templlcp, temprlcp;
    unsigned long long p = 0;    // to avoid compiler warning only

    if (keyLength % CHAR_PER_WORD == 0) {
        packedKeyLength = keyLength / CHAR_PER_WORD;
        keyLengthInBit = packedKeyLength * BITS_IN_WORD;
    } else {
        packedKeyLength = keyLength / CHAR_PER_WORD + 1;
        keyLengthInBit = (keyLength / CHAR_PER_WORD) * BITS_IN_WORD + 
                         (keyLength % CHAR_PER_WORD) * BIT_PER_CHAR;
    }

    // Get the SA index initial range by retrieving cumulative frequency
    firstChar = packedKey[0] >> (BITS_IN_WORD - BIT_PER_CHAR);

    startSaIndex = bwt->cumulativeFreq[firstChar] + 1;
    endSaIndex = bwt->cumulativeFreq[firstChar + 1];

    if (startSaIndex > endSaIndex) {
        // The first character of search pattern does not exists in text
        return 0;
    }

    // Find lcp for left boundary
    saValue = bwt->saValueOnBoundary[firstChar * 2];        // Pre-calculated
    // restriction for positions near the end of text
    maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);

    llcp = 0;
    shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
    index = saValue / CHAR_PER_WORD;

    while (llcp < maxlcp && packedKey[llcp] == 
                    BWTGetWordPackedText(packedText, index + llcp, shift, keyLengthInBit - llcp * BITS_IN_WORD)) {
        llcp++;
    }
    if ((saValue + keyLength > bwt->textLength) && llcp == maxlcp) {
        llcp--;
    }

    // Find lcp for right boundary
    saValue = bwt->saValueOnBoundary[firstChar * 2 + 1];    // Pre-calculated
    shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
    index = saValue / CHAR_PER_WORD;

    // restriction for positions near the end of text
    maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);

    rlcp = 0;
    shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
    index = saValue / CHAR_PER_WORD;

    while (rlcp < maxlcp && packedKey[rlcp] == 
                    BWTGetWordPackedText(packedText, index + rlcp, shift, keyLengthInBit - rlcp * BITS_IN_WORD)) {
        rlcp++;
    }
    if ((saValue + keyLength > bwt->textLength) && rlcp == maxlcp) {
        rlcp--;
    }

    if (llcp == packedKeyLength && rlcp == packedKeyLength) {
        // Probably key is a character only
        *resultSaIndexLeft = startSaIndex;
        *resultSaIndexRight = endSaIndex;
        return 1;
    }

    // Locate in SA index explicitly stored
    saExplicitIndexLeft = startSaIndex / bwt->saInterval;
    saExplicitIndexRight = (endSaIndex - 1) / bwt->saInterval + 1;

    // Help determine where the search for ending boundary starts
    tempSaExplicitIndexLeft = saExplicitIndexLeft;
    tempSaExplicitIndexRight = saExplicitIndexRight;
    templlcp = llcp;
    temprlcp = rlcp;

    // loop until two adjacent SA explicit index is found
    while (saExplicitIndexLeft + 1 < saExplicitIndexRight) {

        saExplicitIndexMiddle = average(saExplicitIndexLeft, saExplicitIndexRight);

        saValue = bwt->saValue[saExplicitIndexMiddle];
        shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
        index = saValue / CHAR_PER_WORD;

        // Try to increase mlcp
        mlcp = min(llcp, rlcp);        // mlcp = the characters (in unit of 16 for DNA) matched so far
        // restriction for positions near the end of text
        maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);
        
        while (mlcp < maxlcp) {
            p = BWTGetWordPackedText(packedText, index + mlcp, shift, keyLengthInBit - mlcp * BITS_IN_WORD);
            if (packedKey[mlcp] != p) {
                break;
            }
            mlcp++;
        }
        if ((saValue + keyLength <= bwt->textLength) || mlcp != maxlcp) {
            if (mlcp < packedKeyLength && packedKey[mlcp] > p) {
                llcp = mlcp;
                saExplicitIndexLeft = saExplicitIndexMiddle;
                tempSaExplicitIndexLeft = max(tempSaExplicitIndexLeft, saExplicitIndexLeft);
                templlcp = llcp;
            } else {
                rlcp = mlcp;
                saExplicitIndexRight = saExplicitIndexMiddle;
                if (mlcp == packedKeyLength) {
                    tempSaExplicitIndexLeft = max(tempSaExplicitIndexLeft, saExplicitIndexRight);
                    templlcp = rlcp;
                } else {
                    tempSaExplicitIndexRight = saExplicitIndexRight;
                    temprlcp = rlcp;
                }
            }
        } else {
            if (packedKey[mlcp-1] >= p) {
                llcp = mlcp - 1;
                saExplicitIndexLeft = saExplicitIndexMiddle;
                tempSaExplicitIndexLeft = max(tempSaExplicitIndexLeft, saExplicitIndexLeft);
                templlcp = llcp;
            } else {
                rlcp = mlcp - 1;
                saExplicitIndexRight = saExplicitIndexMiddle;
                tempSaExplicitIndexRight = saExplicitIndexRight;
                temprlcp = rlcp;
            }
        }

    }

    // Help determine where the search for ending boundary starts
    if (tempSaExplicitIndexLeft == startSaIndex / bwt->saInterval) {
        tempSaIndexLeft = bwt->cumulativeFreq[firstChar] + 1;
    } else {
        tempSaIndexLeft = tempSaExplicitIndexLeft * bwt->saInterval;
    }
    if (tempSaExplicitIndexRight == (endSaIndex - 1) / bwt->saInterval + 1) {
        tempSaIndexRight = bwt->cumulativeFreq[firstChar + 1];
    } else {
        tempSaIndexRight = tempSaExplicitIndexRight * bwt->saInterval;
    }

    // Two adjacent SA explicit index is found, convert back to SA index
    if (saExplicitIndexLeft == startSaIndex / bwt->saInterval) {
        startSaIndex = bwt->cumulativeFreq[firstChar] + 1;
    } else {
        startSaIndex = saExplicitIndexLeft * bwt->saInterval;
    }
    if (saExplicitIndexRight == (endSaIndex - 1) / bwt->saInterval + 1) {
        endSaIndex = bwt->cumulativeFreq[firstChar + 1];
    } else {
        endSaIndex = saExplicitIndexRight * bwt->saInterval;
    }

    // binary search by decoding bwt

    while (startSaIndex < endSaIndex) {

        saIndexMiddle = average(startSaIndex, endSaIndex);

        saValue = BWTSaValue(bwt, saIndexMiddle);
        shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
        index = saValue / CHAR_PER_WORD;

        // Try to increase mlcp
        mlcp = min(llcp, rlcp);        // mlcp = the characters (in unit of 16 for DNA) matched so far
        // restriction for positions near the end of text
        maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);
        
        while (mlcp < maxlcp) {
            p = BWTGetWordPackedText(packedText, index + mlcp, shift, keyLengthInBit - mlcp * BITS_IN_WORD);
            if (packedKey[mlcp] != p) {
                break;
            }
            mlcp++;
        }
        if ((saValue + keyLength <= bwt->textLength) || mlcp != maxlcp) {
            if (mlcp < packedKeyLength && packedKey[mlcp] > p) {
                llcp = mlcp;
                startSaIndex = saIndexMiddle + 1;
                tempSaIndexLeft = max(tempSaIndexLeft, startSaIndex);
                templlcp = llcp;
            } else {
                rlcp = mlcp;
                endSaIndex = saIndexMiddle;
                if (mlcp == packedKeyLength) {
                    tempSaIndexLeft = max(tempSaIndexLeft, endSaIndex);
                    templlcp = rlcp;
                } else {
                    tempSaIndexRight = endSaIndex;
                    temprlcp = rlcp;
                }
            }
        } else {
            if (packedKey[mlcp-1] >= p) {
                llcp = mlcp - 1;
                startSaIndex = saIndexMiddle + 1;
                tempSaIndexLeft = max(tempSaIndexLeft, startSaIndex);
                templlcp = llcp;
            } else {
                rlcp = mlcp - 1;
                endSaIndex = saIndexMiddle;
                tempSaIndexRight = endSaIndex;
                temprlcp = rlcp;
            }
        }

    }

    if (max(llcp, rlcp) < packedKeyLength) {
        // no match found
        return 0;
    }

    // The starting SA index found
    tempResultSaIndexLeft = startSaIndex;

    // search for the ending SA index

    // binary search by decoding bwt

    startSaIndex = tempSaIndexLeft;
    endSaIndex = tempSaIndexRight;
    llcp = templlcp;
    rlcp = temprlcp;

    while (startSaIndex < endSaIndex) {

        saIndexMiddle = average(startSaIndex, endSaIndex + 1);

        saValue = BWTSaValue(bwt, saIndexMiddle);
        shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
        index = saValue / CHAR_PER_WORD;

        // Try to increase mlcp
        mlcp = min(llcp, rlcp);        // mlcp = the characters (in unit of 16 for DNA) matched so far
        // restriction for positions near the end of text
        maxlcp = min(packedKeyLength, (bwt->textLength - saValue + CHAR_PER_WORD - 1) / CHAR_PER_WORD);
        
        while (mlcp < maxlcp) {
            p = BWTGetWordPackedText(packedText, index + mlcp, shift, keyLengthInBit - mlcp * BITS_IN_WORD);
            if (packedKey[mlcp] != p) {
                break;
            }
            mlcp++;
        }
        if ((saValue + keyLength <= bwt->textLength) || mlcp != maxlcp) {
            if (mlcp == packedKeyLength || packedKey[mlcp] > p) {
                llcp = mlcp;
                startSaIndex = saIndexMiddle;
            } else {
                rlcp = mlcp;
                endSaIndex = saIndexMiddle - 1;
            }
        } else {
            if (packedKey[mlcp-1] >= p) {
                llcp = mlcp - 1;
                startSaIndex = saIndexMiddle;
            } else {
                rlcp = mlcp - 1;
                endSaIndex = saIndexMiddle - 1;
            }
        }

    }

    *resultSaIndexLeft = tempResultSaIndexLeft;
    *resultSaIndexRight = endSaIndex;

    return 1;

}

int BWTSaBinarySearch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText, 
                      unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight, unsigned int *tempKey) {    // tempKey = buffer large enough to hold packed key

    unsigned long long saExplicitIndexLeft, saExplicitIndexRight, saExplicitIndexMiddle;
    unsigned long long saIndexLeft, saIndexRight, saIndexMiddle;
    unsigned long long saValue;

    unsigned long long saRangeIndex;
    unsigned long long index, shift;
    unsigned long long llcp, rlcp, mlcp;

    unsigned long long pos;
    unsigned long long cachedNumOfChar;

    unsigned long long i, j;

    unsigned long long numOfWord, numOfFullWord, numOfOddChar;
    unsigned long long text;
    unsigned long long tempKeyATrailing, tempKeyTTrailing;

    unsigned long long initialSaIndexLeft, initialSaIndexRight;
    unsigned long long stage1SaExplicitIndexLeft, stage1SaExplicitIndexRight;
    unsigned long long stage1llcp, stage1rlcp;
    unsigned long long stage2StartSaExplicitIndexLeft, stage2Startllcp, stage2Startrlcp;
    unsigned long long stage2EndSaExplicitIndexLeft, stage2Endllcp, stage2Endrlcp;
    unsigned long long stage3SaIndexLeft, stage3SaIndexRight;
    unsigned long long stage3StartSaIndexLeft, stage3StartSaIndexRight, stage3EndSaIndexLeft, stage3EndSaIndexRight;

    // Get SA index range from cached SA index range

    cachedNumOfChar = min(bwt->cachedSaIndexNumOfChar, keyLength);

    saRangeIndex = 0;
    for (pos = 0; pos < cachedNumOfChar; pos++) {
        saRangeIndex <<= BIT_PER_CHAR;
        saRangeIndex |= convertedKey[pos];
    }

    initialSaIndexLeft = bwt->cachedSaIndex[saRangeIndex << ((bwt->cachedSaIndexNumOfChar - cachedNumOfChar) * BIT_PER_CHAR)];
    initialSaIndexRight = bwt->cachedSaIndex[(saRangeIndex + 1) << ((bwt->cachedSaIndexNumOfChar - cachedNumOfChar) * BIT_PER_CHAR)] - 1;

    if (initialSaIndexLeft > initialSaIndexRight || keyLength == cachedNumOfChar) {
        *resultSaIndexLeft = initialSaIndexLeft;
        *resultSaIndexRight = initialSaIndexRight;
        return (initialSaIndexLeft <= initialSaIndexRight);
    }

    // Pack key into temp
    numOfWord = (keyLength - cachedNumOfChar + CHAR_PER_WORD - 1) / CHAR_PER_WORD;
    numOfFullWord = (keyLength - cachedNumOfChar) / CHAR_PER_WORD;
    numOfOddChar = keyLength - cachedNumOfChar - numOfFullWord * CHAR_PER_WORD;

    for (i=0; i<numOfFullWord; i++) {
        tempKey[i] = 0;
        for (j=0; j<CHAR_PER_WORD; j++) {
            tempKey[i] <<= BIT_PER_CHAR;
            tempKey[i] |= convertedKey[pos];
            pos++;
        }
    }
    if (numOfWord > numOfFullWord) {
        tempKey[i] = 0;
        for (j=0; j<numOfOddChar; j++) {
            tempKey[i] <<= BIT_PER_CHAR;
            tempKey[i] |= convertedKey[pos];
            pos++;
        }
        tempKey[i] <<= BITS_IN_WORD - numOfOddChar * BIT_PER_CHAR;
    }

    tempKeyATrailing = tempKey[numOfWord - 1];
    if (numOfOddChar) {
        tempKeyTTrailing = tempKeyATrailing | (ALL_ONE_MASK >> (numOfOddChar * BIT_PER_CHAR));
    } else {
        tempKeyTTrailing = tempKeyATrailing;
    }

    // Stage 1: search for an SA index where all full words are matched
    saExplicitIndexLeft = initialSaIndexLeft / bwt->saInterval;
    saExplicitIndexRight = (initialSaIndexRight + bwt->saInterval - 1) / bwt->saInterval;

    llcp = 0;
    rlcp = 0;
    mlcp = 0;

    while (mlcp < numOfFullWord && (saExplicitIndexLeft + 1) < saExplicitIndexRight) {

        saExplicitIndexMiddle = average(saExplicitIndexLeft, saExplicitIndexRight);

        saValue = bwt->saValue[saExplicitIndexMiddle] + cachedNumOfChar;
        shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
        index = saValue / CHAR_PER_WORD;

        // Try to increase mlcp
        mlcp = min(llcp, rlcp);        // mlcp = the characters (in unit of 16 for DNA) matched so far
        
        do {
            if (shift != 0) {
                text = (packedText[index + mlcp] << shift) | (packedText[index + mlcp + 1] >> (BITS_IN_WORD - shift));
            } else {
                text = packedText[index + mlcp];
            }
        } while (tempKey[mlcp] == text && ++mlcp < numOfFullWord);

        if (mlcp < numOfFullWord) {
            if (tempKey[mlcp] > text) {
                saExplicitIndexLeft = saExplicitIndexMiddle;
                llcp = mlcp;
            } else {
                saExplicitIndexRight = saExplicitIndexMiddle;
                rlcp = mlcp;
            }
        }

    }

    stage1SaExplicitIndexLeft = saExplicitIndexLeft;
    stage1SaExplicitIndexRight = saExplicitIndexRight;
    stage1llcp = llcp;
    stage1rlcp = rlcp;

    // Store stage 1 result: stage1SaExplicitIndexLeft > key; stage1SaExplicitIndexRight < key
    //                         either (i)  stage1SaExplicitIndexLeft + 1 = stage1SaExplicitIndexRight or
    //                                (ii) all full words are matched somewhere between stage1SaExplicitIndexLeft and stage1SaExplicitIndexRight inclusive

    // Stage 2: locate the starting SA index and ending SA index separately

    // Search for starting SA index

    //saExplicitIndexLeft = stage1SaExplicitIndexLeft;
    //saExplicitIndexRight = stage1SaExplicitIndexRight;
    //llcp = stage1llcp;
    //rlcp = stage1rlcp;
    //tempKey[numOfWord - 1] = tempKeyATrailing;

    mlcp = 0;
    while ((saExplicitIndexLeft + 1) < saExplicitIndexRight) {

        saExplicitIndexMiddle = average(saExplicitIndexLeft, saExplicitIndexRight);

        saValue = bwt->saValue[saExplicitIndexMiddle] + cachedNumOfChar;
        shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
        index = saValue / CHAR_PER_WORD;

        // Try to increase mlcp
        mlcp = min(llcp, rlcp);        // mlcp = the characters (in unit of 16 for DNA) matched so far
        
        do  {
            if (shift != 0) {
                text = (packedText[index + mlcp] << shift) | (packedText[index + mlcp + 1] >> (BITS_IN_WORD - shift));
            } else {
                text = packedText[index + mlcp];
            }
        } while (tempKey[mlcp] == text && ++mlcp < numOfWord);

        if (mlcp < numOfWord && tempKey[mlcp] > text) {
            saExplicitIndexLeft = saExplicitIndexMiddle;
            llcp = mlcp;
        } else {
            saExplicitIndexRight = saExplicitIndexMiddle;
            rlcp = mlcp;
        }
    }

    stage2StartSaExplicitIndexLeft = saExplicitIndexLeft;
    stage2Startllcp = llcp;
    stage2Startrlcp = rlcp;

    // Search for ending SA index

    saExplicitIndexLeft = stage1SaExplicitIndexLeft;
    saExplicitIndexRight = stage1SaExplicitIndexRight;
    llcp = stage1llcp;
    rlcp = stage1rlcp;
    tempKey[numOfWord - 1] = tempKeyTTrailing;

    mlcp = 0;
    while ((saExplicitIndexLeft + 1) < saExplicitIndexRight) {

        saExplicitIndexMiddle = average(saExplicitIndexLeft, saExplicitIndexRight);

        saValue = bwt->saValue[saExplicitIndexMiddle] + cachedNumOfChar;
        shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
        index = saValue / CHAR_PER_WORD;

        // Try to increase mlcp
        mlcp = min(llcp, rlcp);        // mlcp = the characters (in unit of 16 for DNA) matched so far
        
        do {
            if (shift != 0) {
                text = (packedText[index + mlcp] << shift) | (packedText[index + mlcp + 1] >> (BITS_IN_WORD - shift));
            } else {
                text = packedText[index + mlcp];
            }
        } while (tempKey[mlcp] == text && ++mlcp < numOfWord);

        if (mlcp >= numOfWord || tempKey[mlcp] >= text) {
            saExplicitIndexLeft = saExplicitIndexMiddle;
            llcp = mlcp;
        } else {
            saExplicitIndexRight = saExplicitIndexMiddle;
            rlcp = mlcp;
        }

    }

    stage2EndSaExplicitIndexLeft = saExplicitIndexLeft;
    stage2Endllcp = llcp;
    stage2Endrlcp = rlcp;

    // Not found
    if (stage2StartSaExplicitIndexLeft > stage2EndSaExplicitIndexLeft) {
        return 0;
    }

    // Stage 3: search while decoding SA using BWT
    if (stage2StartSaExplicitIndexLeft * bwt->saInterval > initialSaIndexLeft) {
        stage3StartSaIndexLeft = stage2StartSaExplicitIndexLeft * bwt->saInterval;
    } else {
        stage3StartSaIndexLeft = initialSaIndexLeft - 1;
    }
    if ((stage2StartSaExplicitIndexLeft + 1) * bwt->saInterval < initialSaIndexRight) {
        stage3StartSaIndexRight = (stage2StartSaExplicitIndexLeft + 1) * bwt->saInterval;
    } else {
        stage3StartSaIndexRight = initialSaIndexRight + 1;
    }
    if (stage2EndSaExplicitIndexLeft * bwt->saInterval > initialSaIndexLeft) {
        stage3EndSaIndexLeft = stage2EndSaExplicitIndexLeft * bwt->saInterval;
    } else {
        stage3EndSaIndexLeft = initialSaIndexLeft - 1;
    }
    if ((stage2EndSaExplicitIndexLeft + 1) * bwt->saInterval < initialSaIndexRight) {
        stage3EndSaIndexRight = (stage2EndSaExplicitIndexLeft + 1) * bwt->saInterval;
    } else {
        stage3EndSaIndexRight = initialSaIndexRight + 1;
    }

    // Search for starting SA index

    saIndexLeft = stage3StartSaIndexLeft;
    saIndexRight = stage3StartSaIndexRight;
    llcp = stage2Startllcp;
    rlcp = stage2Startrlcp;
    tempKey[numOfWord - 1] = tempKeyATrailing;

    mlcp = 0;
    while ((saIndexLeft + 1) < saIndexRight) {

        saIndexMiddle = average(saIndexLeft, saIndexRight);

        saValue = BWTSaValue(bwt, saIndexMiddle) + cachedNumOfChar;
        shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
        index = saValue / CHAR_PER_WORD;

        // Try to increase mlcp
        mlcp = min(llcp, rlcp);        // mlcp = the characters (in unit of 16 for DNA) matched so far
        
        do  {
            if (shift != 0) {
                text = (packedText[index + mlcp] << shift) | (packedText[index + mlcp + 1] >> (BITS_IN_WORD - shift));
            } else {
                text = packedText[index + mlcp];
            }
        } while (tempKey[mlcp] == text && ++mlcp < numOfWord);

        if (mlcp < numOfWord && tempKey[mlcp] > text) {
            saIndexLeft = saIndexMiddle;
            llcp = mlcp;
        } else {
            saIndexRight = saIndexMiddle;
            rlcp = mlcp;
        }
    }

    stage3SaIndexLeft = saIndexRight;
    
    // Search for ending SA index

    saIndexLeft = stage3EndSaIndexLeft;
    saIndexRight = stage3EndSaIndexRight;
    llcp = stage2Endllcp;
    rlcp = stage2Endrlcp;
    tempKey[numOfWord - 1] = tempKeyTTrailing;

    mlcp = 0;
    while ((saIndexLeft + 1) < saIndexRight) {

        saIndexMiddle = average(saIndexLeft, saIndexRight);

        saValue = BWTSaValue(bwt, saIndexMiddle) + cachedNumOfChar;
        shift = BIT_PER_CHAR * (saValue % CHAR_PER_WORD);
        index = saValue / CHAR_PER_WORD;

        // Try to increase mlcp
        mlcp = min(llcp, rlcp);        // mlcp = the characters (in unit of 16 for DNA) matched so far
        
        do {
            if (shift != 0) {
                text = (packedText[index + mlcp] << shift) | (packedText[index + mlcp + 1] >> (BITS_IN_WORD - shift));
            } else {
                text = packedText[index + mlcp];
            }
        } while (tempKey[mlcp] == text && ++mlcp < numOfWord);
        if (mlcp >= numOfWord || tempKey[mlcp] >= text) {
            saIndexLeft = saIndexMiddle;
            llcp = mlcp;
        } else {
            saIndexRight = saIndexMiddle;
            rlcp = mlcp;
        }
    }

    stage3SaIndexRight = saIndexLeft;

    *resultSaIndexLeft = stage3SaIndexLeft;
    *resultSaIndexRight = stage3SaIndexRight;

    return (stage3SaIndexLeft <= stage3SaIndexRight);

}

int BWTBackwardSearch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, 
                      unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight) {

    unsigned long long pos;
    unsigned int c;
    unsigned long long initialSaRangeIndex;
    unsigned long long direction1, direction2;
    unsigned long long index1, index2;
    unsigned long long occExplicitIndex1, occExplicitIndex2;
    unsigned long long occIndex1, occIndex2;
    unsigned long long tempExplicit1, tempExplicit2;
    unsigned long long tempOccValue1, tempOccValue2;
    unsigned long long estimatedIndex1, estimatedIndex2;
    unsigned long long estimatedOccExplicitIndex1, estimatedOccExplicitIndex2;
    unsigned long long decodeValue;
    unsigned long long cachedNumOfChar;

    cachedNumOfChar = min(bwt->cachedSaIndexNumOfChar, keyLength);

    initialSaRangeIndex = 0;
    for (pos = keyLength; pos > keyLength - cachedNumOfChar; pos--) {
        initialSaRangeIndex |= convertedKey[pos-1] << ((keyLength - pos) * BIT_PER_CHAR);
    }

    index1 = bwt->cachedSaIndex[initialSaRangeIndex << ((bwt->cachedSaIndexNumOfChar - cachedNumOfChar) * BIT_PER_CHAR)];
    index2 = bwt->cachedSaIndex[(initialSaRangeIndex + 1) << ((bwt->cachedSaIndexNumOfChar - cachedNumOfChar) * BIT_PER_CHAR)];

    for (; pos > 0 && index1 < index2; pos--) {

        c = convertedKey[pos-1];
        
        // $ is supposed to be positioned at inverseSa0 but it is not encoded
        // therefore index is subtracted by 1 for adjustment
        // Note that this adjustment is not done when BWTOccValue is used
        index1 -= (index1 > bwt->inverseSa0);
        index2 -= (index2 > bwt->inverseSa0);

        // Calculate index to explicit occurrence
        occExplicitIndex1 = (index1 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;    // Bidirectional encoding
        occExplicitIndex2 = (index2 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;    // Bidirectional encoding
        occIndex1 = occExplicitIndex1 * OCC_INTERVAL;
        occIndex2 = occExplicitIndex2 * OCC_INTERVAL;

        direction1 = -(occIndex1 > index1);
        direction2 = -(occIndex2 > index2);

        tempExplicit1 = BWTOccValueExplicit(bwt, occExplicitIndex1, c);
        tempExplicit2 = BWTOccValueExplicit(bwt, occExplicitIndex2, c);

        // Estimate the SA index before BWT is decoded

        estimatedIndex1 = bwt->cumulativeFreq[c] + tempExplicit1 + (ESTIMATED_OCC_DIFF & ~direction1) - (ESTIMATED_OCC_DIFF & direction1) + 1;
        estimatedIndex2 = bwt->cumulativeFreq[c] + tempExplicit2 + (ESTIMATED_OCC_DIFF & ~direction2) - (ESTIMATED_OCC_DIFF & direction2) + 1;
        estimatedOccExplicitIndex1 = (estimatedIndex1 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;    // Bidirectional encoding
        estimatedOccExplicitIndex2 = (estimatedIndex2 + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;    // Bidirectional encoding

        // Pre-fetch memory to be accessed
        BWTPrefetchBWT(bwt, estimatedIndex1);
        BWTPrefetchBWT(bwt, estimatedIndex2);

        // Pre-fetch memory to be accessed
        BWTPrefetchOccValueExplicit(bwt, estimatedOccExplicitIndex1);
        BWTPrefetchOccValueExplicit(bwt, estimatedOccExplicitIndex2);

        // Decode BWT
        if (occIndex1 != index1) {
            decodeValue = BWTDecode(bwt, occIndex1, index1, c);
            tempOccValue1 = (decodeValue & ~direction1) - (decodeValue & direction1);
        } else {
            tempOccValue1 = 0;
        }

        if (occIndex2 != index2) {
            decodeValue = BWTDecode(bwt, occIndex2, index2, c);
            tempOccValue2 = (decodeValue & ~direction2) - (decodeValue & direction2);
        } else {
            tempOccValue2 = 0;
        }

        index1 = bwt->cumulativeFreq[c] + tempExplicit1 + tempOccValue1 + 1;
        index2 = bwt->cumulativeFreq[c] + tempExplicit2 + tempOccValue2 + 1;

    }

    *resultSaIndexLeft = index1;
    *resultSaIndexRight = index2 - 1;

    return (index1 < index2);

}

int BWTBackwardSearchCheckWithText(const unsigned char *convertedKey, const unsigned int *packedKey, const unsigned int keyLength,
                              const BWT *bwt, const unsigned int *packedText, const unsigned long long textCheckingCostFactor,
                              const unsigned long long maxnumOfTextPosition, HitList* __restrict hitList,
                              unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight) {

    unsigned long long startSaIndex, endSaIndex;
    unsigned long long pos;
    unsigned int c;
    unsigned long long i;
    unsigned long long numOfMatch;
    unsigned long long textPosition;
    unsigned long long wordMatched, wordToMatch;
    unsigned long long shift, index;

    pos = keyLength - 1;
    c = convertedKey[pos];

    #ifdef DEBUG
    if (c >= ALPHABET_SIZE) {
        fprintf(stderr, "BWTBackwardSearchWithText() : invalid key!\n");
        exit(1);
    }
    #endif

    startSaIndex = bwt->cumulativeFreq[c] + 1;
    endSaIndex = bwt->cumulativeFreq[c + 1];

    if (startSaIndex > endSaIndex) {
        // The last character of search pattern does not exists in text
        return 0;
    }

    numOfMatch = endSaIndex - startSaIndex + 1;

    while (pos >= 1 && startSaIndex <= endSaIndex &&        // Search result not determined yet
           !(pos >= numOfMatch * (bwt->saInterval / 2 + pos / textCheckingCostFactor) &&            // Text checking not cheaper
               numOfMatch <= maxnumOfTextPosition && pos + BITS_IN_WORD <= keyLength)) {
                                                             // pos + BITS_IN_WORD <= keyLength -> masking of packed text is not needed
        c = convertedKey[pos - 1];
        startSaIndex = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex, c) + 1;
        endSaIndex = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex + 1, c);
        numOfMatch = endSaIndex - startSaIndex + 1;
        pos--;

    }

    if (pos >= 1 && startSaIndex <= endSaIndex) {
        // Use text checking
        endSaIndex = 0;
        wordToMatch = (pos + CHAR_PER_WORD - 1) / CHAR_PER_WORD;
        for (i=0; i<numOfMatch; i++) {
            textPosition = BWTSaValue(bwt, startSaIndex + i);
            if (textPosition >= pos) {
                textPosition -= pos;
                // check text
                shift = BIT_PER_CHAR * (textPosition % CHAR_PER_WORD);
                index = textPosition / CHAR_PER_WORD;
                wordMatched = 0;
                while (wordMatched < wordToMatch && 
                        packedKey[wordMatched] == BWTGetWordPackedText(packedText, index + wordMatched, 
                                                    shift, BITS_IN_WORD)) {
                    wordMatched++;
                }
                if (wordMatched == wordToMatch) {
                    hitList[endSaIndex].posText = textPosition;
                    endSaIndex++;
                }
            }
        }
        startSaIndex = 1;        // saIndex cannot be returned when text checking is used, so only return the number of occurrences
    }

    *resultSaIndexLeft = startSaIndex;
    *resultSaIndexRight = endSaIndex;

    // Number of occurrence = endSaIndex - startSaIndex + 1
    return (startSaIndex <= endSaIndex);

}


unsigned long long BWTHammingDistMaxSaIndexGroup(const unsigned int keyLength, const unsigned int maxError) {

    unsigned long long a, c, d, e;
    unsigned long long t;

    t = 1;
    a = 1;
    c = 1;
    d = 1;

    for (e=1; e<=maxError; e++) {
        // For error = e, combination = mCe (alphabetSize - 1) ^ e, where m = keyLength
        c *= keyLength + 1 - e;
        d *= e;
        a *= ALPHABET_SIZE - 1;
        t += c / d * a;
    }

    return t;

}

unsigned long long BWTHammingDistCountOcc(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError) {

    // stack
    unsigned long long startSaIndex[MAX_ARPROX_MATCH_LENGTH + 1];
    unsigned long long endSaIndex[MAX_ARPROX_MATCH_LENGTH];
    unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
    unsigned long long errorPos[MAX_APPROX_MATCH_ERROR + 1];

    unsigned long long numOfError;

    unsigned int c, e;
    unsigned long long errorAdded;
    unsigned long long numOfHit = 0;

    #ifdef DEBUG
    if (maxError > keyLength) {
        fprintf(stderr, "BWTHammingDistCountOcc() : maxError > keyLength!\n");
        exit(1);
    }
    if (keyLength > MAX_ARPROX_MATCH_LENGTH) {
        fprintf(stderr, "BWTHammingDistCountOcc() : keyLength > MAX_ARPROX_MATCH_LENGTH!\n");
        exit(1);
    }
    if (maxError > MAX_APPROX_MATCH_ERROR) {
        fprintf(stderr, "BWTHammingDistCountOcc() : maxError > MAX_APPROX_MATCH_ERROR!\n");
        exit(1);
    }
    #endif

    // Set this boundary case so that the last character of generated pattern can be located by backward search
    startSaIndex[keyLength] = 0;
    endSaIndex[keyLength] = bwt->textLength;

    // Exact match
    BWTBackwardSearch(convertedKey, keyLength, bwt, startSaIndex, endSaIndex);
    if (endSaIndex[0] >= startSaIndex[0]) {;
        numOfHit += endSaIndex[0] - startSaIndex[0] + 1;
    }

    // With errors
    for (numOfError = 1; numOfError <= maxError; numOfError++) {

        memcpy(generatedPattern, convertedKey, keyLength);

        // Set initial state
        errorAdded = 1;
        e = keyLength - 1;
        while (e>0) {
            c = generatedPattern[e];
             startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
            endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
            e--;
        }
        errorPos[1] = e;
        generatedPattern[e] = (unsigned char)-1;

        while (errorAdded > 0) {
            
            // Increment the last error
            e = errorPos[errorAdded];
            generatedPattern[e]++;
            if (generatedPattern[e] == convertedKey[e]) {
                generatedPattern[e]++;
            }
            if (generatedPattern[e] >= ALPHABET_SIZE) {
                // Cannot increment the last error; try moving the error to the left
                if (e > numOfError - errorAdded) {
                    c = generatedPattern[e] = convertedKey[e];
                    startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                    endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
                    if (startSaIndex[e] >= endSaIndex[e]) {
                        // Pattern suffix not exists in text
                        errorAdded--;
                        continue;
                    }
                    errorPos[errorAdded]--;
                    e = errorPos[errorAdded];
                    generatedPattern[e] = (convertedKey[e] == 0);    // 1 if convertedKey = 0, 0 otherwise
                } else {
                    // Cannot move error to the left
                    generatedPattern[e] = convertedKey[e];
                    errorAdded--;
                    continue;
                }
            }

            c = generatedPattern[e];
             startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
            endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);

            while (errorAdded < numOfError && startSaIndex[e] <= endSaIndex[e]) {
                // Add errors
                errorAdded++;
                errorPos[errorAdded] = errorPos[errorAdded-1] - 1;
                e = errorPos[errorAdded];
                c = generatedPattern[e] = (convertedKey[e] == 0);    // 1 if convertedKey = 0, 0 otherwise
                startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
            }

            // Search for pattern
            while (e>0 && startSaIndex[e] <= endSaIndex[e]) {
                e--;
                c = generatedPattern[e];
                startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
            }
            if (startSaIndex[e] <= endSaIndex[e]) {
                numOfHit += endSaIndex[0] - startSaIndex[0] + 1;
            }

        }

    }

    return numOfHit;

}
/*
unsigned int BWTHammingDistMatch(const BWT *bwt, const unsigned char *convertedKey, const HitCombination *hitCombination,
                        const SaIndexRange *cachedSaIndex, const unsigned int cachedSaIndexNumOfChar,
                        SaIndexGroupNew* __restrict saIndexGroup, const unsigned int maxSaIndexGroup) {

    // stack
    unsigned int startSaIndex[MAX_ARPROX_MATCH_LENGTH][ALPHABET_SIZE];
    unsigned int endSaIndex[MAX_ARPROX_MATCH_LENGTH][ALPHABET_SIZE];
    unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
    unsigned int errorPos[MAX_APPROX_MATCH_ERROR + 1];

    unsigned int numOfSaGroup;
    unsigned int numOfError;

    unsigned int i, c, e;
    unsigned int errorAdded;
    unsigned int errorVector;

    unsigned int tempSaIndexLeft, tempSaIndexRight;

    unsigned int saRangeIndex;
    unsigned int initialSaRangeIndex;
    unsigned int minE;
    unsigned int mask[16] = { 0xFFFFFFFC, 0xFFFFFFF3, 0xFFFFFFCF, 0xFFFFFF3F,
                     0xFFFFFCFF, 0xFFFFF3FF, 0xFFFFCFFF, 0xFFFF3FFF,
                     0xFFFCFFFF, 0xFFF3FFFF, 0xFFCFFFFF, 0xFF3FFFFF,
                     0xFCFFFFFF, 0xF3FFFFFF, 0xCFFFFFFF, 0x3FFFFFFF };

    #ifdef DEBUG
    if (hitCombination->maxError > hitCombination->keyLength) {
        fprintf(stderr, "BWTHammingDistMatch() : maxError > patternLength!\n");
        exit(1);
    }
    if (hitCombination->keyLength > MAX_ARPROX_MATCH_LENGTH) {
        fprintf(stderr, "BWTHammingDistMatch() : patternLength > MAX_ARPROX_MATCH_LENGTH!\n");
        exit(1);
    }
    if (hitCombination->maxError > MAX_APPROX_MATCH_ERROR) {
        fprintf(stderr, "BWTHammingDistMatch() : maxError > MAX_APPROX_MATCH_ERROR!\n");
        exit(1);
    }
    #endif

    // Setup for looking up SA index range
    initialSaRangeIndex = 0;
    if (cachedSaIndexNumOfChar <= hitCombination->keyLength) {
        for (i=0; i<cachedSaIndexNumOfChar; i++) {
            initialSaRangeIndex |= convertedKey[hitCombination->keyLength - 1 - i] << ((cachedSaIndexNumOfChar - 1 - i) * BIT_PER_CHAR);
        }
        minE = hitCombination->keyLength - 1 - cachedSaIndexNumOfChar;
    } else {
        minE = hitCombination->keyLength;
    }

    // Exact match
    for (i=hitCombination->keyLength; i>0; i++) {
    }


    // Set this boundary case so that the last character of generated pattern can be located by backward search
    for (i=0; i<ALPHABET_SIZE; i++) {
        startSaIndex[hitCombination->keyLength][i] = bwt->cumulativeFreq[i] + 1;
        endSaIndex[hitCombination->keyLength][i] = bwt->cumulativeFreq[i+1];
    }

    // Exact match
    numOfSaGroup = BWTBackwardSearch(convertedKey, keyLength, bwt, &tempSaIndexLeft, &tempSaIndexRight);
    if (numOfSaGroup > 0) {
        saIndexGroup[0].startSaIndex = tempSaIndexLeft;
        saIndexGroup[0].numOfMatch = tempSaIndexRight - tempSaIndexLeft + 1;
        saIndexGroup[0].info = 0;
    }

    // With errors
    for (numOfError = 1; numOfError <= maxError; numOfError++) {

        memcpy(generatedPattern, convertedKey, keyLength);
        saRangeIndex = initialSaRangeIndex;

        // Set initial state
        errorAdded = 1;
        e = keyLength - 1;
        errorVector = FIRST_BIT_MASK >> (keyLength - 1);
        while (e>0 && (errorVector & matchBitVector) != 0) {
            c = generatedPattern[e];

            if (e > minE) {
                saRangeIndex &= mask[e - minE - 1];
                saRangeIndex |= c << ((e - minE - 1) * BIT_PER_CHAR);
            } else {
                if (e < minE) {
                    startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                    endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
                } else {
                    startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].startSaIndex, c) + 1;
                    endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].endSaIndex + 1, c);
                }
            }

            e--;
            errorVector <<= 1;
        }
        if ((errorVector & matchBitVector) == 0) {
            errorPos[1] = e;
            generatedPattern[e] = (unsigned char)-1;
        } else {
            // cannot even place the first error
            return numOfSaGroup;
        }

        while (errorAdded > 0) {
            
            // Increment the last error
            e = errorPos[errorAdded];
            generatedPattern[e]++;
            if (generatedPattern[e] == convertedKey[e]) {
                generatedPattern[e]++;
            }
            if (generatedPattern[e] >= ALPHABET_SIZE) {
                // Cannot increment the last error; try moving the error to the left
                errorVector &= ALL_ONE_MASK >> (e + 1);        // clear the error vector from the error position to the left
                if (e > numOfError - errorAdded) {
                    c = generatedPattern[e] = convertedKey[e];

                    if (e > minE) {
                        saRangeIndex &= mask[e - minE - 1];
                        saRangeIndex |= c << ((e - minE - 1) * BIT_PER_CHAR);
                    } else {
                        if (e < minE) {
                            startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                            endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
                        } else {
                            startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].startSaIndex, c) + 1;
                            endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].endSaIndex + 1, c);
                        }
                        if (startSaIndex[e] >= endSaIndex[e]) {
                            // Pattern suffix not exists in text
                            errorAdded--;
                            continue;
                        }
                    }

                    errorPos[errorAdded]--;
                    e = errorPos[errorAdded];
                    errorVector |= FIRST_BIT_MASK >> e;        // mark the error position with 1 in error vector
                    if ((errorVector & matchBitVector) == 0) {
                        generatedPattern[e] = (convertedKey[e] == 0);    // 1 if convertedKey = 0, 0 otherwise
                    } else {
                        // error in position for match only
                        generatedPattern[e] = (unsigned char)ALPHABET_SIZE - 1;
                        continue;
                    }
                } else {
                    // Cannot move error to the left
                    generatedPattern[e] = convertedKey[e];
                    errorAdded--;
                    continue;
                }
            }

            c = generatedPattern[e];

            if (e > minE) {
                saRangeIndex &= mask[e - minE - 1];
                saRangeIndex |= c << ((e - minE - 1) * BIT_PER_CHAR);
            } else {
                if (e < minE) {
                    startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                    endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
                } else {
                    startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].startSaIndex, c) + 1;
                    endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].endSaIndex + 1, c);
                }
            }

            while (errorAdded < numOfError && (e > minE || startSaIndex[e] <= endSaIndex[e])) {
                // Add errors
                errorAdded++;
                errorPos[errorAdded] = errorPos[errorAdded-1] - 1;
                e = errorPos[errorAdded];
                errorVector |= FIRST_BIT_MASK >> e;        // mark the error position with 1 in error vector
                if ((errorVector & matchBitVector) == 0) {
                    c = generatedPattern[e] = (convertedKey[e] == 0);    // 1 if convertedKey = 0, 0 otherwise
                } else {
                    // error in position for match only
                    generatedPattern[e] = (unsigned char)ALPHABET_SIZE - 1;
                    // so that it skip the remaining of the while loop
                    startSaIndex[e] = 1;        
                    endSaIndex[e] = 0;
                    break;
                }

                if (e > minE) {
                    saRangeIndex &= mask[e - minE - 1];
                    saRangeIndex |= c << ((e - minE - 1) * BIT_PER_CHAR);
                } else {
                    if (e < minE) {
                        startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                        endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
                    } else {
                        startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].startSaIndex, c) + 1;
                        endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].endSaIndex + 1, c);
                    }
                }

            }

            // Search for pattern
            while (e>0 && (e > minE || startSaIndex[e] <= endSaIndex[e])) {
                e--;
                c = generatedPattern[e];

                if (e > minE) {
                    saRangeIndex &= mask[e - minE - 1];
                    saRangeIndex |= c << ((e - minE - 1) * BIT_PER_CHAR);
                } else {
                    if (e < minE) {
                        startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                        endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
                    } else {
                        startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].startSaIndex, c) + 1;
                        endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].endSaIndex + 1, c);
                    }
                }

            }
            if (startSaIndex[e] <= endSaIndex[e]) {
                if (numOfSaGroup < maxSaIndexGroup) {
                    saIndexGroup[numOfSaGroup].startSaIndex = startSaIndex[0];
                    saIndexGroup[numOfSaGroup].numOfMatch = endSaIndex[0] - startSaIndex[0] + 1;
                    saIndexGroup[numOfSaGroup].info = errorVector;
                    numOfSaGroup++;
                } else {
                    fprintf(stderr, "Not enough memory to store all SA index groups!\n");
                    numOfError = maxError;    // to exit the for loop
                    break;
                }
            }

        }

    }

    QSort(saIndexGroup, numOfSaGroup, sizeof(SaIndexGroupNew), SaIndexGroupStartSaIndexOrder);

    return numOfSaGroup;

}
*/

// This is breath-first hamming distance search
// It does not check whether there is enough memory to hold all the saIndexGroup
// Caller needs to invoke BWTHammingDistMaxSaIndexGroup() to determine the maximum no. of saIndexGroup
// maxError must be at least 1
/*
int BWTHammingDistMatchBF(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
                          SaIndexGroupNew* __restrict saIndexGroup, const unsigned int posQuery, const unsigned int info) {

    unsigned int cachedSaNumOfChar;
    unsigned int temp;
    unsigned int pos;

    #ifdef DEBUG
    if (maxError < 1) {
        fprintf(stderr, "BWTHammingDistMatchBF() : maxError must be at least 1!\n");
        exit(1);
    }
    if (maxError > keyLength) {
        fprintf(stderr, "BWTHammingDistMatch() : maxError > keyLength!\n");
        exit(1);
    }
    #endif

    // Setup for looking up cached SA index range
    if (bwt->cachedSaIndexNumOfChar <= keyLength) {
        cachedSaNumOfChar = bwt->cachedSaIndexNumOfChar;
    } else {
        cachedSaNumOfChar = 1;
    }

    ((SaIndexGroupTemp*)saIndexGroup + 0)->startSaIndex1 = 0;
    ((SaIndexGroupTemp*)saIndexGroup + 1)->startSaIndex1 = 1;
    ((SaIndexGroupTemp*)saIndexGroup + 2)->startSaIndex1 = 2;
    ((SaIndexGroupTemp*)saIndexGroup + 3)->startSaIndex1 = 3;

    // Swap exact match to the last entry
    temp = ((SaIndexGroupTemp*)saIndexGroup + convertedKey[keyLength-1])->startSaIndex1;
    ((SaIndexGroupTemp*)saIndexGroup + 0)->startSaIndex1 = ((SaIndexGroupTemp*)saIndexGroup + convertedKey[keyLength-1])->startSaIndex1;
    ((SaIndexGroupTemp*)saIndexGroup + convertedKey[keyLength-1])->startSaIndex1 = temp;

    pos = 1;



}
*/

int BWTHammingDistMatchOld(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
                           SaIndexGroupNew* __restrict saIndexGroup, const unsigned long long maxSaIndexGroup, 
                           const unsigned long long posQuery, const unsigned long long info) {

    // stack
    unsigned long long startSaIndex[MAX_ARPROX_MATCH_LENGTH + 1];
    unsigned long long endSaIndex[MAX_ARPROX_MATCH_LENGTH];
    unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH];
    unsigned long long errorPos[MAX_APPROX_MATCH_ERROR + 1];

    unsigned long long numOfSaGroup;
    unsigned long long numOfError;

    unsigned long long i;
    unsigned c;
    unsigned long long e;
    unsigned long long errorAdded;

    unsigned long long tempSaIndex[2];

    unsigned long long saRangeIndex;
    unsigned long long initialSaRangeIndex;
    unsigned long long minE;
    unsigned int mask[16] = { 0xFFFFFFFC, 0xFFFFFFF3, 0xFFFFFFCF, 0xFFFFFF3F,
                     0xFFFFFCFF, 0xFFFFF3FF, 0xFFFFCFFF, 0xFFFF3FFF,
                     0xFFFCFFFF, 0xFFF3FFFF, 0xFFCFFFFF, 0xFF3FFFFF,
                     0xFCFFFFFF, 0xF3FFFFFF, 0xCFFFFFFF, 0x3FFFFFFF };

    #ifdef DEBUG
    if (maxError > keyLength) {
        fprintf(stderr, "BWTHammingDistMatch() : maxError > keyLength!\n");
        exit(1);
    }
    if (keyLength > MAX_ARPROX_MATCH_LENGTH) {
        fprintf(stderr, "BWTHammingDistMatch() : keyLength > MAX_ARPROX_MATCH_LENGTH!\n");
        exit(1);
    }
    if (maxError > MAX_APPROX_MATCH_ERROR) {
        fprintf(stderr, "BWTHammingDistMatch() : maxError > MAX_APPROX_MATCH_ERROR!\n");
        exit(1);
    }
    #endif

    // Setup for looking up SA index range
    initialSaRangeIndex = 0;
    if (bwt->cachedSaIndexNumOfChar <= keyLength && bwt->cachedSaIndexNumOfChar != 0) {
        for (i = keyLength; i > keyLength - bwt->cachedSaIndexNumOfChar; i--) {
            initialSaRangeIndex |= convertedKey[i-1] << (2 * (keyLength - i));
        }
        minE = keyLength - 1 - bwt->cachedSaIndexNumOfChar;
    } else {
        minE = keyLength;
    }


    // Set this boundary case so that the last character of generated pattern can be located by backward search
    startSaIndex[keyLength] = 0;
    endSaIndex[keyLength] = bwt->textLength;

    // Exact match
    numOfSaGroup = BWTBackwardSearch(convertedKey, keyLength, bwt, tempSaIndex, tempSaIndex + 1);
    if (numOfSaGroup > 0) {
        if (numOfSaGroup <= maxSaIndexGroup) {
            saIndexGroup[0].startSaIndex = tempSaIndex[0];
            saIndexGroup[0].numOfMatch = tempSaIndex[1] - tempSaIndex[0] + 1;
            saIndexGroup[0].posQuery = posQuery;
            saIndexGroup[0].info = info;
        } else {
            return -1;
        }
    }

    // With errors
    for (numOfError = 1; numOfError <= maxError; numOfError++) {

        memcpy(generatedPattern, convertedKey, keyLength);
        saRangeIndex = initialSaRangeIndex;

        // Set initial state
        errorAdded = 1;
        e = keyLength - 1;
        errorPos[1] = e;
        generatedPattern[e] = (unsigned char)-1;

        while (errorAdded > 0) {
            
            // Increment the last error
            e = errorPos[errorAdded];
            generatedPattern[e]++;
            if (generatedPattern[e] == convertedKey[e]) {
                generatedPattern[e]++;
            }
            if (generatedPattern[e] >= ALPHABET_SIZE) {
                // Cannot increment the last error; try moving the error to the left
                if (e > numOfError - errorAdded) {
                    c = generatedPattern[e] = convertedKey[e];

                    if (e > minE) {
                        saRangeIndex &= mask[keyLength - 1 - e];
                        saRangeIndex |= c << ((keyLength - 1 - e) * BIT_PER_CHAR);
                    } else {
                        if (e < minE) {
                            BWTOccValueTwoIndex(bwt, startSaIndex[e+1], endSaIndex[e+1] + 1, c, tempSaIndex);
                        } else {
                            BWTOccValueTwoIndex(bwt, bwt->cachedSaIndex[saRangeIndex], bwt->cachedSaIndex[saRangeIndex + 1], c, tempSaIndex);
                        }
                        startSaIndex[e] = bwt->cumulativeFreq[c] + tempSaIndex[0] + 1;
                        endSaIndex[e] = bwt->cumulativeFreq[c] + tempSaIndex[1];
                        if (startSaIndex[e] >= endSaIndex[e]) {
                            // Pattern suffix not exists in text
                            errorAdded--;
                            continue;
                        }
                    }

                    errorPos[errorAdded]--;
                    e = errorPos[errorAdded];
                    generatedPattern[e] = (convertedKey[e] == 0);    // 1 if convertedKey = 0, 0 otherwise
                } else {
                    // Cannot move error to the left
                    generatedPattern[e] = convertedKey[e];
                    errorAdded--;
                    continue;
                }
            }

            c = generatedPattern[e];

            if (e > minE) {
                saRangeIndex &= mask[keyLength - 1 - e];
                saRangeIndex |= c << ((keyLength - 1 - e) * BIT_PER_CHAR);
            } else {
                if (e < minE) {
                    BWTOccValueTwoIndex(bwt, startSaIndex[e+1], endSaIndex[e+1] + 1, c, tempSaIndex);
                } else {
                    BWTOccValueTwoIndex(bwt, bwt->cachedSaIndex[saRangeIndex], bwt->cachedSaIndex[saRangeIndex + 1], c, tempSaIndex);
                }
                startSaIndex[e] = bwt->cumulativeFreq[c] + tempSaIndex[0] + 1;
                endSaIndex[e] = bwt->cumulativeFreq[c] + tempSaIndex[1];
            }

            while (errorAdded < numOfError && e > 0 && (e > minE || startSaIndex[e] <= endSaIndex[e])) {
                // Add errors
                errorAdded++;
                errorPos[errorAdded] = errorPos[errorAdded-1] - 1;
                e = errorPos[errorAdded];
                c = generatedPattern[e] = (convertedKey[e] == 0);    // 1 if convertedKey = 0, 0 otherwise

                if (e > minE) {
                    saRangeIndex &= mask[keyLength - 1 - e];
                    saRangeIndex |= c << ((keyLength - 1 - e) * BIT_PER_CHAR);
                } else {
                    if (e < minE) {
                        BWTOccValueTwoIndex(bwt, startSaIndex[e+1], endSaIndex[e+1] + 1, c, tempSaIndex);
                    } else {
                        BWTOccValueTwoIndex(bwt, bwt->cachedSaIndex[saRangeIndex], bwt->cachedSaIndex[saRangeIndex + 1], c, tempSaIndex);
                    }
                    startSaIndex[e] = bwt->cumulativeFreq[c] + tempSaIndex[0] + 1;
                    endSaIndex[e] = bwt->cumulativeFreq[c] + tempSaIndex[1];
                }

            }

            // Search for pattern
            while (e>0 && (e > minE || startSaIndex[e] <= endSaIndex[e])) {
                e--;
                c = generatedPattern[e];

                if (e > minE) {
                    saRangeIndex &= mask[keyLength - 1 - e];
                    saRangeIndex |= c << ((keyLength - 1 - e) * BIT_PER_CHAR);
                } else {
                    if (e < minE) {
                        BWTOccValueTwoIndex(bwt, startSaIndex[e+1], endSaIndex[e+1] + 1, c, tempSaIndex);
                    } else {
                        BWTOccValueTwoIndex(bwt, bwt->cachedSaIndex[saRangeIndex], bwt->cachedSaIndex[saRangeIndex + 1], c, tempSaIndex);
                    }
                    startSaIndex[e] = bwt->cumulativeFreq[c] + tempSaIndex[0] + 1;
                    endSaIndex[e] = bwt->cumulativeFreq[c] + tempSaIndex[1];
                }

            }
            if (startSaIndex[e] <= endSaIndex[e]) {
                if (numOfSaGroup < maxSaIndexGroup) {
                    saIndexGroup[numOfSaGroup].startSaIndex = startSaIndex[0];
                    saIndexGroup[numOfSaGroup].numOfMatch = endSaIndex[0] - startSaIndex[0] + 1;
                    saIndexGroup[numOfSaGroup].posQuery = posQuery;
                    saIndexGroup[numOfSaGroup].info = info;
                    numOfSaGroup++;
                } else {
                    return -1;
                }
            }

        }

    }

    return numOfSaGroup;

}

unsigned long long BWTSubPatternHammingDistCountOcc(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
                                     const unsigned int maxError, const unsigned long long skip) {

    unsigned long long i;
    unsigned long long numOfSubPattern, advance;

    unsigned long long numOfHit;
    unsigned long long nextFilteredChar;
    
    // Determine the number of sub-patterns
    advance = skip + 1;
    numOfSubPattern = (keyLength - subPatternLength + advance) / advance;

    numOfHit = 0;
    
    // Find SA index groups for all sub-patterns

    for (nextFilteredChar = keyLength; nextFilteredChar > 0 && convertedKey[nextFilteredChar-1] < ALPHABET_SIZE; nextFilteredChar--) {
    }

    for (i=0; i<numOfSubPattern; i++) {

        if (nextFilteredChar > keyLength - i*advance) {
            // Advance to the next 'N'
            for (; nextFilteredChar > keyLength - i*advance && convertedKey[nextFilteredChar-1] >= ALPHABET_SIZE; nextFilteredChar--) {
            }
            for (; nextFilteredChar > 0 && convertedKey[nextFilteredChar-1] < ALPHABET_SIZE; nextFilteredChar--) {
            }
        }

        if (keyLength - i * advance - subPatternLength >= nextFilteredChar) {
            numOfHit += BWTHammingDistCountOcc(convertedKey + keyLength - i*advance - subPatternLength, subPatternLength,
                                                  bwt, maxError);
        }

    }

    return numOfHit;


}
/*
int BWTSubPatternHammingDistSaIndex(const BWT *bwt, const unsigned char *convertedKey, const int keyLength, const int skip, 
                                    const HitCombination *hitCombination, 
                                    const SaIndexRange *cachedSaIndex, const unsigned int cachedSaIndexNumOfChar,
                                    SaIndexGroupNew* __restrict saIndexGroup, const int maxnumOfSaIndexGroup,
                                    int* __restrict firstSaIndexGroupForSubPattern) {

    int advance, numOfSubPattern;
    int nextFilteredChar;
    int i;
    int numOfSaGroup;

    // Determine the number of sub-patterns
    advance = skip + 1;
    numOfSubPattern = (keyLength - hitCombination->keyLength + advance) / advance;

    for (nextFilteredChar = keyLength; nextFilteredChar > 0 && convertedKey[nextFilteredChar-1] < ALPHABET_SIZE; nextFilteredChar--) {
    }

    firstSaIndexGroupForSubPattern[0] = 0;
    for (i=0; i<numOfSubPattern; i++) {

        if (nextFilteredChar > keyLength - i*advance) {
            // Advance to the next 'N'
            for (; nextFilteredChar > keyLength - i*advance && convertedKey[nextFilteredChar-1] >= ALPHABET_SIZE; nextFilteredChar--) {
            }
            for (; nextFilteredChar > 0 && convertedKey[nextFilteredChar-1] < ALPHABET_SIZE; nextFilteredChar--) {
            }
        }

        if (keyLength - i * advance - hitCombination->keyLength >= nextFilteredChar) {

            numOfSaGroup = BWTHammingDistMatch(bwt, convertedKey + keyLength - i*advance - hitCombination->keyLength, hitCombination,
                                                  cachedSaIndex, cachedSaIndexNumOfChar,
                                                  saIndexGroup + firstSaIndexGroupForSubPattern[i], maxnumOfSaIndexGroup - firstSaIndexGroupForSubPattern[i]);
            firstSaIndexGroupForSubPattern[i+1] = firstSaIndexGroupForSubPattern[i] + numOfSaGroup;

        } else {
            // sub-pattern filtered
            firstSaIndexGroupForSubPattern[i+1] = firstSaIndexGroupForSubPattern[i];
        }

    }
    return firstSaIndexGroupForSubPattern[numOfSubPattern];

}
*/

int BWTSubPatternHammingDistSaIndexOld(const unsigned char *convertedKey, const int keyLength, const int subPatternLength, const BWT *bwt, 
                                    const int maxError, const int skip, const int lengthProcessed, int* __restrict lengthInCurrentRound,
                                    SaIndexGroupNew* __restrict saIndexGroup, const int maxnumOfSaIndexGroup) {

    int advance, numOfSubPattern;
    int nextFilteredChar;
    int i;
    int totalnumOfSaIndexGroup = 0;
    int numOfSaIndexGroup;
    int effectiveKeyLength;

    // Determine the number of sub-patterns
    effectiveKeyLength = keyLength - lengthProcessed;
    advance = skip + 1;
    numOfSubPattern = (effectiveKeyLength - subPatternLength + advance) / advance;

    for (nextFilteredChar = effectiveKeyLength; nextFilteredChar > 0 && convertedKey[nextFilteredChar-1] < ALPHABET_SIZE; nextFilteredChar--) {
    }

    for (i=0; i<numOfSubPattern; i++) {

        if (nextFilteredChar > effectiveKeyLength - i*advance) {
            // Advance to the next 'N'
            for (; nextFilteredChar > effectiveKeyLength - i*advance && convertedKey[nextFilteredChar-1] >= ALPHABET_SIZE; nextFilteredChar--) {
            }
            for (; nextFilteredChar > 0 && convertedKey[nextFilteredChar-1] < ALPHABET_SIZE; nextFilteredChar--) {
            }
        }

        if (effectiveKeyLength - i * advance - subPatternLength >= nextFilteredChar) {

            numOfSaIndexGroup = BWTHammingDistMatchOld(convertedKey + effectiveKeyLength - i*advance - subPatternLength, subPatternLength,
                                                  bwt, maxError, 
                                                  saIndexGroup + totalnumOfSaIndexGroup, maxnumOfSaIndexGroup - totalnumOfSaIndexGroup, 
                                                  effectiveKeyLength - i * advance - subPatternLength,
                                                  effectiveKeyLength - i * advance - subPatternLength);
            if (numOfSaIndexGroup >= 0) {
                totalnumOfSaIndexGroup += numOfSaIndexGroup;
            } else {
                // not enough memory to hold the SA index group for the sub pattern
                break;
            }
        }

    }
    *lengthInCurrentRound = i * advance;

    return totalnumOfSaIndexGroup;

}


// if discardDiagonalHit is set to TRUE, saIndexGroup[].info must = saIndexGroup[].posQuery;
// otherwise the wrong hit may be discarded for saIndexGroups that come from duplicate sub-patterns.

// saIndexGroup[].numOfMatch == 0 is a special flag to indicate that saIndexGroup[].startSaIndex is actually a text position
// ** This arrangement has not been implemented yet **

int BWTDPHit(const BWT *bwt, SaIndexGroupNew* __restrict saIndexGroup, const long long numOfSaIndexGroup, 
             const int firstSaIndexGroupToProcess, int* __restrict saIndexGroupProcessed,
             const int discardDiagonalHit,
             char* workingMemory, const int workingMemorySize,
             BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics) {

    #define MIN_LINKED_FACTOR    2

    HitList* __restrict hitList;

    int i, j, k;

    int* __restrict linkedHitList;

    int numOfHit, numOfLinkedHit;
    int numOfHitForPosQuery, numOfSaIndexForPosQuery;

    int lastDupSaIndexGroupIndex;
    int lastLinkedHitIndex;
    int diagonalHitIndex;
    int maxLinkedHit;

    int workingMemoryUsed;

    int maxnumOfHitPerSaIndexGroup, m, n;
    int saIndexGroupLimit;

    unsigned long long posQuery;
    int undecodedTextPositionIndex;
    unsigned long long undecodedSaIndex;

    int gapFromLastSubPattern[MAX_DIAGONAL_LEVEL+1];
    int numOfUndecodedSa[MAX_DIAGONAL_LEVEL+1];
    int numOfSaProcessed[MAX_DIAGONAL_LEVEL+1];
    int numOfSaIndex[MAX_DIAGONAL_LEVEL+1];
    SaIndexList* __restrict saIndexList[MAX_DIAGONAL_LEVEL + 2];
    int saIndexListIndex[MAX_DIAGONAL_LEVEL+1];    // saIndexListIndex[0] -> current sub pattern, saIndexListIndex[1] -> last sub pattern
    int tempSaIndexListIndex, currentSaIndexListIndex, gapToNextSubPattern;
    int index;

    int diagonalLevel;

    if (firstSaIndexGroupToProcess == 0) {

        // First sort saIndexGroup on startSaIndex + numOfMatch desc + info desc

        QSort(saIndexGroup, numOfSaIndexGroup, sizeof(SaIndexGroupNew), SaIndexGroupDPHitOrder1);

        // Second change posQuery for duplicated or enclosed SA index ranges

        i = 0;
        while (i < numOfSaIndexGroup) {
            j = i + 1;
            while (j < numOfSaIndexGroup && saIndexGroup[j].startSaIndex < saIndexGroup[i].startSaIndex + saIndexGroup[i].numOfMatch) {
                if (saIndexGroup[j].startSaIndex < saIndexGroup[i].startSaIndex) {
                    fprintf(stderr, "BWTDPHit(): Error sorting SA index group!\n");
                    exit(1);
                }
                saIndexGroup[j].posQuery = saIndexGroup[i].posQuery;
                j++;
            }
            i = j;
        }

        // Third sort saIndexGroup on posQuery desc + startSaIndex + numOfMatch desc + info desc
        QSort(saIndexGroup, numOfSaIndexGroup, sizeof(SaIndexGroupNew), SaIndexGroupDPHitOrder2);

    }

    // Fourth sum up the no. of hits and find the maximum no. of hits with the same posQuery

    workingMemoryUsed = 0;
    numOfHit = 0;
    maxnumOfHitPerSaIndexGroup = 0;
    saIndexGroupLimit = firstSaIndexGroupToProcess;
    diagonalLevel = MAX_DIAGONAL_LEVEL;

    i = firstSaIndexGroupToProcess;
    m = 0;
    while (i < numOfSaIndexGroup) {
        numOfHitForPosQuery = 0;
        numOfSaIndexForPosQuery = 0;
        j = i;
        while (j < numOfSaIndexGroup && saIndexGroup[j].posQuery == saIndexGroup[i].posQuery) {
            numOfHitForPosQuery += saIndexGroup[j].numOfMatch;
            numOfSaIndexForPosQuery += saIndexGroup[j].numOfMatch;
            k = j + 1;
            while (k < numOfSaIndexGroup && saIndexGroup[k].posQuery == saIndexGroup[i].posQuery
                                         && saIndexGroup[k].startSaIndex < saIndexGroup[j].startSaIndex + saIndexGroup[j].numOfMatch) {
                numOfHitForPosQuery += saIndexGroup[k].numOfMatch;
                k++;
            }
            j = k;
        }
        i = j;
        if (numOfSaIndexForPosQuery > m) {
            workingMemoryUsed += (numOfSaIndexForPosQuery - m) * (MAX_DIAGONAL_LEVEL+2) * sizeof(SaIndexList);
            m = numOfSaIndexForPosQuery;
        }
        workingMemoryUsed += numOfHitForPosQuery * sizeof(HitList) + (numOfSaIndexForPosQuery / MIN_LINKED_FACTOR) * sizeof(int);
        if (workingMemoryUsed <= workingMemorySize) {
            saIndexGroupLimit = i;
            maxnumOfHitPerSaIndexGroup = m;
            numOfHit += numOfHitForPosQuery;
        } else {
            // Try to decrease diagonal level
            while (workingMemoryUsed > workingMemorySize && diagonalLevel > 0) {
                workingMemoryUsed -= m * sizeof(SaIndexList);
                diagonalLevel--;
            }
            if (workingMemoryUsed <= workingMemorySize) {
                saIndexGroupLimit = i;
                maxnumOfHitPerSaIndexGroup = m;
                numOfHit += numOfHitForPosQuery;
            } else {
                diagonalLevel = MAX_DIAGONAL_LEVEL;
                break;
            }
        }
    }

    // Allocate memory

    if (saIndexGroupLimit == firstSaIndexGroupToProcess) {
        fprintf(stderr, "BWTDPHit(): Not enough memory to decode hits.\n");
        exit(1);
    }

    workingMemoryUsed = 0;
    hitList = (HitList*)(workingMemory);
    workingMemoryUsed += numOfHit * sizeof(HitList);

    for (i=0; i<diagonalLevel+2; i++) {
        saIndexList[i] = (SaIndexList*)(workingMemory + workingMemoryUsed);
        workingMemoryUsed += maxnumOfHitPerSaIndexGroup * sizeof(SaIndexList);
    }

    for (i=0; i<diagonalLevel+1; i++) {
        numOfSaIndex[i] = 0;
        gapFromLastSubPattern[i] = 0;
        saIndexListIndex[i] = i;
    }
    tempSaIndexListIndex = diagonalLevel + 1;

    linkedHitList = (int*)(workingMemory + workingMemoryUsed);
    maxLinkedHit = (workingMemorySize - workingMemoryUsed) / sizeof(int);

    // Start decoding

    numOfHit = 0;
    numOfLinkedHit = 0;
    lastDupSaIndexGroupIndex = -1;
    lastLinkedHitIndex = -1;

    j = firstSaIndexGroupToProcess;
    while (j < saIndexGroupLimit) {

        currentSaIndexListIndex = saIndexListIndex[0];
        numOfSaIndex[0] = 0;
        gapFromLastSubPattern[0] = 0;
        for (k=0; k<diagonalLevel+1; k++) {
            numOfUndecodedSa[k] = 0;
            numOfSaProcessed[k] = 0;
        }

        // determine the gap to the next sub pattern
        posQuery = saIndexGroup[j].posQuery;
        k = j + 1;
        while (k < saIndexGroupLimit && saIndexGroup[k].posQuery == posQuery) {
            k++;
        }
        if (k < saIndexGroupLimit) {
            gapToNextSubPattern = posQuery - saIndexGroup[k].posQuery;
        } else {
            gapToNextSubPattern = ALL_ONE_MASK;
        }

        while (j < saIndexGroupLimit && saIndexGroup[j].posQuery == posQuery) {

            // store in temp variable; the processed SA index group will reuse the space to store other values

            saIndexGroup[j].posQuery = numOfHit;    // place starting text position index in posQuery
    
            for (k=0; k<(int)saIndexGroup[j].numOfMatch; k++) {
                saIndexList[currentSaIndexListIndex][numOfSaIndex[0] + k].saIndex = saIndexGroup[j].startSaIndex + k;
                saIndexList[currentSaIndexListIndex][numOfSaIndex[0] + k].textPositionIndex = numOfHit + k;
                hitList[numOfHit + k].info = posQuery;
            }
            numOfSaIndex[0] += saIndexGroup[j].numOfMatch;
            numOfHit += saIndexGroup[j].numOfMatch;
            
            // Process the undecoded SA index in the last sub-pattern
            for (k=1; k<diagonalLevel+1; k++) {
                index = saIndexListIndex[k];
                while (numOfSaProcessed[k] < numOfSaIndex[k] &&
                       saIndexList[index][numOfSaProcessed[k]].saIndex < saIndexGroup[j].startSaIndex + saIndexGroup[j].numOfMatch) {
                    
                    undecodedTextPositionIndex = saIndexList[index][numOfSaProcessed[k]].textPositionIndex;
                    undecodedSaIndex = saIndexList[index][numOfSaProcessed[k]].saIndex;

                    if (undecodedSaIndex >= saIndexGroup[j].startSaIndex && numOfLinkedHit < maxLinkedHit) {
                        // link the SA
                        hitList[undecodedTextPositionIndex].posText = 
                            saIndexGroup[j].posQuery + undecodedSaIndex - saIndexGroup[j].startSaIndex;

                        // add to linked hit list
                        linkedHitList[numOfLinkedHit] = undecodedTextPositionIndex;
                        numOfLinkedHit++;

                        bwtSaRetrievalStatistics->saDiagonalLinked++;
                    } else {
                        // Pack the undecoded SA to the front
                        if (numOfSaProcessed[k] > numOfUndecodedSa[k]) {
                            saIndexList[index][numOfUndecodedSa[k]].saIndex = saIndexList[index][numOfSaProcessed[k]].saIndex;
                            saIndexList[index][numOfUndecodedSa[k]].textPositionIndex = saIndexList[index][numOfSaProcessed[k]].textPositionIndex;
                        }
                        numOfUndecodedSa[k]++;
                    }
                    numOfSaProcessed[k]++;
                }
                    
            }

            // Skip all duplicated or enclosed saIndexGroup
            k = j;
            j++;
            while (j < saIndexGroupLimit && saIndexGroup[j].posQuery == posQuery
                                         && saIndexGroup[j].startSaIndex < saIndexGroup[k].startSaIndex + saIndexGroup[k].numOfMatch) {
                saIndexGroup[j].posQuery = numOfHit;
                numOfHit += saIndexGroup[j].numOfMatch;
                j++;
            }

        }

        for (k=0; k<diagonalLevel+1; k++) {

            index = saIndexListIndex[k];
                
            // Process the remaining undecoded SA index in the last sub-pattern
            if (numOfSaProcessed[k] > numOfUndecodedSa[k]) {
                while (numOfSaProcessed[k] < numOfSaIndex[k]) {
                    saIndexList[index][numOfUndecodedSa[k]].saIndex = saIndexList[index][numOfSaProcessed[k]].saIndex;
                    saIndexList[index][numOfUndecodedSa[k]].textPositionIndex = saIndexList[index][numOfSaProcessed[k]].textPositionIndex;
                    numOfUndecodedSa[k]++;
                    numOfSaProcessed[k]++;
                }
            } else {
                // No SA was linked
                numOfUndecodedSa[k] = numOfSaIndex[k];
            }

            if (numOfUndecodedSa[k] > 0) {

                if (index != saIndexListIndex[diagonalLevel]) {
                    // Decode to the next sub-pattern
                    numOfUndecodedSa[k] = BWTDecodeTextPosition(bwt, saIndexList[index], saIndexList[tempSaIndexListIndex],
                                                        numOfUndecodedSa[k], gapFromLastSubPattern[k], gapToNextSubPattern,
                                                         (HitList*)hitList, bwtSaRetrievalStatistics);
                    if (gapToNextSubPattern % 2 == 1) {
                        // Swap the index if gap is odd
                        saIndexListIndex[k] = tempSaIndexListIndex;
                        tempSaIndexListIndex = index;
                    }
                } else {
                    // Fully decode the SA index
                    BWTDecodeTextPosition(bwt, saIndexList[index], saIndexList[tempSaIndexListIndex],
                                          numOfUndecodedSa[k], gapFromLastSubPattern[k], ALL_ONE_MASK_64, (HitList*)hitList,
                                          bwtSaRetrievalStatistics);
                }

            }

        }

        currentSaIndexListIndex = saIndexListIndex[diagonalLevel];
        for (k=diagonalLevel; k>0; k--) {
            saIndexListIndex[k] = saIndexListIndex[k-1];
            gapFromLastSubPattern[k] = gapFromLastSubPattern[k-1] + gapToNextSubPattern;
            numOfSaIndex[k] = numOfUndecodedSa[k-1];
        }
        saIndexListIndex[0] = currentSaIndexListIndex;

    }

    // Decode all SA index
    for (k=1; k<diagonalLevel+1; k++) {
        index = saIndexListIndex[k];
        if (numOfSaIndex[index] > 0) {
            // Fully decode the SA index
            BWTDecodeTextPosition(bwt, saIndexList[index], saIndexList[tempSaIndexListIndex],
                                      numOfSaIndex[k], gapFromLastSubPattern[k], ALL_ONE_MASK_64, (HitList*)hitList,
                                      bwtSaRetrievalStatistics);
        }
    }

    // resolve linked hit
    for (i=numOfLinkedHit; i--;) {
        diagonalHitIndex = hitList[linkedHitList[i]].posText;
        posQuery = hitList[diagonalHitIndex].info & ALL_BUT_FIRST_BIT_MASK;    // remove the first bit
        hitList[linkedHitList[i]].posText = hitList[linkedHitList[i]].info - posQuery + hitList[diagonalHitIndex].posText;
        // discard the hit by setting the first bit of posQuery
        hitList[diagonalHitIndex].info |= (discardDiagonalHit << BITS_IN_WORD_MINUS_1);
    }

    n = 0;

    i = firstSaIndexGroupToProcess;
    while (i<saIndexGroupLimit) {
        for (j=0; j<(int)saIndexGroup[i].numOfMatch; j++) {
            if (!discardDiagonalHit || (hitList[saIndexGroup[i].posQuery + j].info & FIRST_BIT_MASK) == 0) {
                hitList[n].posText = hitList[saIndexGroup[i].posQuery + j].posText;
                hitList[n].info = saIndexGroup[i].info;
                n++;
            }
        }
        j = i + 1;
        while (j < saIndexGroupLimit && saIndexGroup[j].startSaIndex >= saIndexGroup[i].startSaIndex
                                     && saIndexGroup[j].startSaIndex < saIndexGroup[i].startSaIndex + saIndexGroup[i].numOfMatch) {
            for (k=0; k<(int)saIndexGroup[j].numOfMatch; k++) {
                hitList[n].posText = hitList[saIndexGroup[i].posQuery + saIndexGroup[j].startSaIndex + k - saIndexGroup[i].startSaIndex].posText;
                hitList[n].info = saIndexGroup[j].info;
                n++;
            }
            j++;
        }
        i = j;
    }

    numOfHit = n;
    *saIndexGroupProcessed = saIndexGroupLimit - firstSaIndexGroupToProcess;

    return numOfHit;

}

/*
unsigned int BWTSubPatternHammingDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
                                  const unsigned int maxError, const unsigned int skip, const unsigned int matchBitVector,
                                  const SaIndexRange *cachedSaIndex, const unsigned int cachedSaIndexNumOfChar,
                                  char *workingMemory, const unsigned int workingMemorySize, const unsigned int sortOption,
                                  BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics) {

    unsigned int i, j, k;
    unsigned int numOfSubPattern, advance;

    unsigned int numOfSaGroup, totalnumOfSaGroup;
    unsigned int numOfHit, totalnumOfHit;
    unsigned int maxnumOfHitPerSaIndexGroup;
    unsigned int nextFilteredChar;

    // Storage to be allocated with workingMemory
    unsigned int* __restrict firstSaIndexGroupForSubPattern;
    SaIndexGroupNew* __restrict saIndexGroup;
    HitListWithPosQuery* __restrict hitList;
    unsigned int* __restrict linkedHitList;
    SaIndexGroupHash* __restrict saIndexGroupHash;

    HitListWithPosQuery* __restrict finalHitList;
    unsigned int bucketCount[NUM_BUCKET + 1];
    unsigned int firstSaIndexGroupMemoryUsed;
    unsigned int saIndexGroupMemoryUsed;
    unsigned int hitListMemoryUsed;
    unsigned int saIndexListMemoryUsed;
    unsigned int linkedHitListMemoryUsed;
    unsigned int saIndexGroupHashMemoryUsed;
    unsigned int maxnumOfSaIndexGroup;

    unsigned int diagonalLevel;
    unsigned int maxnumOfHit, maxnumOfLinkedHit;
    unsigned int numOfLinkedHit;
    unsigned int nextSubPatternIndex;
    unsigned int lastDupSaIndexGroupIndex;
    unsigned int lastLinkedHitIndex;
    unsigned int diagonalHitIndex;
    unsigned int posQuery, error, startSaIndex, numOfMatch;
    unsigned int undecodedTextPositionIndex, undecodedSaIndex;
    unsigned int hashedSaIndexGroupIndex;
    unsigned int originalSaIndexGroupIndex, dupTextPositionIndex, originalTextPositionIndex;

    unsigned int gapFromLastSubPattern[MAX_DIAGONAL_LEVEL+1];
    unsigned int numOfUndecodedSa[MAX_DIAGONAL_LEVEL+1];
    unsigned int numOfSaProcessed[MAX_DIAGONAL_LEVEL+1];
    unsigned int numOfSaIndex[MAX_DIAGONAL_LEVEL+1];
    SaIndexList* __restrict saIndexList[MAX_DIAGONAL_LEVEL + 2];
    unsigned int saIndexListIndex[MAX_DIAGONAL_LEVEL+1];    // saIndexListIndex[0] -> current sub pattern, saIndexListIndex[1] -> last sub pattern
    unsigned int tempSaIndexListIndex, currentSaIndexListIndex, gapToNextSubPattern;
    unsigned int index;

    // Hash variables
    float hashLoadFactor = 0.5;
    unsigned int hashBuffer = 100;
    unsigned int hashCoefficient1, hashCoefficient2, hashCoefficient3;
    unsigned int hashValue;

    // Determine the number of sub-patterns
    advance = skip + 1;
    numOfSubPattern = (keyLength - subPatternLength + advance) / advance;

    // Allocate memory for temporary variables
    firstSaIndexGroupForSubPattern = (unsigned int*)workingMemory;
    firstSaIndexGroupMemoryUsed = (numOfSubPattern + 1) * sizeof(unsigned int);
    saIndexGroup = (SaIndexGroupNew*)(workingMemory + firstSaIndexGroupMemoryUsed);
    maxnumOfSaIndexGroup = (workingMemorySize - firstSaIndexGroupMemoryUsed) / sizeof(SaIndexGroupNew);

    firstSaIndexGroupForSubPattern[0] = 0;
    maxnumOfHitPerSaIndexGroup = 0;
    
    // Find SA index groups for all sub-patterns

    for (nextFilteredChar = keyLength; nextFilteredChar > 0 && convertedKey[nextFilteredChar-1] < ALPHABET_SIZE; nextFilteredChar--) {
    }

    totalnumOfHit = 0;

    for (i=0; i<numOfSubPattern; i++) {

        if (nextFilteredChar > keyLength - i*advance) {
            // Advance to the next 'N'
            for (; nextFilteredChar > keyLength - i*advance && convertedKey[nextFilteredChar-1] >= ALPHABET_SIZE; nextFilteredChar--) {
            }
            for (; nextFilteredChar > 0 && convertedKey[nextFilteredChar-1] < ALPHABET_SIZE; nextFilteredChar--) {
            }
        }

        if (keyLength - i * advance - subPatternLength >= nextFilteredChar) {

            numOfSaGroup = BWTHammingDistMatch(convertedKey + keyLength - i*advance - subPatternLength, subPatternLength,
                                                  bwt, maxError, matchBitVector,
                                                  cachedSaIndex, cachedSaIndexNumOfChar,
                                                  saIndexGroup + firstSaIndexGroupForSubPattern[i], maxnumOfSaIndexGroup - firstSaIndexGroupForSubPattern[i]);
            firstSaIndexGroupForSubPattern[i+1] = firstSaIndexGroupForSubPattern[i] + numOfSaGroup;

            // determine the number of low occ SA and total number of hits
            numOfHit = 0;
            for (j=0; j<numOfSaGroup; j++) {
                numOfHit += saIndexGroup[firstSaIndexGroupForSubPattern[i]+j].numOfMatch;
            }
            if (numOfHit > maxnumOfHitPerSaIndexGroup) {
                maxnumOfHitPerSaIndexGroup = numOfHit;
            }
            totalnumOfHit += numOfHit;

        } else {
            // sub-pattern filtered
            firstSaIndexGroupForSubPattern[i+1] = firstSaIndexGroupForSubPattern[i];
        }

    }
    totalnumOfSaGroup = firstSaIndexGroupForSubPattern[numOfSubPattern];

    if (totalnumOfSaGroup == 0) {
        return 0;
    }

    // set up hash variables
    hashCoefficient1 = r250();
    hashCoefficient2 = r250();
    hashCoefficient3 = nextPrime((unsigned int)((float)totalnumOfSaGroup / hashLoadFactor));
    saIndexGroupHashMemoryUsed = (hashCoefficient3 + hashBuffer) * sizeof(SaIndexGroupHash);

    // Determine whether there is enough memory to store the hits
    // diagonalLevel must be at least 0, i.e. there must be enough memory for 2 saIndexList to be allocated

    saIndexGroupMemoryUsed = totalnumOfSaGroup * sizeof(SaIndexGroupNew);

    if (workingMemorySize < firstSaIndexGroupMemoryUsed + saIndexGroupMemoryUsed + saIndexGroupHashMemoryUsed
                            + sizeof(SaIndexList) * maxnumOfHitPerSaIndexGroup * 2) {
        maxnumOfHit = 0;
    } else {
        maxnumOfHit = (workingMemorySize - firstSaIndexGroupMemoryUsed - saIndexGroupMemoryUsed - saIndexGroupHashMemoryUsed
                           - sizeof(SaIndexList) * maxnumOfHitPerSaIndexGroup * 2) / sizeof(HitList);
    }
    if ((firstSaIndexGroupMemoryUsed + saIndexGroupMemoryUsed) * 2 > workingMemorySize || totalnumOfHit > maxnumOfHit) {
        fprintf(stderr, "Not enough memory allocated for hit generation!\n");
        return 0;
    }

    hitListMemoryUsed = sizeof(HitListWithPosQuery) * totalnumOfHit;

    if (sortOption == SORT_16_BIT && workingMemorySize > firstSaIndexGroupMemoryUsed + saIndexGroupMemoryUsed + saIndexGroupHashMemoryUsed
                                    + sizeof(SaIndexList) * maxnumOfHitPerSaIndexGroup * 2 + hitListMemoryUsed * 2) {
        // enough memory for bucket sort
        // allocate memory at end of working memory for hits
        hitList = (HitListWithPosQuery*)(workingMemory + workingMemorySize - hitListMemoryUsed);

        // allocate memory for saIndexGroupHash
        saIndexGroupHash = (SaIndexGroupHash*)(workingMemory + firstSaIndexGroupMemoryUsed + saIndexGroupMemoryUsed);
        for (i=0; i<hashCoefficient3 + hashBuffer; i++) {
            saIndexGroupHash[i].startSaIndex = 0;
        }

        // allocate memory for saIndexList
        diagonalLevel = (workingMemorySize - firstSaIndexGroupMemoryUsed - saIndexGroupMemoryUsed - saIndexGroupHashMemoryUsed - hitListMemoryUsed)
                        / (sizeof(SaIndexList) * maxnumOfHitPerSaIndexGroup) - 2;
        if (diagonalLevel > MAX_DIAGONAL_LEVEL) {
            diagonalLevel = MAX_DIAGONAL_LEVEL;
        }

        saIndexListMemoryUsed = 0;
        for (i=0; i<diagonalLevel+2; i++) {
            saIndexList[i] = (SaIndexList*)(workingMemory + firstSaIndexGroupMemoryUsed + saIndexGroupMemoryUsed + saIndexGroupHashMemoryUsed + saIndexListMemoryUsed);
            saIndexListMemoryUsed += sizeof(SaIndexList) * maxnumOfHitPerSaIndexGroup;
        }

        // Allocate the remaining memory to linkedHitList
        linkedHitListMemoryUsed = workingMemorySize - firstSaIndexGroupMemoryUsed - saIndexGroupMemoryUsed
                                  - saIndexGroupHashMemoryUsed - hitListMemoryUsed - saIndexListMemoryUsed;
        maxnumOfLinkedHit = linkedHitListMemoryUsed / sizeof(unsigned int);
        linkedHitList = (unsigned int*)(workingMemory + firstSaIndexGroupMemoryUsed + saIndexGroupMemoryUsed + saIndexGroupHashMemoryUsed + saIndexListMemoryUsed);

    } else {

        if (sortOption == SORT_16_BIT) {
            // not enough memory for bucket sort
            fprintf(stderr, "Not enough memory for faster sorting!\n");
        }

        // move firstSaIndexGroupForSubPattern and saIndexGroup to end of working memory
        memcpy(workingMemory + workingMemorySize - firstSaIndexGroupMemoryUsed, firstSaIndexGroupForSubPattern, firstSaIndexGroupMemoryUsed);
        firstSaIndexGroupForSubPattern = (unsigned int*)(workingMemory + workingMemorySize - firstSaIndexGroupMemoryUsed);
        memcpy(workingMemory - firstSaIndexGroupMemoryUsed - saIndexGroupMemoryUsed, saIndexGroup, saIndexGroupMemoryUsed);
        saIndexGroup = (SaIndexGroupNew*)(workingMemory - firstSaIndexGroupMemoryUsed - saIndexGroupMemoryUsed);

        // allocate memory for saIndexGroupHash
        saIndexGroupHash = (SaIndexGroupHash*)(workingMemory - firstSaIndexGroupMemoryUsed - saIndexGroupMemoryUsed - saIndexGroupHashMemoryUsed);
        for (i=0; i<hashCoefficient3 + hashBuffer; i++) {
            saIndexGroupHash[i].startSaIndex = 0;
        }

        // allocate memory at start of working memory for hits
        hitList = (HitListWithPosQuery*)workingMemory;

        // allocate memory for saIndexList
        diagonalLevel = (workingMemorySize - firstSaIndexGroupMemoryUsed - saIndexGroupMemoryUsed - saIndexGroupHashMemoryUsed - hitListMemoryUsed)
                        / (sizeof(SaIndexList) * maxnumOfHitPerSaIndexGroup) - 2;
        if (diagonalLevel > MAX_DIAGONAL_LEVEL) {
            diagonalLevel = MAX_DIAGONAL_LEVEL;
        }

        saIndexListMemoryUsed = 0;
        for (i=0; i<diagonalLevel+2; i++) {
            saIndexList[i] = (SaIndexList*)(workingMemory + hitListMemoryUsed + saIndexListMemoryUsed);
            saIndexListMemoryUsed += sizeof(SaIndexList) * maxnumOfHitPerSaIndexGroup;
        }

        // Allocate the remaining memory to linkedHitList
        linkedHitListMemoryUsed = workingMemorySize - firstSaIndexGroupMemoryUsed - saIndexGroupMemoryUsed
                                  - saIndexGroupHashMemoryUsed - hitListMemoryUsed - saIndexListMemoryUsed;
        maxnumOfLinkedHit = linkedHitListMemoryUsed / sizeof(unsigned int);
        linkedHitList = (unsigned int*)(workingMemory + hitListMemoryUsed + saIndexListMemoryUsed);

    }


    for (i=0; i<diagonalLevel+1; i++) {
        numOfSaIndex[i] = 0;
        gapFromLastSubPattern[i] = 0;
        saIndexListIndex[i] = i;
    }
    tempSaIndexListIndex = diagonalLevel + 1;

    totalnumOfHit = 0;
    numOfLinkedHit = 0;
    lastDupSaIndexGroupIndex = ALL_ONE_MASK;
    lastLinkedHitIndex = ALL_ONE_MASK;

    i = 0;
    while (i<numOfSubPattern) {

        currentSaIndexListIndex = saIndexListIndex[0];
        numOfSaIndex[0] = 0;
        gapFromLastSubPattern[0] = 0;
        for (k=0; k<diagonalLevel+1; k++) {
            numOfUndecodedSa[k] = 0;
            numOfSaProcessed[k] = 0;
        }

        posQuery = keyLength - i*advance - subPatternLength;

        // determine the gap to the next sub pattern
        nextSubPatternIndex = i + 1;
        while (nextSubPatternIndex < numOfSubPattern && firstSaIndexGroupForSubPattern[nextSubPatternIndex] == firstSaIndexGroupForSubPattern[nextSubPatternIndex+1]) {
            nextSubPatternIndex++;
        }
        if (nextSubPatternIndex < numOfSubPattern) {
            gapToNextSubPattern = (nextSubPatternIndex - i) * advance;
        } else {
            // there is no next sub pattern
            gapToNextSubPattern = ALL_ONE_MASK;
        }

        for (j=firstSaIndexGroupForSubPattern[i]; j<firstSaIndexGroupForSubPattern[i+1]; j++) {

            // store in temp variable; the processed SA index group will reuse the space to store other values
            error = saIndexGroup[j].info;
            startSaIndex = saIndexGroup[j].startSaIndex;
            numOfMatch = saIndexGroup[j].numOfMatch;

            ((SaIndexGroupProcessed*)saIndexGroup + j)->textPositionIndex = totalnumOfHit;

            // first hash to see if it's duplicate
            hashValue = (startSaIndex * hashCoefficient1 + (startSaIndex >> 16) * hashCoefficient2) % hashCoefficient3;
            while (saIndexGroupHash[hashValue].startSaIndex != 0 && saIndexGroupHash[hashValue].startSaIndex != startSaIndex) {
                hashValue++;
                if (hashValue >= hashCoefficient3 + hashBuffer) {
                    fprintf(stderr, "BWTSubPatternHammingDistMatch(): hashValue >= maxHashTableSize!\n");
                    exit(1);
                }
            }
            if (saIndexGroupHash[hashValue].startSaIndex == 0) {
                saIndexGroupHash[hashValue].startSaIndex = startSaIndex;
                saIndexGroupHash[hashValue].saIndexGroupIndex = j;
            } else {

                // duplicate
                hashedSaIndexGroupIndex = saIndexGroupHash[hashValue].saIndexGroupIndex;

                // Store posQuery in hit list
                for (k=0; k<numOfMatch; k++) {
                    hitList[totalnumOfHit + k].posQuery = posQuery;
                }
                ((DupSaIndexGroup*)saIndexGroup + j)->saIndexGroupIndex = saIndexGroupHash[hashValue].saIndexGroupIndex;
                ((DupSaIndexGroup*)saIndexGroup + j)->textPositionIndex = totalnumOfHit;

                // A linked list for duplicates SA index group
                ((DupSaIndexGroup*)saIndexGroup + j)->lastDupSaIndexGroupIndex = lastDupSaIndexGroupIndex;
                lastDupSaIndexGroupIndex = j;

                bwtSaRetrievalStatistics->saDuplicated += numOfMatch;

                totalnumOfHit += numOfMatch;

                continue;
            }
    
            for (k=0; k<numOfMatch; k++) {
                saIndexList[currentSaIndexListIndex][numOfSaIndex[0] + k].saIndex = startSaIndex + k;
                saIndexList[currentSaIndexListIndex][numOfSaIndex[0] + k].textPositionIndex = totalnumOfHit + k;
                // Store posQuery in hit list
                hitList[totalnumOfHit + k].posQuery = posQuery;
            }
            numOfSaIndex[0] += numOfMatch;
            totalnumOfHit += numOfMatch;
            
            // Process the undecoded SA index in the last sub-pattern
            for (k=1; k<diagonalLevel+1; k++) {
                index = saIndexListIndex[k];
                while (numOfSaProcessed[k] < numOfSaIndex[k] &&
                       saIndexList[index][numOfSaProcessed[k]].saIndex < startSaIndex + numOfMatch) {
                    
                    undecodedTextPositionIndex = saIndexList[index][numOfSaProcessed[k]].textPositionIndex;
                    undecodedSaIndex = saIndexList[index][numOfSaProcessed[k]].saIndex;

                    if (undecodedSaIndex >= startSaIndex && numOfLinkedHit < maxnumOfLinkedHit) {
                        // link the SA
                        hitList[undecodedTextPositionIndex].posText = 
                            ((SaIndexGroupProcessed*)saIndexGroup + j)->textPositionIndex + undecodedSaIndex - startSaIndex;

                        // add to linked hit list
                        linkedHitList[numOfLinkedHit] = undecodedTextPositionIndex;
                        numOfLinkedHit++;

                        bwtSaRetrievalStatistics->saDiagonalLinked++;
                    } else {
                        // Pack the undecoded SA to the front
                        if (numOfSaProcessed[k] > numOfUndecodedSa[k]) {
                            saIndexList[index][numOfUndecodedSa[k]].saIndex = saIndexList[index][numOfSaProcessed[k]].saIndex;
                            saIndexList[index][numOfUndecodedSa[k]].textPositionIndex = saIndexList[index][numOfSaProcessed[k]].textPositionIndex;
                        }
                        numOfUndecodedSa[k]++;
                    }
                    numOfSaProcessed[k]++;
                }
                    
            }

        }

        for (k=0; k<diagonalLevel+1; k++) {

            index = saIndexListIndex[k];
                
            // Process the remaining undecoded SA index in the last sub-pattern
            if (numOfSaProcessed[k] > numOfUndecodedSa[k]) {
                while (numOfSaProcessed[k] < numOfSaIndex[k]) {
                    saIndexList[index][numOfUndecodedSa[k]].saIndex = saIndexList[index][numOfSaProcessed[k]].saIndex;
                    saIndexList[index][numOfUndecodedSa[k]].textPositionIndex = saIndexList[index][numOfSaProcessed[k]].textPositionIndex;
                    numOfUndecodedSa[k]++;
                    numOfSaProcessed[k]++;
                }
            } else {
                // No SA was linked
                numOfUndecodedSa[k] = numOfSaIndex[k];
            }

            if (numOfUndecodedSa[k] > 0) {

                if (index != saIndexListIndex[diagonalLevel]) {
                    // Decode to the next sub-pattern
                    numOfUndecodedSa[k] = BWTDecodeTextPosition(bwt, saIndexList[index], saIndexList[tempSaIndexListIndex],
                                                        numOfUndecodedSa[k], gapFromLastSubPattern[k], gapToNextSubPattern,
                                                         (HitList*)hitList, bwtSaRetrievalStatistics);
                    if (gapToNextSubPattern % 2 == 1) {
                        // Swap the index if gap is odd
                        saIndexListIndex[k] = tempSaIndexListIndex;
                        tempSaIndexListIndex = index;
                    }
                } else {
                    // Fully decode the SA index
                    BWTDecodeTextPosition(bwt, saIndexList[index], saIndexList[tempSaIndexListIndex],
                                          numOfUndecodedSa[k], gapFromLastSubPattern[k], ALL_ONE_MASK, (HitList*)hitList,
                                          bwtSaRetrievalStatistics);
                }

            }

        }

        i = nextSubPatternIndex;

        currentSaIndexListIndex = saIndexListIndex[diagonalLevel];
        for (k=diagonalLevel; k>0; k--) {
            saIndexListIndex[k] = saIndexListIndex[k-1];
            gapFromLastSubPattern[k] = gapFromLastSubPattern[k-1] + gapToNextSubPattern;
            numOfSaIndex[k] = numOfUndecodedSa[k-1];
        }
        saIndexListIndex[0] = currentSaIndexListIndex;

    }

    // Decode all SA index
    for (k=1; k<diagonalLevel+1; k++) {
        index = saIndexListIndex[k];
        if (numOfSaIndex[index] > 0) {
            // Fully decode the SA index
            BWTDecodeTextPosition(bwt, saIndexList[index], saIndexList[tempSaIndexListIndex],
                                      numOfSaIndex[k], gapFromLastSubPattern[k], ALL_ONE_MASK, (HitList*)hitList,
                                      bwtSaRetrievalStatistics);
        }
    }

    // resolve linked hit
    for (i=numOfLinkedHit; i--;) {
        diagonalHitIndex = hitList[linkedHitList[i]].posText;
        posQuery = hitList[diagonalHitIndex].posQuery & ALL_BUT_FIRST_BIT_MASK;    // remove the first bit
        hitList[linkedHitList[i]].posText = hitList[linkedHitList[i]].posQuery - posQuery + hitList[diagonalHitIndex].posText;
        // discard the hit by setting the first bit of posQuery
        hitList[diagonalHitIndex].posQuery |= FIRST_BIT_MASK;
    }

    // resolve duplicate SA index groups
    while (lastDupSaIndexGroupIndex != ALL_ONE_MASK) {
        
        originalSaIndexGroupIndex = ((DupSaIndexGroup*)saIndexGroup + lastDupSaIndexGroupIndex)->saIndexGroupIndex;
        dupTextPositionIndex = ((DupSaIndexGroup*)saIndexGroup + lastDupSaIndexGroupIndex)->textPositionIndex;
        originalTextPositionIndex = ((SaIndexGroupProcessed*)saIndexGroup + originalSaIndexGroupIndex)->textPositionIndex;
        numOfMatch = ((SaIndexGroupProcessed*)saIndexGroup + originalSaIndexGroupIndex)->numOfMatch;

        for (i=0; i<numOfMatch; i++) {
            hitList[dupTextPositionIndex + i].posText = hitList[originalTextPositionIndex + i].posText;
        }

        lastDupSaIndexGroupIndex = ((DupSaIndexGroup*)saIndexGroup + lastDupSaIndexGroupIndex)->lastDupSaIndexGroupIndex;

    }

    if (sortOption != SORT_NONE) {

        if (hitList != (HitListWithPosQuery*)workingMemory) {

            // bucket sort the first 16 bit of posText
            for (i=0; i<=NUM_BUCKET; i++) {
                bucketCount[i] = 0;
            }

            for (i=0; i<totalnumOfHit; i++) {
                if ((hitList[i].posQuery & FIRST_BIT_MASK) == 0) {
                    // the hit has not been discarded
                    bucketCount[(hitList[i].posText >> BUCKET_BIT) + 1]++;
                }
            }
            for (i=0; i<NUM_BUCKET; i++) {
                bucketCount[i+1] += bucketCount[i];
            }

            finalHitList = (HitListWithPosQuery*)workingMemory;

            for (i=0; i<totalnumOfHit; i++) {
                if ((hitList[i].posQuery & FIRST_BIT_MASK) == 0) {
                    j = hitList[i].posText >> BUCKET_BIT;
                    finalHitList[bucketCount[j]].posQuery = hitList[i].posQuery;
                    finalHitList[bucketCount[j]].posText = hitList[i].posText;
                    bucketCount[j]++;
                }
            }

            totalnumOfHit = bucketCount[NUM_BUCKET];

            if (sortOption == SORT_ALL) {
                if (bucketCount[0] > 0) {
                    QSort(hitList, bucketCount[0], sizeof(HitListWithPosQuery), HitListPosTextOrder);
                }
                for (i=1; i<NUM_BUCKET; i++) {
                    if (bucketCount[i] > bucketCount[i-1]) {
                        QSort(hitList + bucketCount[i-1], bucketCount[i] - bucketCount[i-1], sizeof(HitListWithPosQuery), HitListPosTextOrder);
                    }
                }
            }

        } else {

            // Compress hitList by removing discarded hits
            j = 0;
            for (i=0; i<totalnumOfHit; i++) {
                if ((hitList[i].posQuery & FIRST_BIT_MASK) == 0) {
                    hitList[j].posQuery = hitList[i].posQuery;
                    hitList[j].posText = hitList[i].posText;
                    j++;
                }
            }

            totalnumOfHit = j;

            if (sortOption != SORT_16_BIT) {
                QSort(hitList, totalnumOfHit, sizeof(HitListWithPosQuery), HitListPosText16BitOrder);
            } else {
                QSort(hitList, totalnumOfHit, sizeof(HitListWithPosQuery), HitListPosTextOrder);
            }
            
        }

    } else {
        // Compress hitList by removing discarded hits
        j = 0;
        for (i=0; i<totalnumOfHit; i++) {
            if ((hitList[i].posQuery & FIRST_BIT_MASK) == 0) {
                hitList[j].posQuery = hitList[i].posQuery;
                hitList[j].posText = hitList[i].posText;
                j++;
            }
        }
        totalnumOfHit = j;
    }

    return totalnumOfHit;

}
*/
unsigned long long BWTEditDistMaxSaIndexGroup(const unsigned int keyLength, const unsigned int maxError) {

    unsigned long long a, c, d, e, i;
    unsigned long long t;

    unsigned long long deleteHammingCombination[MAX_APPROX_MATCH_ERROR];

    t = 1;
    a = 1;
    c = 1;
    d = 1;

    deleteHammingCombination[0] = 1;

    for (e=1; e<=maxError; e++) {
        // For error = e, combination for delete + change = mCe (alphabetSize) ^ e, where m = keyLength
        c *= keyLength + 1 - e;
        d *= e;
        a *= ALPHABET_SIZE;
        deleteHammingCombination[e] = c / d * a;
    }

    for (e=1; e<=maxError; e++) {
        // For insertion, approximate by multiplying with (m * alphabetSize) ^ error
        c = 1;
        for (i=0; i<=e; i++) {
            t += deleteHammingCombination[e - i] * c;
            c *= keyLength * ALPHABET_SIZE;
        }

    }

    // The number of combination for insert is over-estimated and that hopefully is enough for spliting the delete groups
    return t;

}

unsigned long long BWTEditDistMatchOld(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
                     SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned long long maxSaIndexGroup) {

    // stack
    unsigned long long startSaIndex[MAX_ARPROX_MATCH_LENGTH + 1];
    unsigned long long endSaIndex[MAX_ARPROX_MATCH_LENGTH];
    unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH+1];
    unsigned long long errorPos[MAX_APPROX_MATCH_ERROR + 1];

    // for insertion
    unsigned long long startSaIndexBefore[MAX_APPROX_MATCH_ERROR + 1];
    unsigned long long endSaIndexBefore[MAX_APPROX_MATCH_ERROR];
    unsigned char charBefore[MAX_APPROX_MATCH_ERROR + 1];

    unsigned long long numOfSaGroup;
    unsigned long long numOfError;

    unsigned int c, e;
    unsigned int errorAdded;
    unsigned char maxPatternValue;
    unsigned long long tempSaIndexLeft, tempSaIndexRight;

    unsigned long long length;

    unsigned int exceedMaxSaIndexGroupPrinted = FALSE;


    #ifdef DEBUG
    if (maxError > keyLength) {
        fprintf(stderr, "BWTEditDistMatch() : maxError > keyLength!\n");
        exit(1);
    }
    if (keyLength > MAX_ARPROX_MATCH_LENGTH) {
        fprintf(stderr, "BWTEditDistMatch() : keyLength > MAX_ARPROX_MATCH_LENGTH!\n");
        exit(1);
    }
    if (maxError > MAX_APPROX_MATCH_ERROR) {
        fprintf(stderr, "BWTEditDistMatch() : maxError > MAX_APPROX_MATCH_ERROR!\n");
        exit(1);
    }
    #endif

    // Set this boundary case so that the last character of generated pattern can be located by backward search
    startSaIndex[keyLength] = 0;
    endSaIndex[keyLength] = bwt->textLength;
    maxPatternValue = (unsigned char)(ALPHABET_SIZE * 2);    // 0 to alphabetSize - 1 = edit; alphabetSize = delete; alphabetSize + 1 to alphabetsize * 2 = insert

    length = keyLength;

    // Exact match
    numOfSaGroup = BWTBackwardSearch(convertedKey, keyLength, bwt, &tempSaIndexLeft, &tempSaIndexRight);
    if (numOfSaGroup > 0) {
        saIndexGroup[0].startSaIndex = tempSaIndexLeft;
        saIndexGroup[0].numOfMatch = tempSaIndexRight - tempSaIndexLeft + 1;
        saIndexGroup[0].length = length;
        saIndexGroup[0].error = 0;
    }

    // With errors
    for (numOfError = 1; numOfError <= maxError; numOfError++) {

        memcpy(generatedPattern, convertedKey, keyLength);
        length = keyLength;

        // Set initial state
        e = errorPos[1] = keyLength - 1;
        generatedPattern[e] = (unsigned char)-1;
        errorAdded = 1;

        // store before change data
        c = charBefore[1] = convertedKey[e];
        startSaIndexBefore[1] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
        endSaIndexBefore[1] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);

        while (errorAdded > 0) {
            
            // Increment the last error
            e = errorPos[errorAdded];
            generatedPattern[e]++;
            if (generatedPattern[e] == convertedKey[e]) {
                generatedPattern[e]++;
            }
            if (generatedPattern[e] > maxPatternValue) {

                length--;

                // Cannot increment the last error; try moving the error to the left
                if (e > 1 || (e == 1 && numOfError - errorAdded == 0)) {    // only the last error can move to the leftmost position as insertion is not done on pattern boundary

                    // restore before change data
                    c = generatedPattern[e] = charBefore[errorAdded];
                    startSaIndex[e] = startSaIndexBefore[errorAdded];
                    endSaIndex[e] = endSaIndexBefore[errorAdded];

                    // move error
                    errorPos[errorAdded]--;
                    e = errorPos[errorAdded];
                    generatedPattern[e] = (convertedKey[e] == 0);    // 1 if convertedKey = 0, 0 otherwise

                    // store before change data
                    c = charBefore[errorAdded] = convertedKey[e];
                    startSaIndexBefore[errorAdded] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                    endSaIndexBefore[errorAdded] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);

                } else {
                    // Cannot move error to the left
                    // restore before change data
                    c = generatedPattern[e] = charBefore[errorAdded];
                    startSaIndex[e] = startSaIndexBefore[errorAdded];
                    endSaIndex[e] = endSaIndexBefore[errorAdded];

                    errorAdded--;
                    continue;
                }
            }

            if (generatedPattern[e] < ALPHABET_SIZE) {
                // edit
                c = generatedPattern[e];
                 startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
            } else {
                if (generatedPattern[e] > ALPHABET_SIZE) {
                    if (e > 0) {
                        // insert
                        c = generatedPattern[e] - ALPHABET_SIZE - 1;    // get back the character to be inserted
                         startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndexBefore[errorAdded], c) + 1;
                        endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndexBefore[errorAdded] + 1, c);
                        if (c == 0) {
                            // insertion following deletion
                            length += 2;
                        }
                    } else {
                        // do not insert on pattern boundary
                        generatedPattern[e] = maxPatternValue;
                        length += 2;    // it will be subtracted when error is shifted
                        continue;
                    }
                } else {
                    // delete
                    startSaIndex[e] = startSaIndex[e+1];
                    endSaIndex[e] = endSaIndex[e+1];
                    length--;
                }
            }

            while (errorAdded < numOfError && startSaIndex[e] <= endSaIndex[e]) {

                // Add insertions
                errorAdded++;
                errorPos[errorAdded] = errorPos[errorAdded-1];

                // store before change data
                charBefore[errorAdded] = generatedPattern[e];
                startSaIndexBefore[errorAdded] = startSaIndex[e];
                endSaIndexBefore[errorAdded] = endSaIndex[e];

                // Assign code for inserting the char = 0
                generatedPattern[e] = (unsigned char)(ALPHABET_SIZE + 1);
                startSaIndex[e] = BWTOccValue(bwt, startSaIndexBefore[errorAdded], 0) + 1;
                endSaIndex[e] = BWTOccValue(bwt, endSaIndexBefore[errorAdded] + 1, 0);

                length++;

            }

            // Search for pattern
            while (e>0 && startSaIndex[e] <= endSaIndex[e]) {
                e--;
                c = generatedPattern[e];
                startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
            }
            if (startSaIndex[e] <= endSaIndex[e]) {
                if (numOfSaGroup < maxSaIndexGroup) {
                    saIndexGroup[numOfSaGroup].startSaIndex = startSaIndex[0];
                    saIndexGroup[numOfSaGroup].numOfMatch = endSaIndex[0] - startSaIndex[0] + 1;
                    saIndexGroup[numOfSaGroup].length = length;
                    saIndexGroup[numOfSaGroup].error = numOfError;
                    numOfSaGroup++;
                } else {
                    if (exceedMaxSaIndexGroupPrinted == FALSE) {
                        fprintf(stderr, "Not enough memory to store all SA index groups!\n");
                        exceedMaxSaIndexGroupPrinted = TRUE;
                    }
                }
            }

        }

    }

    numOfSaGroup = BWTEliminateDupSaIndexGroup(saIndexGroup, numOfSaGroup);

    return numOfSaGroup;

}
/*

// This is a version making use of cachedSaIndex and is with bug
unsigned int BWTEditDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
                     SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned int maxSaIndexGroup) {

    // stack
    unsigned int startSaIndex[MAX_ARPROX_MATCH_LENGTH + 1];
    unsigned int endSaIndex[MAX_ARPROX_MATCH_LENGTH+1];
    unsigned char generatedPattern[MAX_ARPROX_MATCH_LENGTH+1];
    unsigned int errorPos[MAX_APPROX_MATCH_ERROR + 1];
    unsigned int errorPosLength[MAX_APPROX_MATCH_ERROR + 1];

    // for insertion
    unsigned int startSaIndexBefore[MAX_APPROX_MATCH_ERROR + 1];
    unsigned int endSaIndexBefore[MAX_APPROX_MATCH_ERROR+1];
    unsigned char charBefore[MAX_APPROX_MATCH_ERROR + 1];
    unsigned int saRangeIndexBefore[MAX_APPROX_MATCH_ERROR + 1];

    unsigned int numOfSaGroup;
    unsigned int numOfError;

    unsigned int c, e;
    unsigned int errorAdded;
    unsigned char maxPatternValue;
    unsigned int tempSaIndexLeft, tempSaIndexRight;

    unsigned int length;

    unsigned int exceedMaxSaIndexGroupPrinted = FALSE;

    unsigned int saRangeIndex;
    unsigned int saIndexRangeLength;

    #ifdef DEBUG
    if (maxError > keyLength) {
        fprintf(stderr, "BWTEditDistMatch() : maxError > keyLength!\n");
        exit(1);
    }
    if (keyLength > MAX_ARPROX_MATCH_LENGTH) {
        fprintf(stderr, "BWTEditDistMatch() : keyLength > MAX_ARPROX_MATCH_LENGTH!\n");
        exit(1);
    }
    if (maxError > MAX_APPROX_MATCH_ERROR) {
        fprintf(stderr, "BWTEditDistMatch() : maxError > MAX_APPROX_MATCH_ERROR!\n");
        exit(1);
    }
    #endif

    // Setup for looking up SA index range
    if (cachedSaIndexNumOfChar <= keyLength - maxError) {
        saIndexRangeLength = cachedSaIndexNumOfChar;
    } else {
        saIndexRangeLength = 0;
    }

    // Set this boundary case so that the last character of generated pattern can be located by backward search
    startSaIndex[keyLength] = 0;
    endSaIndex[keyLength] = bwt->textLength;
    maxPatternValue = (unsigned char)(ALPHABET_SIZE * 2);    // 0 to alphabetSize - 1 = edit; alphabetSize = delete; alphabetSize + 1 to alphabetsize * 2 = insert

    // Exact match
    numOfSaGroup = BWTBackwardSearch(convertedKey, keyLength, bwt, &tempSaIndexLeft, &tempSaIndexRight);
    if (numOfSaGroup > 0) {
        saIndexGroup[0].startSaIndex = tempSaIndexLeft;
        saIndexGroup[0].numOfMatch = tempSaIndexRight - tempSaIndexLeft + 1;
        saIndexGroup[0].length = keyLength;
        saIndexGroup[0].error = 0;
    }

    // With errors
    for (numOfError = 1; numOfError <= maxError; numOfError++) {

        memcpy(generatedPattern, convertedKey, keyLength);

        // Set initial state
        e = errorPos[1] = keyLength - 1;
        errorPosLength[1] = 1;
        generatedPattern[e] = (unsigned char)-1;
        errorAdded = 1;
        saRangeIndex = 0;

        // store before change data
        c = charBefore[1] = convertedKey[e];
        startSaIndexBefore[1] = bwt->cumulativeFreq[c] + 1;
        endSaIndexBefore[1] = bwt->cumulativeFreq[c+1];
        //startSaIndexBefore[1] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
        //endSaIndexBefore[1] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
        saRangeIndexBefore[1] = c;

        while (errorAdded > 0) {
            
            // Increment the last error
            e = errorPos[errorAdded];
            length = errorPosLength[errorAdded];
            if (length < cachedSaIndexNumOfChar + 1) {
                saRangeIndex >>= (cachedSaIndexNumOfChar - length + 1) * BIT_PER_CHAR;
            }
            generatedPattern[e]++;
            if (generatedPattern[e] == convertedKey[e]) {
                generatedPattern[e]++;
            }
            if (generatedPattern[e] > maxPatternValue) {

                // Cannot increment the last error; try moving the error to the left
                if (e > 1 || (e == 1 && numOfError - errorAdded == 0)) {    // only the last error can move to the leftmost position as insertion is not done on pattern boundary

                    // restore before change data
                    c = generatedPattern[e] = charBefore[errorAdded];
                    if (startSaIndexBefore[errorAdded] <= endSaIndexBefore[errorAdded] || length <= saIndexRangeLength) {
                        startSaIndex[e] = startSaIndexBefore[errorAdded];
                        endSaIndex[e] = endSaIndexBefore[errorAdded];
                        saRangeIndex = saRangeIndexBefore[errorAdded];
                    } else {
                        // Pattern suffix not found
                        errorAdded--;
                        continue;
                    }

                    // move error
                    errorPos[errorAdded]--;
                    e = errorPos[errorAdded];
                    generatedPattern[e] = (convertedKey[e] == 0);    // 1 if convertedKey = 0, 0 otherwise

                    // store before change data
                    c = charBefore[errorAdded] = convertedKey[e];

                    length++;
                    errorPosLength[errorAdded] = length;
                    if (length <= saIndexRangeLength) {
                        saRangeIndexBefore[errorAdded] = (saRangeIndex << BIT_PER_CHAR) | c;
                    } else {
                        if (length > saIndexRangeLength + 1) {
                            startSaIndexBefore[errorAdded] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                            endSaIndexBefore[errorAdded] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
                        } else {
                            startSaIndexBefore[errorAdded] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].startSaIndex, c) + 1;
                            endSaIndexBefore[errorAdded] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].endSaIndex + 1, c);
                        }
                    }

                } else {
                    // Cannot move error to the left
                    // restore before change data
                    generatedPattern[e] = charBefore[errorAdded];
                    //c = generatedPattern[e] = charBefore[errorAdded];
                    //startSaIndex[e] = startSaIndexBefore[errorAdded];
                    //endSaIndex[e] = endSaIndexBefore[errorAdded];

                    errorAdded--;
                    continue;
                }
            }

            if (generatedPattern[e] < ALPHABET_SIZE) {

                // edit
                c = generatedPattern[e];

                if (length <= saIndexRangeLength) {
                    saRangeIndex = (saRangeIndex << BIT_PER_CHAR) | c;
                } else {
                    if (length > saIndexRangeLength + 1) {
                         startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                        endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
                    } else {
                         startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].startSaIndex, c) + 1;
                        endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].endSaIndex + 1, c);
                    }
                }

            } else {
                if (generatedPattern[e] > ALPHABET_SIZE) {
                    if (e > 0 && (startSaIndexBefore[errorAdded] <= endSaIndexBefore[errorAdded]
                                  || length <= saIndexRangeLength)) {
                        // insert
                        c = generatedPattern[e] - ALPHABET_SIZE - 1;    // get back the character to be inserted
                        length++;
                        if (length <= saIndexRangeLength) {
                            saRangeIndex = (saRangeIndexBefore[errorAdded] << BIT_PER_CHAR) | c;
                        } else {
                            if (length > saIndexRangeLength + 1) {
                                 startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndexBefore[errorAdded], c) + 1;
                                endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndexBefore[errorAdded] + 1, c);
                            } else {
                                 startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndexBefore[errorAdded]].startSaIndex, c) + 1;
                                endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndexBefore[errorAdded]].endSaIndex + 1, c);
                            }
                        }
                    } else {
                        // do not insert on pattern boundary
                        generatedPattern[e] = maxPatternValue;
                        continue;
                    }
                } else {
                    // delete
                    startSaIndex[e] = startSaIndex[e+1];
                    endSaIndex[e] = endSaIndex[e+1];
                    //saRangeIndex >>= BIT_PER_CHAR;
                    length--;
                }
            }

            while (errorAdded < numOfError && startSaIndex[e] <= endSaIndex[e]) {

                // Add insertions
                errorAdded++;
                errorPos[errorAdded] = errorPos[errorAdded-1];
                errorPosLength[errorAdded] = length;

                // store before change data
                charBefore[errorAdded] = generatedPattern[e];
                startSaIndexBefore[errorAdded] = startSaIndex[e];
                endSaIndexBefore[errorAdded] = endSaIndex[e];
                saRangeIndexBefore[errorAdded] = saRangeIndex;

                // Assign code for inserting the char = 0
                generatedPattern[e] = (unsigned char)(ALPHABET_SIZE + 1);
                length++;
                if (length <= saIndexRangeLength) {
                    saRangeIndex <<= BIT_PER_CHAR;
                } else {
                    if (length > saIndexRangeLength + 1) {
                         startSaIndex[e] = BWTOccValue(bwt, startSaIndex[e], 0) + 1;
                        endSaIndex[e] = BWTOccValue(bwt, endSaIndex[e] + 1, 0);
                    } else {
                         startSaIndex[e] = BWTOccValue(bwt, cachedSaIndex[saRangeIndex].startSaIndex, 0) + 1;
                        endSaIndex[e] = BWTOccValue(bwt, cachedSaIndex[saRangeIndex].endSaIndex + 1, 0);
                    }
                }
            }

            // Search for pattern
            while (e>0 && (length <= saIndexRangeLength || startSaIndex[e] <= endSaIndex[e])) {
                e--;
                length++;
                c = generatedPattern[e];
                if (length <= saIndexRangeLength) {
                    saRangeIndex = (saRangeIndex << BIT_PER_CHAR) | c;
                } else {
                    if (length > saIndexRangeLength + 1) {
                        startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, startSaIndex[e+1], c) + 1;
                        endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, endSaIndex[e+1] + 1, c);
                    } else {
                         startSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].startSaIndex, c) + 1;
                        endSaIndex[e] = bwt->cumulativeFreq[c] + BWTOccValue(bwt, cachedSaIndex[saRangeIndex].endSaIndex + 1, c);
                    }
                }
            }
            if (startSaIndex[e] <= endSaIndex[e]) {
                if (numOfSaGroup < maxSaIndexGroup) {
                    saIndexGroup[numOfSaGroup].startSaIndex = startSaIndex[0];
                    saIndexGroup[numOfSaGroup].numOfMatch = endSaIndex[0] - startSaIndex[0] + 1;
                    saIndexGroup[numOfSaGroup].length = length;
                    saIndexGroup[numOfSaGroup].error = numOfError;
                    numOfSaGroup++;
                } else {
                    if (exceedMaxSaIndexGroupPrinted == FALSE) {
                        fprintf(stderr, "Not enough memory to store all SA index groups!\n");
                        exceedMaxSaIndexGroupPrinted = TRUE;
                    }
                }
            }

        }

    }

    numOfSaGroup = BWTEliminateDupSaIndexGroup(saIndexGroup, numOfSaGroup);

    return numOfSaGroup;

}
*/
unsigned long long BWTEliminateDupSaIndexGroup(SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned long long numOfSaGroup) {

    unsigned long long resultnumOfGroup;
    unsigned long long i;

    QSort(saIndexGroup, numOfSaGroup, sizeof(SaIndexGroupNew), SaIndexGroupStartSaIndexLengthErrorOrder);

    resultnumOfGroup = 0;
    for (i=0; i<numOfSaGroup; i++) {
        if (i > resultnumOfGroup) {
            saIndexGroup[resultnumOfGroup].startSaIndex = saIndexGroup[i].startSaIndex;
            saIndexGroup[resultnumOfGroup].numOfMatch = saIndexGroup[i].numOfMatch;
            saIndexGroup[resultnumOfGroup].length = saIndexGroup[i].length;
            saIndexGroup[resultnumOfGroup].error = saIndexGroup[i].error;
        }
        while (i<numOfSaGroup && saIndexGroup[i].startSaIndex == saIndexGroup[i+1].startSaIndex &&
                                  saIndexGroup[i].length == saIndexGroup[i+1].length) {
            // the group with the least error will be taken for the same saIndex range and length
            i++;
        }
        resultnumOfGroup++;
    }
    
    return resultnumOfGroup;
}

unsigned long long BWTSubPatternEditDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
                               const unsigned int maxError, const unsigned long long skip, const unsigned long long maxnumOfHit, 
                               HitListWithPosQueryLengthError* __restrict hitList, BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics,
                               const unsigned long long eliminateDuplicateStartingPos) {

    unsigned long long i, j;
    unsigned long long numOfSubPattern;
    unsigned long long advance;

    unsigned long long maxSaIndexGroup;
    SaIndexGroupWithLengthError *saIndexGroup;

    unsigned long long totalnumOfSaGroup, numOfSaGroup;

    unsigned long long nextFilteredChar;

    unsigned long long maxnumOfHitPerSaIndexGroup, numOfHit, totalnumOfHit;
    SaIndexList *tempSaIndexList1, *tempSaIndexList2;

    maxnumOfHitPerSaIndexGroup = 0;

    advance = skip + 1;
    numOfSubPattern = (keyLength - subPatternLength + advance) / advance;

    maxSaIndexGroup = BWTEditDistMaxSaIndexGroup(subPatternLength, maxError) * numOfSubPattern;
    saIndexGroup = (SaIndexGroupWithLengthError*) MMUnitAllocate(maxSaIndexGroup * sizeof(SaIndexGroupWithLengthError));
    totalnumOfSaGroup = 0;

    for (nextFilteredChar = 0; nextFilteredChar<keyLength && convertedKey[nextFilteredChar] < ALPHABET_SIZE; nextFilteredChar++) {
    }

    totalnumOfHit = 0;

    for (i=0; i<numOfSubPattern; i++) {

        if (nextFilteredChar < i*advance) {
            // Advance to the next 'N'
            for (; nextFilteredChar<i*advance && convertedKey[nextFilteredChar] >= ALPHABET_SIZE; nextFilteredChar++) {
            }
            for (; nextFilteredChar<keyLength && convertedKey[nextFilteredChar] < ALPHABET_SIZE; nextFilteredChar++) {
            }
        }

        if (i * advance + subPatternLength <= nextFilteredChar) {
            numOfSaGroup = BWTEditDistMatchOld(convertedKey + i*advance, subPatternLength, bwt, maxError, 
                                                saIndexGroup + totalnumOfSaGroup, maxSaIndexGroup - totalnumOfSaGroup);
            // Preserve position query information in error
            numOfHit = 0;
            for (j=0; j<numOfSaGroup; j++) {
                saIndexGroup[totalnumOfSaGroup + j].posQuery = i*advance;
                // determine the number of low occ SA and total number of hits
                numOfHit += saIndexGroup[totalnumOfSaGroup + j].numOfMatch;
            }
            if (numOfHit > maxnumOfHitPerSaIndexGroup) {
                maxnumOfHitPerSaIndexGroup = numOfHit;
            }

            totalnumOfSaGroup += numOfSaGroup;
            totalnumOfHit += numOfHit;
        }

    }

    if (totalnumOfHit > maxnumOfHit) {
        fprintf(stderr, "BWTSubPatternEditDistMatch(): Not enough memory allocated for hit generation!\n");
        exit(1);
    }

    if (totalnumOfSaGroup > 0) {

        tempSaIndexList1 = (SaIndexList*) MMUnitAllocate(maxnumOfHitPerSaIndexGroup * sizeof(SaIndexList));
        tempSaIndexList2 = (SaIndexList*) MMUnitAllocate(maxnumOfHitPerSaIndexGroup * sizeof(SaIndexList));

        // SA index group must be sorted by SA index + length + error
        QSort(saIndexGroup, totalnumOfSaGroup, sizeof(SaIndexGroupNew), SaIndexGroupStartSaIndexLengthErrorOrder);
        numOfHit = BWTTextPosition(bwt, (SaIndexGroupNew*)saIndexGroup, totalnumOfSaGroup, (HitList*)hitList,
                                      tempSaIndexList1, tempSaIndexList2,
                                      bwtSaRetrievalStatistics, eliminateDuplicateStartingPos);
    } else {
        numOfHit = 0;
    }

     MMUnitFree(saIndexGroup, maxSaIndexGroup * sizeof(SaIndexGroup));

     return numOfHit;

}

int BWTGappedDPDBvsQuery(BWT *bwt, const unsigned char *convertedKey, const int queryPatternLength, 
                         char* __restrict workingMemory, const int workingMemorySize, int* __restrict totalNumOfQueryPos,
                         const int matchScore, const int mismatchScore,
                         const int gapOpenScore, const int gapExtendScore,
                         const int cutoffScore,
                         BWTDPStatistics* __restrict bwtDPStatistics,
                         const int printProgressDepth) {


    #define NEG_INFINITY        -1073741824        // use -(2^30) to leave room for decreasing the value without overflow
    #define SCAN_DEPTH            4
    #define MIN_NUM_SOURCE_BIT    20
    #define MAX_NUM_SOURCE_BIT    22                // 2^22 * 4 bytes = 16M will be allocated for scanDepthP; 16M x 2 = 32M will be allocated for nonGapB
    #define DPCELL_BUFFER_SIZE    4194304            // 48M will be allocated for dpCell

    DPText* __restrict dpText;
    DPCell* __restrict dpCell;
    DPScanDepth* __restrict dpScanDepth;    // Logical position for depth = scanDepth
    int* __restrict nonGapB;                // Best score for scanDepth < depth <= nonGapDepth
    int M;
    int D;                // insert and delete wrt query

    int numSourceBit;
    int nonGapDepth;    // no gap when depth <= nonGapDepth

    int scoringMatrix[DNA_CHAR_SIZE][DNA_CHAR_SIZE];

    int cutoffScoreShifted, gapOpenScoreShifted, gapExtendScoreShifted;
    int minPositiveScore;

    int depth;
    int sourceBitMask;

    int numOfSaIndexGroup, numOfHit;

    int pathAccepted;

    int i, j, k;

    int workingMemoryUsed = 0;
    SaIndexGroupNew* __restrict saIndexGroup;
    unsigned int* __restrict posQueryList;
    int* __restrict numOfPosQuery;
    int maxNumOfPosQuery;
    unsigned long long minPosQuery;

    unsigned long long withAmbiguity;

    int dpScanDepthIndex, nonGapBIndex, dpCellIndex, lastDpCellIndex;
    int scanDepthScore;

    int dpScanDepthIndexUsed;
    int nonGapBIndexUsed;
    
    unsigned t, q;
    unsigned long long lastQPos, qPos;

    int insertScore;

    int maxNumOfSaIndexGroup, depthBit;

    // BWT DP statistics
    int maxDepth;
    int maxDPCell;
    int maxDPMemoryInWord;

    if ((1 - mismatchScore + matchScore - 1) / matchScore < SCAN_DEPTH) {
        fprintf(stderr, "Mismatch penalty is not large enough!\n");
        exit(1);
    }

    nonGapDepth = (1 - gapOpenScore - gapExtendScore + matchScore - 1) / matchScore - 1;
    if (nonGapDepth >= cutoffScore) {
        fprintf(stderr, "Cutoff score is too short for gaps!\n");
        exit(1);
    }
    if (nonGapDepth > CHAR_PER_WORD) {
        nonGapDepth = CHAR_PER_WORD;
    }

    numSourceBit = ceilLog2(queryPatternLength) - 4;    // heuristic; composition not more baised than 1/16 for 4 char substring
    if (numSourceBit > MAX_NUM_SOURCE_BIT) {
        numSourceBit = MAX_NUM_SOURCE_BIT;
    }
    if (numSourceBit < MIN_NUM_SOURCE_BIT) {
        numSourceBit = MIN_NUM_SOURCE_BIT;
    }

    HSPFillScoringMatrix(scoringMatrix, matchScore, mismatchScore, numSourceBit);
    scanDepthScore = SCAN_DEPTH * matchScore * (1 << numSourceBit);

    sourceBitMask = ~(ALL_ONE_MASK << numSourceBit);
    gapOpenScoreShifted = gapOpenScore * (1 << numSourceBit);
    gapExtendScoreShifted = gapExtendScore * (1 << numSourceBit);
    cutoffScoreShifted = cutoffScore * (1 << numSourceBit);
    minPositiveScore = 1 << numSourceBit;

    depthBit = ceilLog2(BWTDP_MAX_SUBSTRING_LENGTH);
    maxNumOfSaIndexGroup = 1 << (BITS_IN_WORD - depthBit);

    // Initialize variables to avoid compiler warnings
    init(dpScanDepthIndexUsed);
    init(nonGapBIndex);
    init(nonGapBIndexUsed);

    // allocate working memory
    dpText = (DPText*) MMUnitAllocate((BWTDP_MAX_SUBSTRING_LENGTH+1) * sizeof(DPText));
    dpScanDepth = (DPScanDepth*) MMUnitAllocate((1 << numSourceBit) * sizeof(DPScanDepth));
    nonGapB = (int*) MMUnitAllocate((1 << numSourceBit) * 2 * sizeof(int));
    dpCell = (DPCell*) MMUnitAllocate(DPCELL_BUFFER_SIZE * sizeof(DPCell));
    
    numOfSaIndexGroup = 0;
    *totalNumOfQueryPos = 0;
    numOfHit = 0;

    dpText[0].charBeingProcessed = 0;
    dpText[0].saIndexLeft[0] = 0;
    dpText[0].saIndexRight[0] = bwt->textLength;

    depth = 1;

    dpText[1].saIndexLeft[0] = bwt->cumulativeFreq[0] + 1;
    dpText[1].saIndexLeft[1] = bwt->cumulativeFreq[1] + 1;
    dpText[1].saIndexLeft[2] = bwt->cumulativeFreq[2] + 1;
    dpText[1].saIndexLeft[3] = bwt->cumulativeFreq[3] + 1;
    dpText[1].saIndexRight[0] = bwt->cumulativeFreq[1];
    dpText[1].saIndexRight[1] = bwt->cumulativeFreq[2];
    dpText[1].saIndexRight[2] = bwt->cumulativeFreq[3];
    dpText[1].saIndexRight[3] = bwt->cumulativeFreq[4];

    bwtDPStatistics->totalNode[0]++;
    maxDepth = 0;
    maxDPCell = 0;
    maxDPMemoryInWord = 0;
    if (SCAN_DEPTH > maxDepth) {
        maxDepth = SCAN_DEPTH;
    }

    dpText[depth].charBeingProcessed = -1;

    dpText[SCAN_DEPTH+1].dpCellIndex = 0;    // This is always zero because it is the start of another chunk of memory

    if (printProgressDepth != 0) {
        printf("DP Progress:\n");
    }

    while (depth > 0) {

        while (depth > 0 && depth < SCAN_DEPTH) {

            dpText[depth].charBeingProcessed++;
            while (dpText[depth].charBeingProcessed < ALPHABET_SIZE && 
                   dpText[depth].saIndexLeft[dpText[depth].charBeingProcessed] > dpText[depth].saIndexRight[dpText[depth].charBeingProcessed]) {
                dpText[depth].charBeingProcessed++;
            }

            if (dpText[depth].charBeingProcessed < ALPHABET_SIZE) {

                bwtDPStatistics->totalNode[depth]++;

                if (depth <= printProgressDepth) {
                    for (i=0; i<depth; i++) {
                        printf(" ");
                    }
                    printf("%d\n", dpText[depth].charBeingProcessed);
                    fflush(stdout);
                }

                // DP not done yet; just generate the string
                depth++;
                BWTAllOccValueTwoIndex(bwt, dpText[depth-1].saIndexLeft[dpText[depth-1].charBeingProcessed],
                                            dpText[depth-1].saIndexRight[dpText[depth-1].charBeingProcessed] + 1,                        
                                            dpText[depth].saIndexLeft,
                                            dpText[depth].saIndexRight);
                dpText[depth].saIndexLeft[0] += bwt->cumulativeFreq[0] + 1;
                dpText[depth].saIndexLeft[1] += bwt->cumulativeFreq[1] + 1;
                dpText[depth].saIndexLeft[2] += bwt->cumulativeFreq[2] + 1;
                dpText[depth].saIndexLeft[3] += bwt->cumulativeFreq[3] + 1;
                dpText[depth].saIndexRight[0] += bwt->cumulativeFreq[0];
                dpText[depth].saIndexRight[1] += bwt->cumulativeFreq[1];
                dpText[depth].saIndexRight[2] += bwt->cumulativeFreq[2];
                dpText[depth].saIndexRight[3] += bwt->cumulativeFreq[3];
                dpText[depth].charBeingProcessed = -1;

            } else {
                depth--;
            }

        }

        if (depth == SCAN_DEPTH) {

            dpText[depth].charBeingProcessed++;
            while (dpText[depth].charBeingProcessed < ALPHABET_SIZE && 
                   dpText[depth].saIndexLeft[dpText[depth].charBeingProcessed] > dpText[depth].saIndexRight[dpText[depth].charBeingProcessed]) {
                dpText[depth].charBeingProcessed++;
            }

            if (dpText[depth].charBeingProcessed < ALPHABET_SIZE) {

                bwtDPStatistics->totalNode[depth]++;

                // scan the query

                // pack generated substring into a word
                t = dpText[1].charBeingProcessed;
                for (i=2; i<=SCAN_DEPTH; i++) {
                    t <<= BITS_IN_BYTE;
                    t |= dpText[i].charBeingProcessed;
                }

                // pack query substring into a word
                q = convertedKey[queryPatternLength - 1];
                withAmbiguity = (convertedKey[queryPatternLength - 1] >= ALPHABET_SIZE);
                for (i=1; i<SCAN_DEPTH; i++) {
                    q <<= BITS_IN_BYTE;
                    q |= convertedKey[queryPatternLength - 1 - i];
                    withAmbiguity <<= BITS_IN_BYTE;
                    withAmbiguity |= (convertedKey[queryPatternLength - 1 - i] >= ALPHABET_SIZE);
                }

                dpScanDepthIndex = 0;
                for (i=queryPatternLength-1; i>nonGapDepth; i--) {

                    if ((withAmbiguity && !(withAmbiguity >> (BITS_IN_WORD - BITS_IN_BYTE))) || (!withAmbiguity && q == t)) {
                        if (dpScanDepthIndex >= minPositiveScore) {
                            fprintf(stderr, "BWTGappedDPDBvsQuery(): Not enough source bit!\n");
                            exit(1);
                        }
                        dpScanDepth[dpScanDepthIndex].P = i;
                        dpScanDepth[dpScanDepthIndex].withAmbiguity = (withAmbiguity != 0);
                        dpScanDepthIndex++;
                    }

                    q <<= BITS_IN_BYTE;
                    q |= convertedKey[i - SCAN_DEPTH];
                    withAmbiguity <<= BITS_IN_BYTE;
                    withAmbiguity |= (convertedKey[i - SCAN_DEPTH] >= ALPHABET_SIZE);

                }

                if (dpScanDepthIndex > 0) {

                    bwtDPStatistics->totalDPCell[depth] += dpScanDepthIndex;

                    dpScanDepthIndexUsed = dpScanDepthIndex;
                    if (dpScanDepthIndexUsed > maxDPCell) {
                        maxDPCell = dpScanDepthIndexUsed;
                    }
                    if (dpScanDepthIndexUsed > maxDPMemoryInWord) {
                        maxDPMemoryInWord = dpScanDepthIndexUsed;
                    }

                    depth++;
                    if (depth > maxDepth) {
                        maxDepth = depth;
                    }
                    BWTAllOccValueTwoIndex(bwt, dpText[depth-1].saIndexLeft[dpText[depth-1].charBeingProcessed],
                                                dpText[depth-1].saIndexRight[dpText[depth-1].charBeingProcessed] + 1,                        
                                                dpText[depth].saIndexLeft,
                                                dpText[depth].saIndexRight);
                    dpText[depth].saIndexLeft[0] += bwt->cumulativeFreq[0] + 1;
                    dpText[depth].saIndexLeft[1] += bwt->cumulativeFreq[1] + 1;
                    dpText[depth].saIndexLeft[2] += bwt->cumulativeFreq[2] + 1;
                    dpText[depth].saIndexLeft[3] += bwt->cumulativeFreq[3] + 1;
                    dpText[depth].saIndexRight[0] += bwt->cumulativeFreq[0];
                    dpText[depth].saIndexRight[1] += bwt->cumulativeFreq[1];
                    dpText[depth].saIndexRight[2] += bwt->cumulativeFreq[2];
                    dpText[depth].saIndexRight[3] += bwt->cumulativeFreq[3];
                    dpText[depth].charBeingProcessed = -1;

                } else {
                    bwtDPStatistics->rejectedPath++;
                    bwtDPStatistics->rejectedPathDepth += depth;
                    bwtDPStatistics->rejectedNode[depth]++;
                }

            } else {
                depth--;
            }

        }

        if (depth == SCAN_DEPTH + 1) {

            dpText[depth].charBeingProcessed++;
            while (dpText[depth].charBeingProcessed < ALPHABET_SIZE && 
                   dpText[depth].saIndexLeft[dpText[depth].charBeingProcessed] > dpText[depth].saIndexRight[dpText[depth].charBeingProcessed]) {
                dpText[depth].charBeingProcessed++;
            }

            if (dpText[depth].charBeingProcessed < ALPHABET_SIZE) {

                bwtDPStatistics->totalNode[depth]++;

                // scan for positive alignment and put to P

                nonGapBIndex = 0;
                for (i=0; i<dpScanDepthIndexUsed; i++) {

                    if (!dpScanDepth[i].withAmbiguity) {
                        t = dpText[depth].charBeingProcessed;
                        q = convertedKey[dpScanDepth[i].P - SCAN_DEPTH];
                        nonGapB[nonGapBIndex] = scanDepthScore + scoringMatrix[t][q] + i;
                    } else {
                        // Calculate score for the whole substring
                        t = dpText[1].charBeingProcessed;
                        q = convertedKey[dpScanDepth[i].P];
                        nonGapB[nonGapBIndex] = scoringMatrix[t][q] + i;
                        for (j=1; j<=SCAN_DEPTH && nonGapB[nonGapBIndex] >= minPositiveScore; j++) {
                            t = dpText[j+1].charBeingProcessed;
                            q = convertedKey[dpScanDepth[i].P - j];
                            nonGapB[nonGapBIndex] += scoringMatrix[t][q];
                        }
                    }
                    if (nonGapB[nonGapBIndex] >= minPositiveScore) {
                        nonGapBIndex++;
                    }
    
                }

                if (nonGapBIndex > 0) {

                    bwtDPStatistics->totalDPCell[depth] += nonGapBIndex;
                    if (dpScanDepthIndexUsed + nonGapBIndex > maxDPCell) {
                        maxDPCell = dpScanDepthIndexUsed + nonGapBIndex;
                    }
                    if (dpScanDepthIndexUsed + nonGapBIndex > maxDPMemoryInWord) {
                        maxDPMemoryInWord = dpScanDepthIndexUsed + nonGapBIndex;
                    }

                    dpText[depth+1].dpCellIndex = nonGapBIndex;

                    depth++;
                    if (depth > maxDepth) {
                        maxDepth = depth;
                    }
                    BWTAllOccValueTwoIndex(bwt, dpText[depth-1].saIndexLeft[dpText[depth-1].charBeingProcessed],
                                                dpText[depth-1].saIndexRight[dpText[depth-1].charBeingProcessed] + 1,                        
                                                dpText[depth].saIndexLeft,
                                                dpText[depth].saIndexRight);
                    dpText[depth].saIndexLeft[0] += bwt->cumulativeFreq[0] + 1;
                    dpText[depth].saIndexLeft[1] += bwt->cumulativeFreq[1] + 1;
                    dpText[depth].saIndexLeft[2] += bwt->cumulativeFreq[2] + 1;
                    dpText[depth].saIndexLeft[3] += bwt->cumulativeFreq[3] + 1;
                    dpText[depth].saIndexRight[0] += bwt->cumulativeFreq[0];
                    dpText[depth].saIndexRight[1] += bwt->cumulativeFreq[1];
                    dpText[depth].saIndexRight[2] += bwt->cumulativeFreq[2];
                    dpText[depth].saIndexRight[3] += bwt->cumulativeFreq[3];
                    dpText[depth].charBeingProcessed = -1;

                } else {
                    bwtDPStatistics->rejectedPath++;
                    bwtDPStatistics->rejectedPathDepth += depth;
                    bwtDPStatistics->rejectedNode[depth]++;
                }

            } else {
                depth--;
            }

        }

        while (depth > SCAN_DEPTH + 1 && depth < nonGapDepth) {

            dpText[depth].charBeingProcessed++;
            while (dpText[depth].charBeingProcessed < ALPHABET_SIZE && 
                   dpText[depth].saIndexLeft[dpText[depth].charBeingProcessed] > dpText[depth].saIndexRight[dpText[depth].charBeingProcessed]) {
                dpText[depth].charBeingProcessed++;
            }

            if (dpText[depth].charBeingProcessed < ALPHABET_SIZE) {

                bwtDPStatistics->totalNode[depth]++;

                // Gap is not considered

                t = dpText[depth].charBeingProcessed;

                nonGapBIndex = dpText[depth].dpCellIndex;
                for (i=dpText[depth-1].dpCellIndex; i<dpText[depth].dpCellIndex; i++) {
                    q = convertedKey[dpScanDepth[nonGapB[i] & sourceBitMask].P - depth + 1];
                    nonGapB[nonGapBIndex] = nonGapB[i] + scoringMatrix[t][q];
                    if (nonGapB[nonGapBIndex] >= minPositiveScore) {
                        nonGapBIndex++;
                    }
                }

                if (nonGapBIndex > dpText[depth].dpCellIndex) {

                    bwtDPStatistics->totalDPCell[depth] += nonGapBIndex - dpText[depth].dpCellIndex;
                    if (dpScanDepthIndexUsed + nonGapBIndex > maxDPCell) {
                        maxDPCell = dpScanDepthIndexUsed + nonGapBIndex;
                    }
                    if (dpScanDepthIndexUsed + nonGapBIndex > maxDPMemoryInWord) {
                        maxDPMemoryInWord = dpScanDepthIndexUsed + nonGapBIndex;
                    }

                    dpText[depth+1].dpCellIndex = nonGapBIndex;

                    // The following two variables will be used if depth = nonGapDepth - 1
                    nonGapBIndexUsed = nonGapBIndex;    // Use separate variable coz another chunk of memory is used
                    dpText[nonGapDepth].dpCellIndex = 0;

                    depth++;
                    if (depth > maxDepth) {
                        maxDepth = depth;
                    }
                    BWTAllOccValueTwoIndex(bwt, dpText[depth-1].saIndexLeft[dpText[depth-1].charBeingProcessed],
                                                dpText[depth-1].saIndexRight[dpText[depth-1].charBeingProcessed] + 1,                        
                                                dpText[depth].saIndexLeft,
                                                dpText[depth].saIndexRight);
                    dpText[depth].saIndexLeft[0] += bwt->cumulativeFreq[0] + 1;
                    dpText[depth].saIndexLeft[1] += bwt->cumulativeFreq[1] + 1;
                    dpText[depth].saIndexLeft[2] += bwt->cumulativeFreq[2] + 1;
                    dpText[depth].saIndexLeft[3] += bwt->cumulativeFreq[3] + 1;
                    dpText[depth].saIndexRight[0] += bwt->cumulativeFreq[0];
                    dpText[depth].saIndexRight[1] += bwt->cumulativeFreq[1];
                    dpText[depth].saIndexRight[2] += bwt->cumulativeFreq[2];
                    dpText[depth].saIndexRight[3] += bwt->cumulativeFreq[3];
                    dpText[depth].charBeingProcessed = -1;

                } else {
                    bwtDPStatistics->rejectedPath++;
                    bwtDPStatistics->rejectedPathDepth += depth;
                    bwtDPStatistics->rejectedNode[depth]++;
                }

            } else {
                depth--;
            }

        }

        if (depth == nonGapDepth) {

            dpText[depth].charBeingProcessed++;
            while (dpText[depth].charBeingProcessed < ALPHABET_SIZE && 
                   dpText[depth].saIndexLeft[dpText[depth].charBeingProcessed] > dpText[depth].saIndexRight[dpText[depth].charBeingProcessed]) {
                dpText[depth].charBeingProcessed++;
            }

            if (dpText[depth].charBeingProcessed < ALPHABET_SIZE) {

                bwtDPStatistics->totalNode[depth]++;

                // Gap is not considered

                t = dpText[depth].charBeingProcessed;

                dpCellIndex = 0;
                for (i=dpText[depth-1].dpCellIndex; i<nonGapBIndexUsed; i++) {
                    q = convertedKey[dpScanDepth[nonGapB[i] & sourceBitMask].P - depth + 1];
                    dpCell[dpCellIndex].B = nonGapB[i] + scoringMatrix[t][q];
                    if (dpCell[dpCellIndex].B >= minPositiveScore) {
                        dpCell[dpCellIndex].I = NEG_INFINITY;
                        dpCell[dpCellIndex].P = dpScanDepth[nonGapB[i] & sourceBitMask].P - depth + 1;
                        dpCellIndex++;
                    }
                }

                if (dpCellIndex > 0) {

                    bwtDPStatistics->totalDPCell[depth] += dpCellIndex;
                    if (dpScanDepthIndexUsed + nonGapBIndex + dpCellIndex > maxDPCell) {
                        maxDPCell = dpScanDepthIndexUsed + nonGapBIndex + dpCellIndex;
                    }
                    if (dpScanDepthIndexUsed + nonGapBIndex + dpCellIndex * (int)(sizeof(DPCell) / BYTES_IN_WORD) > maxDPMemoryInWord) {
                        maxDPMemoryInWord = dpScanDepthIndexUsed + nonGapBIndex + dpCellIndex * (sizeof(DPCell) / BYTES_IN_WORD);
                    }

                    depth++;
                    if (depth > maxDepth) {
                        maxDepth = depth;
                    }

                    dpText[depth].dpCellIndex = dpCellIndex;
                    dpText[depth].numOfDpCellSegment = (dpCellIndex + 4 - 1) / 4;        // each cell can at most produce 1 extra cell

                    BWTAllOccValueTwoIndex(bwt, dpText[depth-1].saIndexLeft[dpText[depth-1].charBeingProcessed],
                                                dpText[depth-1].saIndexRight[dpText[depth-1].charBeingProcessed] + 1,                        
                                                dpText[depth].saIndexLeft,
                                                dpText[depth].saIndexRight);
                    dpText[depth].saIndexLeft[0] += bwt->cumulativeFreq[0] + 1;
                    dpText[depth].saIndexLeft[1] += bwt->cumulativeFreq[1] + 1;
                    dpText[depth].saIndexLeft[2] += bwt->cumulativeFreq[2] + 1;
                    dpText[depth].saIndexLeft[3] += bwt->cumulativeFreq[3] + 1;
                    dpText[depth].saIndexRight[0] += bwt->cumulativeFreq[0];
                    dpText[depth].saIndexRight[1] += bwt->cumulativeFreq[1];
                    dpText[depth].saIndexRight[2] += bwt->cumulativeFreq[2];
                    dpText[depth].saIndexRight[3] += bwt->cumulativeFreq[3];
                    dpText[depth].charBeingProcessed = -1;

                } else {
                    bwtDPStatistics->rejectedPath++;
                    bwtDPStatistics->rejectedPathDepth += depth;
                    bwtDPStatistics->rejectedNode[depth]++;
                }

            } else {
                depth--;
            }

        }

        while (depth > nonGapDepth) {

            dpText[depth].charBeingProcessed++;
            while (dpText[depth].charBeingProcessed < ALPHABET_SIZE && 
                   dpText[depth].saIndexLeft[dpText[depth].charBeingProcessed] > dpText[depth].saIndexRight[dpText[depth].charBeingProcessed]) {
                dpText[depth].charBeingProcessed++;
            }

            if (dpText[depth].charBeingProcessed < ALPHABET_SIZE) {

                bwtDPStatistics->totalNode[depth]++;

                // DP

                D = NEG_INFINITY;
                lastQPos = ALL_ONE_MASK;
                t = dpText[depth].charBeingProcessed;

                pathAccepted = FALSE;

                dpCellIndex = dpText[depth].dpCellIndex;
                // check buffer
                if (dpCellIndex + dpText[depth].numOfDpCellSegment * 4 + dpText[depth].dpCellIndex - dpText[depth-1].dpCellIndex >= DPCELL_BUFFER_SIZE) {
                    // each segment can produce a maximum of 4 extra dpCell after 2 level
                    fprintf(stderr, "BWTGappedDPDBvsQuery(): Not enough DPCell buffer!\n");
                    exit(1);
                }
                dpText[depth+1].numOfDpCellSegment = 0;

                for (lastDpCellIndex=dpText[depth-1].dpCellIndex; lastDpCellIndex<dpText[depth].dpCellIndex; lastDpCellIndex++) {

                    qPos = dpCell[lastDpCellIndex].P - 1;

                    if (qPos != lastQPos - 1) {

                        dpText[depth+1].numOfDpCellSegment++;

                        // The cell on the left is non-positive; handle the insertion + deletion only
                        if (dpCell[lastDpCellIndex].I >= minPositiveScore || D >= minPositiveScore) {

                            dpCell[dpCellIndex].B = max(dpCell[lastDpCellIndex].I, D);
                            dpCell[dpCellIndex].I = dpCell[lastDpCellIndex].I + gapExtendScoreShifted;    // insert cannot follow delete
                            D = D + gapExtendScoreShifted;    // delete cannot follow insert 
                            dpCell[dpCellIndex].P = qPos + 1;

                            dpCellIndex++;

                        }
                    }
                
                    // DP - Match/mismatch
                    M = dpCell[lastDpCellIndex].B + scoringMatrix[t][convertedKey[qPos]];

                    if (lastDpCellIndex+1 < dpText[depth].dpCellIndex && dpCell[lastDpCellIndex+1].P == qPos) {
                        // insert score available
                        insertScore = dpCell[lastDpCellIndex+1].I;
                    } else {
                        insertScore = DP_NEG_INFINITY;
                    }

                    // Determine the best score
                    if (M >= insertScore && M >= D) {
                        // match score is maximum
                        dpCell[dpCellIndex].B = M;
                        dpCell[dpCellIndex].I = max(M + gapOpenScoreShifted, insertScore) + gapExtendScoreShifted;
                        D = max(M + gapOpenScoreShifted, D) + gapExtendScoreShifted;
                    } else {
                        dpCell[dpCellIndex].B = max(insertScore, D);
                        dpCell[dpCellIndex].I = insertScore + gapExtendScoreShifted;
                        D = D + gapExtendScoreShifted;
                    }

                    // check if the best score is positive
                    if (dpCell[dpCellIndex].B >= minPositiveScore) {

                        dpCell[dpCellIndex].P = qPos;
                        if (dpCell[dpCellIndex].B >= cutoffScoreShifted) {
                            pathAccepted = TRUE;    // Do not break out of the loop so that the whole row will be processed
                        }

                        dpCellIndex++;

                    }

                    lastQPos = qPos;

                    if ((lastDpCellIndex+1 >= dpText[depth].dpCellIndex || dpCell[lastDpCellIndex+1].P < (qPos-1)) && qPos > 0 && D >= minPositiveScore) {

                        // delete score is still positive;
                        qPos--;
                        dpCell[dpCellIndex].B = D;
                        dpCell[dpCellIndex].I = DP_NEG_INFINITY;
                        dpCell[dpCellIndex].P = qPos;
                        dpCellIndex++;
                        D = D + gapExtendScoreShifted;

                    }
                    if (qPos == 0) {
                        // Reached the start of query
                        if (dpCellIndex > 0 && dpCell[dpCellIndex-1].P == 0 && pathAccepted == FALSE) {
                            // No need to process the cell on the start of query
                            dpCellIndex--;
                        }
                        break;
                    }

                }

                dpText[depth+1].dpCellIndex = dpCellIndex;

                if (dpText[depth+1].dpCellIndex == dpText[depth].dpCellIndex) {

                    bwtDPStatistics->rejectedPath++;
                    bwtDPStatistics->rejectedPathDepth += depth;
                    bwtDPStatistics->rejectedNode[depth]++;

                } else {

                    bwtDPStatistics->totalDPCell[depth] += dpText[depth+1].dpCellIndex - dpText[depth].dpCellIndex;

                    if (dpScanDepthIndexUsed + nonGapBIndex + dpText[depth+1].dpCellIndex > maxDPCell) {
                        maxDPCell = dpScanDepthIndexUsed + nonGapBIndex + dpText[depth+1].dpCellIndex;
                    }
                    if (dpScanDepthIndexUsed + nonGapBIndex + dpText[depth+1].dpCellIndex * (int)(sizeof(DPCell) / BYTES_IN_WORD) > maxDPMemoryInWord) {
                        maxDPMemoryInWord = dpScanDepthIndexUsed + nonGapBIndex + dpText[depth+1].dpCellIndex * (sizeof(DPCell) / BYTES_IN_WORD);
                    }
                
                    if (pathAccepted || depth >= BWTDP_MAX_SUBSTRING_LENGTH) {

                        bwtDPStatistics->acceptedPath++;
                        bwtDPStatistics->acceptedPathDepth += depth;
                
                        if (workingMemorySize - workingMemoryUsed < sizeof(SaIndexGroupNew) + sizeof(int) * 2) {
                            fprintf(stderr, "BWTGappedDPDBvsQuery(): Not enough working memory allocated!\n");
                            exit(1);
                        }

                        saIndexGroup = (SaIndexGroupNew*)(workingMemory + workingMemoryUsed);
                        workingMemoryUsed += sizeof(SaIndexGroupNew);
                        saIndexGroup->startSaIndex = dpText[depth].saIndexLeft[dpText[depth].charBeingProcessed];
                        saIndexGroup->numOfMatch = dpText[depth].saIndexRight[dpText[depth].charBeingProcessed] - dpText[depth].saIndexLeft[dpText[depth].charBeingProcessed] + 1;
                        if (numOfSaIndexGroup >= maxNumOfSaIndexGroup) {
                            fprintf(stderr, "BWTGappedDPDBvsQuery(): Too many SA index group!\n");
                            exit(1);
                        }
                        saIndexGroup->info = (numOfSaIndexGroup << depthBit) | depth;
                        numOfHit += saIndexGroup->numOfMatch;
                        numOfSaIndexGroup++;

                        numOfPosQuery = (int*)(workingMemory + workingMemoryUsed);
                        workingMemoryUsed += sizeof(int);
                        posQueryList = (unsigned int*)(workingMemory + workingMemoryUsed);
                        maxNumOfPosQuery = (workingMemorySize - workingMemoryUsed) / sizeof(int);

                        // Scan for positive posQuery
                        *numOfPosQuery = 0;
                        minPosQuery = queryPatternLength + 1;
                        for (i=dpText[depth].dpCellIndex; i<dpText[depth+1].dpCellIndex; i++) {

                            if (*numOfPosQuery + 1 > maxNumOfPosQuery) {
                                fprintf(stderr, "BWTGappedDP2(): Not enough working memory allocated!\n");
                                exit(1);
                            }
                            // posQuery points to the character on the right of the ending position of alignment
                            qPos = dpScanDepth[dpCell[i].B & sourceBitMask].P + 1;
                            if (qPos <= minPosQuery) {
                                if (qPos < minPosQuery) {
                                    posQueryList[*numOfPosQuery] = qPos;
                                    (*numOfPosQuery)++;
                                    minPosQuery = qPos;
                                }
                            } else {
                                for (j=0; qPos < posQueryList[j]; j++) {    // some posQuery in the list should <= posQuery
                                }
                                if (qPos > posQueryList[j]) {
                                    k = *numOfPosQuery;
                                    (*numOfPosQuery)++;
                                    for (; k > j; k--) {
                                        posQueryList[k] = posQueryList[k-1];
                                    }
                                    posQueryList[j] = qPos;
                                }
                            }
                        }
                        saIndexGroup->posQuery = posQueryList[0];    // For detecting diagonal hits only
                        *totalNumOfQueryPos += *numOfPosQuery;
                        workingMemoryUsed += *numOfPosQuery * sizeof(int);

                    } else {

                        depth++;
                        if (depth > maxDepth) {
                            maxDepth = depth;
                        }
                        BWTAllOccValueTwoIndex(bwt, dpText[depth-1].saIndexLeft[dpText[depth-1].charBeingProcessed],
                                                    dpText[depth-1].saIndexRight[dpText[depth-1].charBeingProcessed] + 1,                        
                                                    dpText[depth].saIndexLeft,
                                                    dpText[depth].saIndexRight);
                        dpText[depth].saIndexLeft[0] += bwt->cumulativeFreq[0] + 1;
                        dpText[depth].saIndexLeft[1] += bwt->cumulativeFreq[1] + 1;
                        dpText[depth].saIndexLeft[2] += bwt->cumulativeFreq[2] + 1;
                        dpText[depth].saIndexLeft[3] += bwt->cumulativeFreq[3] + 1;
                        dpText[depth].saIndexRight[0] += bwt->cumulativeFreq[0];
                        dpText[depth].saIndexRight[1] += bwt->cumulativeFreq[1];
                        dpText[depth].saIndexRight[2] += bwt->cumulativeFreq[2];
                        dpText[depth].saIndexRight[3] += bwt->cumulativeFreq[3];
                        dpText[depth].charBeingProcessed = -1;

                    }

                }

            } else {
                depth--;
            }

        }

    }

    if (maxDepth > bwtDPStatistics->maxDepth) {
        bwtDPStatistics->maxDepth = maxDepth;
    }
    if (maxDPCell > bwtDPStatistics->maxDPCell) {
        bwtDPStatistics->maxDPCell = maxDPCell;
    }
    if (maxDPMemoryInWord > bwtDPStatistics->maxDPMemoryInWord) {
        bwtDPStatistics->maxDPMemoryInWord = maxDPMemoryInWord;
    }
    bwtDPStatistics->totalMaxDepth += maxDepth;
    bwtDPStatistics->totalMaxDPCell += maxDPCell;
    bwtDPStatistics->totalMaxDPMemoryInWord += maxDPMemoryInWord;

    if (printProgressDepth != 0) {
        printf("\n");
    }

    // free working memory
    MMUnitFree(dpText, (BWTDP_MAX_SUBSTRING_LENGTH+1) * sizeof(DPText));
    MMUnitFree(dpScanDepth, (1 << numSourceBit) * sizeof(DPScanDepth));
    MMUnitFree(nonGapB, (1 << numSourceBit) * 2 * sizeof(int));
    MMUnitFree(dpCell, DPCELL_BUFFER_SIZE * sizeof(DPCell));

    return numOfSaIndexGroup;

}


void BWTPrintDPStatistics(FILE * outFile, const BWTDPStatistics* bwtDPStatistics) {

    int i;

    Socketfprintf(outFile, "Depth     No. of nodes     No. of rejected nodes    No. of DP Cells\n");
    for (i=0; i<=BWTDP_MAX_SUBSTRING_LENGTH && bwtDPStatistics->totalNode[i] > 0; i++) {
        Socketfprintf(outFile, " %3d   %10lld            %10lld             %10lld\n", i, bwtDPStatistics->totalNode[i], bwtDPStatistics->rejectedNode[i], bwtDPStatistics->totalDPCell[i]);
    }

    Socketfprintf(outFile, "Max depth               : %d (%d)\n", bwtDPStatistics->maxDepth, bwtDPStatistics->totalMaxDepth);
    Socketfprintf(outFile, "Max no. of DP Cell      : %d (%d)\n", bwtDPStatistics->maxDPCell, bwtDPStatistics->totalMaxDPCell);
    Socketfprintf(outFile, "Max amount of DP Memory : %d (%d)\n", bwtDPStatistics->maxDPMemoryInWord * BYTES_IN_WORD, bwtDPStatistics->totalMaxDPMemoryInWord * BYTES_IN_WORD);
    Socketfprintf(outFile, "No. of accepted path    : %lld (%lld)\n", bwtDPStatistics->acceptedPath, bwtDPStatistics->acceptedPathDepth);
    Socketfprintf(outFile, "No. of rejected path    : %lld (%lld)\n", bwtDPStatistics->rejectedPath, bwtDPStatistics->rejectedPathDepth);
    Socketfprintf(outFile, "\n");

}


// saIndexGroup must be sorted in startSaIndex + length + error order; there must be no overlapping groups except that one group can completely enclose another

unsigned long long BWTTextPosition(const BWT *bwt, const SaIndexGroupNew *saIndexGroup, const unsigned long long numOfSaIndexGroups, 
                    HitList* __restrict hitList, 
                    SaIndexList* __restrict tempSaIndexList1, SaIndexList* __restrict tempSaIndexList2, 
                    BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics, const unsigned long long eliminateDuplicateStartingPos) {

    unsigned long long i, j, k;
    unsigned long long startIndex;
    unsigned long long retreivedStartSaIndex, retrievedEndSaIndexPlus1;
    unsigned long long numOfMatch;
    unsigned long long retreivedTextPositionIndex;
    unsigned long long left, rightPlus1;
    unsigned long long totalnumOfHit, numOfHitRetrieved;

    totalnumOfHit = 0;

    i = 0;
    while(i < numOfSaIndexGroups) {

        startIndex = i;
        retreivedStartSaIndex = saIndexGroup[startIndex].startSaIndex;
        retreivedTextPositionIndex = totalnumOfHit;
        numOfMatch = saIndexGroup[startIndex].numOfMatch;
        retrievedEndSaIndexPlus1 = retreivedStartSaIndex + numOfMatch;
        numOfHitRetrieved = 0;

        // Retrieve text positions for the SA index group with the shortest length for SA index groups with the same startSaIndex
        if (tempSaIndexList1 != NULL) {
            for (k=0; k<numOfMatch; k++) {
                tempSaIndexList1[k].saIndex = saIndexGroup[startIndex].startSaIndex + k;
                tempSaIndexList1[k].textPositionIndex = totalnumOfHit + k;
                hitList[totalnumOfHit + k].info = saIndexGroup[startIndex].info;
            }
            BWTDecodeTextPosition(bwt, tempSaIndexList1, tempSaIndexList2, numOfMatch, 0, ALL_ONE_MASK_64,
                                  hitList, bwtSaRetrievalStatistics);
        } else {
            for (k=0; k<numOfMatch; k++) {
                hitList[totalnumOfHit + k].posText = BWTSaValue(bwt, saIndexGroup[startIndex].startSaIndex + k);
                hitList[totalnumOfHit + k].info = saIndexGroup[startIndex].info;
            }
            bwtSaRetrievalStatistics->bwtSaRetrieved += numOfMatch;
        }
        numOfHitRetrieved += numOfMatch;
        totalnumOfHit += numOfMatch;

        i++;

        if (numOfHitRetrieved > 0) {
            
            // fill in text positions for enclosed groups

            while (i < numOfSaIndexGroups && saIndexGroup[i].startSaIndex + saIndexGroup[i].numOfMatch <= retrievedEndSaIndexPlus1) {
                left = saIndexGroup[i].startSaIndex;
                rightPlus1 = left + saIndexGroup[i].numOfMatch;
                if (left < retreivedStartSaIndex) {
                    fprintf(stderr, "BWTTextPosition(): SA index group not sorted!\n");
                    exit(1);
                }
                for (k=left; k<rightPlus1; k++) {
                    hitList[totalnumOfHit].posText = hitList[retreivedTextPositionIndex + k - retreivedStartSaIndex].posText;
                    hitList[totalnumOfHit].info = saIndexGroup[i].info;
                    totalnumOfHit++;
                }
                bwtSaRetrievalStatistics->saDuplicated += saIndexGroup[i].numOfMatch;
                i++;
            }

        }

    }

    if (eliminateDuplicateStartingPos) {

        QSort(hitList, totalnumOfHit, sizeof(HitList), HitListPosTextErrorLengthOrder);

        i = 0;
        j = 0;
        while (i < totalnumOfHit) {
            if (i > j) {
                hitList[j].posText = hitList[i].posText;
                hitList[j].info = hitList[i].info;
            }
            i++;
            while (i < totalnumOfHit && hitList[i].posText == hitList[j].posText) {
                i++;
            }
            j++;
        }
        totalnumOfHit = j;
    }


    return totalnumOfHit;

}

// Input SA index is in evenIterationSaIndexList and must be sorted in increasing SA index order
// Input SA index must not be zero
// If maxnumOfIteration is set, undecoded SA index is in evenIterationSaIndexList if maxnumOfIteration is even or 
//                                                            oddIterationSaIndexList if maxnumOfIteration is odd
// The number of undecoded SA is returned
unsigned long long BWTDecodeTextPosition(const BWT *bwt, SaIndexList* __restrict evenIterationSaIndexList, SaIndexList* __restrict oddIterationSaIndexList,
                          const unsigned long long numOfSaIndex, const unsigned long long numOfIterationProcessed, const unsigned long long maxnumOfIteration, 
                          HitList* __restrict hitList, BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics) {

    #define    CHAR_GAP_THRESHOLD    16
    #define    WORD_GAP_THRESHOLD    128

    unsigned long long iteration;
    unsigned long long numOfSaToDecode;
    unsigned long long numOfUndecodedSa, numOfSaProcessed;
    unsigned long long ALIGN_16 occValue[ALPHABET_SIZE];
    unsigned long long undecodedSaBwtCharCount[ALPHABET_SIZE + 1];

    unsigned long long firstSaIndex, lastSaIndexInPrevGroup;
    unsigned long long saIndexAdjustment;
    unsigned long long numOfSaToProcess, maxnumOfSaToProcess, saIndexLimit;
    unsigned long long lastBwtWordOffset, lastBwtBitOffset;
    unsigned long long bwtWordOffset, bwtBitOffset;

    unsigned long long i, c, sum, bwtChar;
    unsigned long long firstBwtChar;
    unsigned long long saIndex;


    firstBwtChar = (bwt->bwtCode[0] & 0x77777777) >> (BITS_IN_WORD - BIT_PER_CHAR);

    // First decode SA index with saIndex % saInterval == 0
    if (numOfIterationProcessed == 0) {
        numOfSaProcessed = 0;
        numOfUndecodedSa = 0;
        while (numOfSaProcessed < numOfSaIndex) {
            if (evenIterationSaIndexList[numOfSaProcessed].saIndex % bwt->saInterval != 0) {
                // Pack undecoded SA to the front
                if (numOfSaProcessed > numOfUndecodedSa) {
                    evenIterationSaIndexList[numOfUndecodedSa].saIndex = evenIterationSaIndexList[numOfSaProcessed].saIndex;
                    evenIterationSaIndexList[numOfUndecodedSa].textPositionIndex = evenIterationSaIndexList[numOfSaProcessed].textPositionIndex;
                }
                numOfUndecodedSa++;
            } else {
                // decode the text position
                hitList[evenIterationSaIndexList[numOfSaProcessed].textPositionIndex].posText = numOfIterationProcessed
                                                                                                   + bwt->saValue[evenIterationSaIndexList[numOfSaProcessed].saIndex / bwt->saInterval];
                bwtSaRetrievalStatistics->bwtSaRetrieved++;
            }
            numOfSaProcessed++;
        }
    } else {
        numOfUndecodedSa = numOfSaIndex;
    }


    // This while loop processes 2 iterations at a time
    // The codes for both odd and even iterations are the same except for the source and target SA index list (to avoid aliasing)

    iteration = 0;

    while (numOfUndecodedSa > 0 && iteration < maxnumOfIteration) {

        // Odd iteration; input in evenIterationSaIndexList; output to oddIterationSaIndexList

        numOfSaProcessed = 0;
        lastSaIndexInPrevGroup = 0;
        occValue[0] = 0;
        occValue[1] = 0;
        occValue[2] = 0;
        occValue[3] = 0;
        occValue[firstBwtChar] = 1;
        undecodedSaBwtCharCount[0] = 0;
        undecodedSaBwtCharCount[1] = 0;
        undecodedSaBwtCharCount[2] = 0;
        undecodedSaBwtCharCount[3] = 0;
        undecodedSaBwtCharCount[4] = 0;
        lastBwtWordOffset = 0;
        lastBwtBitOffset = BIT_PER_CHAR;

        while (numOfSaProcessed < numOfUndecodedSa) {

            firstSaIndex = evenIterationSaIndexList[numOfSaProcessed].saIndex;

            // keep SA index in a group either < or > inverseSa0
            if (firstSaIndex < bwt->inverseSa0) {
                saIndexLimit = bwt->inverseSa0;
                saIndexAdjustment = 0;
            } else {
                if (firstSaIndex > bwt->inverseSa0) {
                    saIndexLimit = ALL_ONE_MASK_64;
                    saIndexAdjustment = (unsigned int)-1;    // $ in BWT is not encoded
                } else {
                    // Decode the SA when SA index = inverseSa0
                    hitList[evenIterationSaIndexList[numOfSaProcessed].textPositionIndex].posText = numOfIterationProcessed + iteration;
                    evenIterationSaIndexList[numOfSaProcessed].saIndex = 0;    // mark as decoded
                    lastSaIndexInPrevGroup = bwt->inverseSa0;
                    numOfSaProcessed++;
                    bwtSaRetrievalStatistics->bwtSaRetrieved++;
                    continue;
                }
            }

            // Find a group of near consecutive SA index

            numOfSaToProcess = 1;
            maxnumOfSaToProcess = numOfUndecodedSa - numOfSaProcessed;

            while (numOfSaToProcess < maxnumOfSaToProcess &&
                   evenIterationSaIndexList[numOfSaProcessed + numOfSaToProcess].saIndex
                    - evenIterationSaIndexList[numOfSaProcessed + numOfSaToProcess - 1].saIndex < CHAR_GAP_THRESHOLD &&
                   evenIterationSaIndexList[numOfSaProcessed + numOfSaToProcess].saIndex < saIndexLimit) {
                numOfSaToProcess++;
            }

            // Find the occurrence count of the first SA index (up to the character before firstSaIndex)
            bwtWordOffset = (firstSaIndex + saIndexAdjustment) / CHAR_PER_WORD;
            bwtBitOffset = ((firstSaIndex + saIndexAdjustment) % CHAR_PER_WORD) * BIT_PER_CHAR;

            if (firstSaIndex - lastSaIndexInPrevGroup >= WORD_GAP_THRESHOLD ||
                (lastSaIndexInPrevGroup <= bwt->inverseSa0 && firstSaIndex > bwt->inverseSa0)) {

                // retrieve occurrence count without making use of last occurrence value information
                BWTAllOccValue(bwt, firstSaIndex, occValue);    // SA index is adjusted in BWTALLOccValue()

            } else {

                // Advance by retrieving BWT only
                // The BWT character pointed by wordOffset+bitOffset is not counted

                #ifdef DEBUG
                if (bwtWordOffset < lastBwtWordOffset) {
                    fprintf(stderr, "BWTDecodeTextPosition(): bwtIndex < lastBwtIndex!\n");
                    exit(1);
                }
                #endif

                // first advance to the next word boundary
                sum = 0;
                // The whole word = lastBwtWordOffset is counted
                if (lastBwtBitOffset > 0) {
                    c = bwt->bwtCode[lastBwtWordOffset] & 0x77777777 & (ALL_ONE_MASK >> lastBwtBitOffset);
                    sum += (bwt->decodeTable[c >> 16] + bwt->decodeTable[c & 0xFFFF]) - lastBwtBitOffset / BIT_PER_CHAR;
                    lastBwtBitOffset = 0;
                    lastBwtWordOffset++;
                }
                while (lastBwtWordOffset <= bwtWordOffset) {
                    c = bwt->bwtCode[lastBwtWordOffset] & 0x77777777;
                    sum += bwt->decodeTable[c >> 16] + bwt->decodeTable[c & 0xFFFF];
                    lastBwtWordOffset++;
                }
                // Subtract the remaining in bwtWordOffset
                c = bwt->bwtCode[bwtWordOffset] & 0x77777777 & (ALL_ONE_MASK >> bwtBitOffset);
                sum -= (bwt->decodeTable[c >> 16] + bwt->decodeTable[c & 0xFFFF]) - bwtBitOffset / BIT_PER_CHAR;

                occValue[0] += sum & 0x000000FF;    sum >>= 8;
                occValue[1] += sum & 0x000000FF;    sum >>= 8;
                occValue[2] += sum & 0x000000FF;    sum >>= 8;
                occValue[3] += sum;

            }

            // Decode the SA range for the group
            lastSaIndexInPrevGroup = evenIterationSaIndexList[numOfSaProcessed + numOfSaToProcess - 1].saIndex;
            numOfSaToDecode = lastSaIndexInPrevGroup - firstSaIndex + 1;

            saIndex = firstSaIndex;
            c = (bwt->bwtCode[bwtWordOffset] & 0x77777777) << bwtBitOffset;

            for (i=0; i<numOfSaToDecode; i++) {

                bwtChar = c >> (BITS_IN_WORD - BIT_PER_CHAR);
                occValue[bwtChar]++;

                if (saIndex == evenIterationSaIndexList[numOfSaProcessed].saIndex) {
                    if ((occValue[bwtChar] + bwt->cumulativeFreq[bwtChar]) % bwt->saInterval != 0) {
                        evenIterationSaIndexList[numOfSaProcessed].textPositionIndex |= bwtChar << (BITS_IN_WORD - BIT_PER_CHAR);
                        evenIterationSaIndexList[numOfSaProcessed].saIndex = occValue[bwtChar] + bwt->cumulativeFreq[bwtChar];
                        undecodedSaBwtCharCount[bwtChar + 1]++;
                    } else {
                        // decode the text position
                        hitList[evenIterationSaIndexList[numOfSaProcessed].textPositionIndex].posText = numOfIterationProcessed + iteration + 1
                                                                                                           + bwt->saValue[(occValue[bwtChar] + bwt->cumulativeFreq[bwtChar]) / bwt->saInterval];
                        evenIterationSaIndexList[numOfSaProcessed].saIndex = 0;    // mark as decoded
                        bwtSaRetrievalStatistics->bwtSaRetrieved++;
                    }
                    numOfSaProcessed++;
                }

                saIndex++;
                c <<= BIT_PER_CHAR;
                bwtBitOffset += BIT_PER_CHAR;
                if (bwtBitOffset >= BITS_IN_WORD) {
                    bwtBitOffset = 0;
                    bwtWordOffset++;
                    c = bwt->bwtCode[bwtWordOffset] & 0x77777777;
                }

            }

            lastBwtWordOffset = bwtWordOffset;
            lastBwtBitOffset = bwtBitOffset;

        }

        undecodedSaBwtCharCount[2] += undecodedSaBwtCharCount[1];
        undecodedSaBwtCharCount[3] += undecodedSaBwtCharCount[2];
        undecodedSaBwtCharCount[4] += undecodedSaBwtCharCount[3];

        // Stable bucket sort the decoded SA into olddIterationSaIndexList
        for (i=0; i<numOfSaProcessed; i++) {
            if (evenIterationSaIndexList[i].saIndex != 0) {
                bwtChar = evenIterationSaIndexList[i].textPositionIndex >> (BITS_IN_WORD - BIT_PER_CHAR);
                oddIterationSaIndexList[undecodedSaBwtCharCount[bwtChar]].saIndex = evenIterationSaIndexList[i].saIndex;
                oddIterationSaIndexList[undecodedSaBwtCharCount[bwtChar]].textPositionIndex = evenIterationSaIndexList[i].textPositionIndex & (ALL_ONE_MASK >> BIT_PER_CHAR);
                undecodedSaBwtCharCount[bwtChar]++;
            }
        }

        iteration++;
        numOfUndecodedSa = undecodedSaBwtCharCount[4];

        if (numOfUndecodedSa == 0 || iteration >= maxnumOfIteration) {
            break;
        }


        // Even iteration; input in oddIterationSaIndexList; output to evenIterationSaIndexList

        numOfSaProcessed = 0;
        lastSaIndexInPrevGroup = 0;
        occValue[0] = 0;
        occValue[1] = 0;
        occValue[2] = 0;
        occValue[3] = 0;
        occValue[firstBwtChar] = 1;
        undecodedSaBwtCharCount[0] = 0;
        undecodedSaBwtCharCount[1] = 0;
        undecodedSaBwtCharCount[2] = 0;
        undecodedSaBwtCharCount[3] = 0;
        undecodedSaBwtCharCount[4] = 0;
        lastBwtWordOffset = 0;
        lastBwtBitOffset = BIT_PER_CHAR;

        while (numOfSaProcessed < numOfUndecodedSa) {

            firstSaIndex = oddIterationSaIndexList[numOfSaProcessed].saIndex;

            // keep SA index in a group either < or > inverseSa0
            if (firstSaIndex < bwt->inverseSa0) {
                saIndexLimit = bwt->inverseSa0;
                saIndexAdjustment = 0;
            } else {
                if (firstSaIndex > bwt->inverseSa0) {
                    saIndexLimit = ALL_ONE_MASK;
                    saIndexAdjustment = (unsigned int)-1;    // $ in BWT is not encoded
                } else {
                    // Decode the SA when SA index = inverseSa0
                    hitList[oddIterationSaIndexList[numOfSaProcessed].textPositionIndex].posText = numOfIterationProcessed + iteration;
                    oddIterationSaIndexList[numOfSaProcessed].saIndex = 0;    // mark as decoded
                    lastSaIndexInPrevGroup = bwt->inverseSa0;
                    numOfSaProcessed++;
                    bwtSaRetrievalStatistics->bwtSaRetrieved++;
                    continue;
                }
            }

            // Find a group of near consecutive SA index

            numOfSaToProcess = 1;
            maxnumOfSaToProcess = numOfUndecodedSa - numOfSaProcessed;

            while (numOfSaToProcess < maxnumOfSaToProcess &&
                   oddIterationSaIndexList[numOfSaProcessed + numOfSaToProcess].saIndex
                    - oddIterationSaIndexList[numOfSaProcessed + numOfSaToProcess - 1].saIndex < CHAR_GAP_THRESHOLD &&
                   oddIterationSaIndexList[numOfSaProcessed + numOfSaToProcess].saIndex < saIndexLimit) {
                numOfSaToProcess++;
            }

            // Find the occurrence count of the first SA index (up to the character before firstSaIndex)
            bwtWordOffset = (firstSaIndex + saIndexAdjustment) / CHAR_PER_WORD;
            bwtBitOffset = ((firstSaIndex + saIndexAdjustment) % CHAR_PER_WORD) * BIT_PER_CHAR;

            if (firstSaIndex - lastSaIndexInPrevGroup >= WORD_GAP_THRESHOLD ||
                (lastSaIndexInPrevGroup <= bwt->inverseSa0 && firstSaIndex > bwt->inverseSa0)) {

                // retrieve occurrence count without making use of last occurrence value information
                BWTAllOccValue(bwt, firstSaIndex, occValue);    // SA index is adjusted in BWTALLOccValue()

            } else {

                // Advance by retrieving BWT only
                // The BWT character pointed by wordOffset+bitOffset is not counted

                #ifdef DEBUG
                if (bwtWordOffset < lastBwtWordOffset) {
                    fprintf(stderr, "BWTDecodeTextPosition(): bwtIndex < lastBwtIndex!\n");
                    exit(1);
                }
                #endif

                // first advance to the next word boundary
                sum = 0;
                // The whole word = lastBwtWordOffset is counted
                if (lastBwtBitOffset > 0) {
                    c = bwt->bwtCode[lastBwtWordOffset] & 0x77777777 & (ALL_ONE_MASK >> lastBwtBitOffset);
                    sum += (bwt->decodeTable[c >> 16] + bwt->decodeTable[c & 0xFFFF]) - lastBwtBitOffset / BIT_PER_CHAR;
                    lastBwtBitOffset = 0;
                    lastBwtWordOffset++;
                }
                while (lastBwtWordOffset <= bwtWordOffset) {
                    c = bwt->bwtCode[lastBwtWordOffset] & 0x77777777;
                    sum += bwt->decodeTable[c >> 16] + bwt->decodeTable[c & 0xFFFF];
                    lastBwtWordOffset++;
                }
                // Subtract the remaining in bwtWordOffset
                c = bwt->bwtCode[bwtWordOffset] & 0x77777777 & (ALL_ONE_MASK >> bwtBitOffset);
                sum -= (bwt->decodeTable[c >> 16] + bwt->decodeTable[c & 0xFFFF]) - bwtBitOffset / BIT_PER_CHAR;

                occValue[0] += sum & 0x000000FF;    sum >>= 8;
                occValue[1] += sum & 0x000000FF;    sum >>= 8;
                occValue[2] += sum & 0x000000FF;    sum >>= 8;
                occValue[3] += sum;

            }

            // Decode the SA range for the group to bufferDecodedSaIndex
            lastSaIndexInPrevGroup = oddIterationSaIndexList[numOfSaProcessed + numOfSaToProcess - 1].saIndex;
            numOfSaToDecode = lastSaIndexInPrevGroup - firstSaIndex + 1;

            saIndex = firstSaIndex;
            c = (bwt->bwtCode[bwtWordOffset] & 0x77777777) << bwtBitOffset;

            for (i=0; i<numOfSaToDecode; i++) {

                bwtChar = c >> (BITS_IN_WORD - BIT_PER_CHAR);
                occValue[bwtChar]++;

                if (saIndex == oddIterationSaIndexList[numOfSaProcessed].saIndex) {
                    if ((occValue[bwtChar] + bwt->cumulativeFreq[bwtChar]) % bwt->saInterval != 0) {
                        oddIterationSaIndexList[numOfSaProcessed].textPositionIndex |= bwtChar << (BITS_IN_WORD - BIT_PER_CHAR);
                        oddIterationSaIndexList[numOfSaProcessed].saIndex = occValue[bwtChar] + bwt->cumulativeFreq[bwtChar];
                        undecodedSaBwtCharCount[bwtChar + 1]++;
                    } else {
                        // decode the text position
                        hitList[oddIterationSaIndexList[numOfSaProcessed].textPositionIndex].posText = numOfIterationProcessed + iteration + 1
                                                                                                   + bwt->saValue[(occValue[bwtChar] + bwt->cumulativeFreq[bwtChar]) / bwt->saInterval];
                        oddIterationSaIndexList[numOfSaProcessed].saIndex = 0;    // mark as decoded
                        bwtSaRetrievalStatistics->bwtSaRetrieved++;
                    }
                    numOfSaProcessed++;
                }

                saIndex++;
                c <<= BIT_PER_CHAR;
                bwtBitOffset += BIT_PER_CHAR;
                if (bwtBitOffset >= BITS_IN_WORD) {
                    bwtBitOffset = 0;
                    bwtWordOffset++;
                    c = bwt->bwtCode[bwtWordOffset] & 0x77777777;
                }

            }

            lastBwtWordOffset = bwtWordOffset;
            lastBwtBitOffset = bwtBitOffset;

        }

        undecodedSaBwtCharCount[2] += undecodedSaBwtCharCount[1];
        undecodedSaBwtCharCount[3] += undecodedSaBwtCharCount[2];
        undecodedSaBwtCharCount[4] += undecodedSaBwtCharCount[3];

        // Stable bucket sort the decoded SA into olddIterationSaIndexList
        for (i=0; i<numOfSaProcessed; i++) {
            if (oddIterationSaIndexList[i].saIndex != 0) {
                bwtChar = oddIterationSaIndexList[i].textPositionIndex >> (BITS_IN_WORD - BIT_PER_CHAR);
                evenIterationSaIndexList[undecodedSaBwtCharCount[bwtChar]].saIndex = oddIterationSaIndexList[i].saIndex;
                evenIterationSaIndexList[undecodedSaBwtCharCount[bwtChar]].textPositionIndex = oddIterationSaIndexList[i].textPositionIndex & (ALL_ONE_MASK >> BIT_PER_CHAR);
                undecodedSaBwtCharCount[bwtChar]++;
            }
        }

        iteration++;
        numOfUndecodedSa = undecodedSaBwtCharCount[4];

    }

    return numOfUndecodedSa;

}

unsigned long long BWTDecompressText(const BWT *bwt, const unsigned long long endSaIndex, const unsigned long long length, const unsigned char *reverseCharMap,
                      unsigned char *decompressedText) {

    unsigned long long i, j;
    unsigned long long saIndex;

    saIndex = endSaIndex;

    for (i=length; i>0 && saIndex != 0; i--) {
        j = (saIndex > bwt->cumulativeFreq[1]) + (saIndex > bwt->cumulativeFreq[2])
                                               + (saIndex > bwt->cumulativeFreq[3]);
        decompressedText[i-1] = reverseCharMap[j];
        saIndex = BWTPsiMinusValue(bwt, saIndex);
    }

    return i;

}

unsigned long long BWTDecompressTextAsWordPacked(const BWT *bwt, const unsigned long long endSaIndex, const unsigned long long length, unsigned int *decompressedText) {

    unsigned long long i, j, k;
    unsigned long long saIndex;
    unsigned long long copyLeftShift, copyRightShift;
    unsigned long long copyLeft, copyRight;

    saIndex = endSaIndex;

    if (length / CHAR_PER_WORD * CHAR_PER_WORD == length) {

        for (i=length/CHAR_PER_WORD; i>0 && saIndex!=0; i--) {
            decompressedText[i-1] = 0;
            for (j=CHAR_PER_WORD; j>0 && saIndex!=0; j--) {
                k = (saIndex > bwt->cumulativeFreq[1]) + (saIndex > bwt->cumulativeFreq[2])
                                                       + (saIndex > bwt->cumulativeFreq[3]);
                decompressedText[i-1] |= (k << (BITS_IN_WORD - BIT_PER_CHAR)) >> ((j-1) * BIT_PER_CHAR);
                saIndex = BWTPsiMinusValue(bwt, saIndex);
            }
        }

    } else {

        for (i=length/CHAR_PER_WORD; i>0 && saIndex!=0; i--) {
            decompressedText[i] = 0;
            for (j=CHAR_PER_WORD; j>0 && saIndex!=0; j--) {
                k = (saIndex > bwt->cumulativeFreq[1]) + (saIndex > bwt->cumulativeFreq[2])
                                                       + (saIndex > bwt->cumulativeFreq[3]);
                decompressedText[i] |= (k << (BITS_IN_WORD - BIT_PER_CHAR)) >> ((j-1) * BIT_PER_CHAR);
                saIndex = BWTPsiMinusValue(bwt, saIndex);
            }
        }

        decompressedText[i] = 0;
        for (j=0; i*CHAR_PER_WORD+j<length && saIndex!=0; j++) {
            k = (saIndex > bwt->cumulativeFreq[1]) + (saIndex > bwt->cumulativeFreq[2])
                                                   + (saIndex > bwt->cumulativeFreq[3]);
            decompressedText[i] |= (k << (BITS_IN_WORD - BIT_PER_CHAR)) >> (j * BIT_PER_CHAR);
            saIndex = BWTPsiMinusValue(bwt, saIndex);
        }

        //Shift leftward
        copyLeftShift = BIT_PER_CHAR * length % CHAR_PER_WORD;
        copyRightShift = BITS_IN_WORD - copyLeftShift;
        for (i=0; i<length/CHAR_PER_WORD; i++) {
            copyLeft = decompressedText[i+1] >> copyLeftShift;
            copyRight = decompressedText[i] << copyRightShift;
            decompressedText[i] = copyLeft | copyRight;
        }
    }

    return i;

}

void BWTInitializeSaRetrievalStatistics(BWTSaRetrievalStatistics *bwtSaRetrievalStatistics) {

    bwtSaRetrievalStatistics->bwtSaRetrieved = 0;
    bwtSaRetrievalStatistics->saDiagonalLinked = 0;
    bwtSaRetrievalStatistics->cachedSaRetrieved = 0;
    bwtSaRetrievalStatistics->saDuplicated = 0;

}

void BWTAllocateDPStatistics(BWTDPStatistics *bwtDPStatistics) {

    bwtDPStatistics->totalNode = (LONG*) MMUnitAllocate((BWTDP_MAX_SUBSTRING_LENGTH+1) * sizeof(LONG));
    bwtDPStatistics->rejectedNode = (LONG*) MMUnitAllocate((BWTDP_MAX_SUBSTRING_LENGTH+1) * sizeof(LONG));
    bwtDPStatistics->totalDPCell = (LONG*) MMUnitAllocate((BWTDP_MAX_SUBSTRING_LENGTH+1) * sizeof(LONG));

}

void BWTInitializeDPStatistics(BWTDPStatistics *bwtDPStatistics) {

    int i;

    bwtDPStatistics->maxDepth = 0;
    bwtDPStatistics->maxDPCell = 0;
    bwtDPStatistics->maxDPMemoryInWord = 0;
    bwtDPStatistics->totalMaxDepth = 0;
    bwtDPStatistics->totalMaxDPCell = 0;
    bwtDPStatistics->totalMaxDPMemoryInWord = 0;
    bwtDPStatistics->acceptedPathDepth = 0;
    bwtDPStatistics->acceptedPath = 0;
    bwtDPStatistics->rejectedPathDepth = 0;
    bwtDPStatistics->rejectedPath = 0;
    for (i=0; i<=BWTDP_MAX_SUBSTRING_LENGTH; i++) {
        bwtDPStatistics->totalNode[i] = 0;
        bwtDPStatistics->rejectedNode[i] = 0;
        bwtDPStatistics->totalDPCell[i] = 0;
    }

}

void BWTFreeDPStatistics(BWTDPStatistics *bwtDPStatistics) {

    if (bwtDPStatistics->totalNode != NULL) {
        MMUnitFree(bwtDPStatistics->totalNode, (BWTDP_MAX_SUBSTRING_LENGTH+1) * sizeof(LONG));
    }
    if (bwtDPStatistics->rejectedNode != NULL) {
        MMUnitFree(bwtDPStatistics->rejectedNode, (BWTDP_MAX_SUBSTRING_LENGTH+1) * sizeof(LONG));
    }
    if (bwtDPStatistics->totalDPCell != NULL) {
        MMUnitFree(bwtDPStatistics->totalDPCell, (BWTDP_MAX_SUBSTRING_LENGTH+1) * sizeof(LONG));
    }

}

int SaIndexGroupStartSaIndexOrder(const void *saIndexGroup, const long long index1, const long long index2) {

    if (((SaIndexGroupNew*)saIndexGroup + index1)->startSaIndex != ((SaIndexGroupNew*)saIndexGroup + index2)->startSaIndex) {
        if (((SaIndexGroupNew*)saIndexGroup + index1)->startSaIndex > ((SaIndexGroupNew*)saIndexGroup + index2)->startSaIndex) {
            return 1;
        } else {
            return -1;
        }
    } else {
        return 0;
    }

}

int SaIndexGroupStartSaIndexLengthErrorOrder(const void *saIndexGroup, const long long index1, const long long index2) {

    if (((SaIndexGroupWithLengthError*)saIndexGroup + index1)->startSaIndex != ((SaIndexGroupWithLengthError*)saIndexGroup + index2)->startSaIndex) {
        if (((SaIndexGroupWithLengthError*)saIndexGroup + index1)->startSaIndex > ((SaIndexGroupWithLengthError*)saIndexGroup + index2)->startSaIndex) {
            return 1;
        } else {
            return -1;
        }
    } else {
        if (((SaIndexGroupWithLengthError*)saIndexGroup + index1)->length != ((SaIndexGroupWithLengthError*)saIndexGroup + index2)->length) {
            return ((SaIndexGroupWithLengthError*)saIndexGroup + index1)->length - ((SaIndexGroupWithLengthError*)saIndexGroup + index2)->length;
        } else {
            return ((SaIndexGroupWithLengthError*)saIndexGroup + index1)->error - ((SaIndexGroupWithLengthError*)saIndexGroup + index2)->error;
        }
    }

}

int HitListPosTextErrorLengthOrder(const void *hitList, const long long index1, const long long index2) {

    if (((HitListWithPosQueryLengthError*)hitList + index1)->posText != ((HitListWithPosQueryLengthError*)hitList + index2)->posText) {
        if (((HitListWithPosQueryLengthError*)hitList + index1)->posText > ((HitListWithPosQueryLengthError*)hitList + index2)->posText) {
            return 1;
        } else {
            return -1;
        }
    } else {
        if (((HitListWithPosQueryLengthError*)hitList + index1)->error != ((HitListWithPosQueryLengthError*)hitList + index2)->error) {
            return ((HitListWithPosQueryLengthError*)hitList + index1)->error - ((HitListWithPosQueryLengthError*)hitList + index2)->error;
        } else {
            return ((HitListWithPosQueryLengthError*)hitList + index1)->length - ((HitListWithPosQueryLengthError*)hitList + index2)->length;
        }
    }

}

int HitListPosText16BitOrder(const void *hitList, const long long index1, const long long index2) {

    return (((HitList*)hitList + index1)->posText >> 16) - (((HitList*)hitList + index2)->posText >> 16);

}

int HitListPosTextOrder(const void *hitList, const long long index1, const long long index2) {

    if (((HitList*)hitList + index1)->posText != ((HitList*)hitList + index2)->posText) {
        if (((HitList*)hitList + index1)->posText > ((HitList*)hitList + index2)->posText) {
            return 1;
        } else {
            return -1;
        }
    } else {
        return 0;
    }

}

int GappedHitListScorePosTextOrder(const void *gappedHitList, const long long index1, const long long index2) {

    if (((GappedHitList*)gappedHitList + index1)->score != ((GappedHitList*)gappedHitList + index2)->score) {
        return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
    } else {
        if (((GappedHitList*)gappedHitList + index1)->posText != ((GappedHitList*)gappedHitList + index2)->posText) {
            if (((GappedHitList*)gappedHitList + index1)->posText > ((GappedHitList*)gappedHitList + index2)->posText) {
                return 1;
            } else {
                return -1;
            }
        } else {
            return 0;
        }
    }

}

int GappedHitListDbSeqIndexScorePosTextOrder(const void *gappedHitList, const long long index1, const long long index2) {

    if (((GappedHitList*)gappedHitList + index1)->dbSeqIndex != ((GappedHitList*)gappedHitList + index2)->dbSeqIndex) {
        return ((GappedHitList*)gappedHitList + index1)->dbSeqIndex - ((GappedHitList*)gappedHitList + index2)->dbSeqIndex;
    } else {
        if (((GappedHitList*)gappedHitList + index1)->score != ((GappedHitList*)gappedHitList + index2)->score) {
            return ((GappedHitList*)gappedHitList + index2)->score - ((GappedHitList*)gappedHitList + index1)->score;
        } else {
            if (((GappedHitList*)gappedHitList + index1)->posText != ((GappedHitList*)gappedHitList + index2)->posText) {
                if (((GappedHitList*)gappedHitList + index1)->posText > ((GappedHitList*)gappedHitList + index2)->posText) {
                    return 1;
                } else {
                    return -1;
                }
            } else {
                return 0;
            }
        }
    }

}


int SaIndexGroupDPHitOrder1(const void *saIndexGroup, const long long index1, const long long index2) {

    if (((SaIndexGroupNew*)saIndexGroup + index1)->startSaIndex != ((SaIndexGroupNew*)saIndexGroup + index2)->startSaIndex) {
        if (((SaIndexGroupNew*)saIndexGroup + index1)->startSaIndex > ((SaIndexGroupNew*)saIndexGroup + index2)->startSaIndex) {
            return 1;
        } else {
            return -1;
        }
    } else {
        if (((SaIndexGroupNew*)saIndexGroup + index1)->numOfMatch != ((SaIndexGroupNew*)saIndexGroup + index2)->numOfMatch) {
            if (((SaIndexGroupNew*)saIndexGroup + index1)->numOfMatch < ((SaIndexGroupNew*)saIndexGroup + index2)->numOfMatch) {
                return 1;
            } else {
                return -1;
            }
        } else {
            if (((SaIndexGroupNew*)saIndexGroup + index1)->info != ((SaIndexGroupNew*)saIndexGroup + index2)->info) {
                if (((SaIndexGroupNew*)saIndexGroup + index1)->info < ((SaIndexGroupNew*)saIndexGroup + index2)->info) {
                    return 1;
                } else {
                    return -1;
                }
            } else {
                return 0;
            }
        }
    }

}

int SaIndexGroupDPHitOrder2(const void *saIndexGroup, const long long index1, const long long index2) {

    if (((SaIndexGroupNew*)saIndexGroup + index1)->posQuery != ((SaIndexGroupNew*)saIndexGroup + index2)->posQuery) {
        if (((SaIndexGroupNew*)saIndexGroup + index1)->posQuery < ((SaIndexGroupNew*)saIndexGroup + index2)->posQuery) {
            return 1;
        } else {
            return -1;
        }
    } else {
        if (((SaIndexGroupNew*)saIndexGroup + index1)->startSaIndex != ((SaIndexGroupNew*)saIndexGroup + index2)->startSaIndex) {
            if (((SaIndexGroupNew*)saIndexGroup + index1)->startSaIndex > ((SaIndexGroupNew*)saIndexGroup + index2)->startSaIndex) {
                return 1;
            } else {
                return -1;
            }
        } else {
            if (((SaIndexGroupNew*)saIndexGroup + index1)->numOfMatch != ((SaIndexGroupNew*)saIndexGroup + index2)->numOfMatch) {
                if (((SaIndexGroupNew*)saIndexGroup + index1)->numOfMatch < ((SaIndexGroupNew*)saIndexGroup + index2)->numOfMatch) {
                    return 1;
                } else {
                    return -1;
                }
            } else {
                if (((SaIndexGroupNew*)saIndexGroup + index1)->info != ((SaIndexGroupNew*)saIndexGroup + index2)->info) {
                    if (((SaIndexGroupNew*)saIndexGroup + index1)->info < ((SaIndexGroupNew*)saIndexGroup + index2)->info) {
                        return 1;
                    } else {
                        return -1;
                    }
                } else {
                    return 0;
                }
            }
        }
    }

}


