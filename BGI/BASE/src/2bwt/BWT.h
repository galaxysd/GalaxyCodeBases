/*

   BWT.h    BWT-Index

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
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef __BWT_H__
#define __BWT_H__

#include "TypeNLimit.h"
#include "MemManager.h"
#include "TextConverter.h"
#include "HSP.h"

#define BITS_PER_OCC_VALUE            16
#define OCC_VALUE_PER_WORD            2
#define OCC_VALUE_PER_LONG            4
#define OCC_INTERVAL                  256
#define WORD_BETWEEN_OCC              16
#define OCC_INTERVAL_MAJOR            65536

#define SORT_ALL                      0
#define SORT_16_BIT                   1
#define SORT_NONE                     2

#define BUCKET_BIT                    16
#define NUM_BUCKET                    65536

#define MAX_APPROX_MATCH_ERROR        7
#define MAX_ARPROX_MATCH_LENGTH       64

#define BWTDP_MAX_SUBSTRING_LENGTH    512

#define ESTIMATED_OCC_DIFF            32    // 128 / 4
#define MAX_OCC_DIFF                  128



typedef struct BWT {
    unsigned long long textLength;            // length of the text
    unsigned long long saInterval;            // interval between two SA values stored explicitly
    unsigned long long inverseSaInterval;        // interval between two inverse SA stored explicitly
    unsigned long long inverseSa0;            // SA-1[0]
    unsigned long long *cumulativeFreq;        // cumulative frequency
    unsigned int *bwtCode;                // BWT code
    unsigned int *occValue;                // Occurrence values stored explicitly
    unsigned long long *occValueMajor;        // Occurrence values stored explicitly
    unsigned long long *saValue;                // SA values stored explicitly
    unsigned long long *inverseSa;            // Inverse SA stored explicitly
    unsigned long long *cachedSaIndex;        // Cached SA index
    unsigned long long cachedSaIndexNumOfChar;    // Number of characters indexed in SA index range
    unsigned long long *saValueOnBoundary;    // Pre-calculated frequently referred data
    unsigned long long *decodeTable;            // For decoding BWT by table lookup
    unsigned int decodeTableGenerated;    // == TRUE if decode table is generated on load and will be freed
    unsigned long long bwtSizeInWord;            // Temporary variable to hold the memory allocated
    unsigned long long occSizeInWord;            // Temporary variable to hold the memory allocated
    unsigned long long occMajorSizeInWord;    // Temporary variable to hold the memory allocated
    unsigned long long saValueSizeInWord;        // Temporary variable to hold the memory allocated
    unsigned long long inverseSaSizeInWord;    // Temporary variable to hold the memory allocated
    unsigned long long cachedSaIndexSizeInWord;    // Temporary variable to hold the memory allocated
} BWT;

#define MAX_DIAGONAL_LEVEL 4                // Number of sub-pattern to keep for detecting diagonal hit

// Error information is stored as:
// 1. bitVector
//      After hamming distance match
// 2. count
//    After edit distance match
// 3. score
//    After the hits are processed with scoring functions

typedef struct SaIndexGroupNew {    // SA index range and information of a particular error arrangement of a matched sub-pattern
    unsigned long long startSaIndex;            // starting SA index
    unsigned int numOfMatch;                // number of match
    unsigned int posQuery;                // position in query; used for detecting diagonal hits
    unsigned int info;                    // extra hit information; to be copied to hitList.info
} SaIndexGroupNew;

typedef struct SaIndexGroupTemp {    // SA index range and information of a particular error arrangement of a matched sub-pattern
    unsigned long long startSaIndex1;            // starting SA index
    unsigned int numOfMatch1;            // number of match
    unsigned long long startSaIndex2;            // position in query; used for detecting diagonal hits
    unsigned int numOfMatch2;            // extra hit information; to be copied to hitList.info
} SaIndexGroupTemp;

typedef struct SaIndexGroupOld {    // SA index range and information of a particular error arrangement of a matched sub-pattern
    unsigned long long startSaIndex;            // starting SA index
    unsigned int numOfMatch;                // number of match
    unsigned int info;                    // extra hit information; to be copied to hitList.info
} SaIndexGroupOld;

typedef struct SaIndexGroup {    // SA index range and information of a particular error arrangement of a matched sub-pattern
    unsigned long long startSaIndex;            // starting SA index
    unsigned int numOfMatch;            // number of match
    unsigned int info;                    // extra hit information
} SaIndexGroup;

typedef struct SaIndexGroupWithErrorBitVector {    // SA index range and information of a particular error arrangement of a matched sub-pattern
    unsigned long long startSaIndex;            // starting SA index
    unsigned int numOfMatch;            // number of match
    unsigned int errorBitVector;            // error bit vector
} SaIndexGroupWithErrorBitVector;

typedef struct SaIndexGroupWithLengthError {    // SA index range and information of a particular error arrangement of a matched sub-pattern
    unsigned long long startSaIndex;            // starting SA index
    unsigned int numOfMatch;            // number of match
    unsigned posQuery : 16;        // position in query
    unsigned length   : 8;        // length of hit
    unsigned error    : 8;        // error in hit
} SaIndexGroupWithLengthError;

typedef struct SaIndexGroupProcessed {    // Alternative usage of SaIndexGroup - once processed, error bit vector is replaced by index to text position
    unsigned long long startSaIndex;            // starting SA index
    unsigned int numOfMatch;            // number of match
    unsigned long long textPositionIndex;        // storing the pointer to text position
} SaIndexGroupProcessed;

typedef struct DupSaIndexGroup {    // Alternative usage of SaIndexGroup - the group duplicates another group
    unsigned long long lastDupSaIndexGroupIndex;    // index to last duplicated group
    unsigned long long saIndexGroupIndex;            // index to the first SA into group among the duplicates
    unsigned long long textPositionIndex;            // storing the pointer to text position
} DupSaIndexGroup;

typedef struct SaIndexGroupHash {    // Hash table for checking duplicate SA index group
    unsigned long long startSaIndex;
    unsigned long long saIndexGroupIndex;
} SaIndexGroupHash;

typedef struct BWTSaRetrievalStatistics {
    unsigned long long bwtSaRetrieved;
    unsigned long long saDiagonalLinked;
    unsigned long long saDuplicated;
    unsigned long long cachedSaRetrieved;
} BWTSaRetrievalStatistics;

typedef struct BWTDPStatistics {
    int maxDepth;
    int maxDPCell;
    int maxDPMemoryInWord;
    int totalMaxDepth;
    int totalMaxDPCell;
    int totalMaxDPMemoryInWord;
    LONG acceptedPathDepth;
    LONG acceptedPath;
    LONG rejectedPathDepth;
    LONG rejectedPath;
    LONG* __restrict totalNode;
    LONG* __restrict rejectedNode;
    LONG* __restrict totalDPCell;
} BWTDPStatistics;

typedef struct SaIndexList {
    unsigned long long saIndex;
    unsigned long long textPositionIndex;
} SaIndexList;

typedef struct HitCombination {
    int numOfCombination;
    int maxError;
    int keyLength;
    int skipTableWidth;
    int *errorPos;
    int *skip;
    int *skipErrorIndex;
} HitCombination;

typedef struct DPText {
    int charBeingProcessed;
    int dpCellIndex;
    int numOfDpCellSegment;
    unsigned long long dummy1;    // Must not be removed; so that saIndexLeft and saIndexRight are aligned to 16 byte boundary
    unsigned long long saIndexLeft[ALPHABET_SIZE];
    unsigned long long saIndexRight[ALPHABET_SIZE];
} DPText;

typedef struct DPScanDepth {
    unsigned P                :    31;
    unsigned withAmbiguity    :    1;
} DPScanDepth;


// Load / unload functions
BWT *BWTCreate(MMPool *mmPool, const unsigned long long textLength, unsigned long long *decodeTable);
BWT *BWTLoad(MMPool *mmPool, const char *bwtCodeFileName, const char *occValueFileName, 
             const char *saValueFileName, const char *inverseSaFileName, const char *saIndexRangeFileName,
             unsigned long long *decodeTable);
void BWTFree(MMPool *mmPool, BWT *bwt);
void BWTPrintMemoryUsage(const BWT *bwt, FILE *output, const unsigned long long packedDNASize);

// Precalculate frequenctly accessed data
void BWTGenerateSaValueOnBoundary(MMPool *mmPool, BWT *bwt);

// Core functions
// The following must be customized for differenet compression schemes ***
unsigned long long BWTDecode(const BWT *bwt, const unsigned long long index1, const unsigned long long index2, const unsigned int character);
void BWTDecodeAll(const BWT *bwt, const unsigned long long index1, const unsigned long long index2, unsigned long long* __restrict occValue);
unsigned long long BWTOccValue(const BWT *bwt, unsigned long long index, const unsigned int character);
void BWTOccValueTwoIndex(const BWT *bwt, unsigned long long index1, unsigned long long index2, const unsigned int character, unsigned long long* __restrict occValue);
void BWTAllOccValue(const BWT *bwt, unsigned long long index, unsigned long long* __restrict occValue);
void BWTAllOccValueTwoIndex(const BWT *bwt, unsigned long long index1, unsigned long long index2, unsigned long long* __restrict occValue1, unsigned long long* __restrict occValue2);
unsigned long long BWTOccValueOnSpot(const BWT *bwt, unsigned long long index, unsigned long long* __restrict character);
unsigned long long BWTSearchOccValue(const BWT *bwt, const unsigned int character, const unsigned long long searchOccValue);


// Utility functions for no compression only
unsigned long long BWTResidentSizeInWord(const unsigned long long numChar);
unsigned long long BWTFileSizeInWord(const unsigned long long numChar);
void BWTClearTrailingBwtCode(BWT *bwt);

// These are generic to different compression schemes (and generic to no compression as well)
unsigned long long BWTPsiMinusValue(const BWT *bwt, const unsigned long long index);
unsigned long long BWTPsiPlusValue(const BWT *bwt, const unsigned long long index);
unsigned long long BWTSaValue(const BWT *bwt, unsigned long long index);
unsigned long long BWTInverseSa(const BWT *bwt, unsigned long long saValue);
unsigned long long BWTOccIntervalMajor(const unsigned long long occInterval);
unsigned long long BWTOccValueMinorSizeInWord(const unsigned long long numChar);
unsigned long long BWTOccValueMajorSizeInWord(const unsigned long long numChar);

// Search functions
// packedText should be allocated with at least 1 Word buffer initialized to zero
int BWTForwardSearch(const unsigned int *packedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText);
int BWTForwardSearchSaIndex(const unsigned int *packedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText, 
                     unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight);
int BWTSaBinarySearch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int *packedText, 
                      unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight, unsigned int *tempKey);    // tempKey = buffer large enough to hold packed key
int BWTBackwardSearch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, 
                      unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight);
int BWTBackwardSearchCheckWithText(const unsigned char *convertedKey, const unsigned int *packedKey, const unsigned int keyLength,
                              const BWT *bwt, const unsigned int *packedText, const unsigned int textCheckingCostFactor,
                              const unsigned int maxnumOfTextPosition, HitList* __restrict hitList,
                              unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight);

// Approximate match functions - brute force deep first search by backward search is used
unsigned long long BWTHammingDistMaxSaIndexGroup(const unsigned int keyLength, const unsigned int maxError);
unsigned long long BWTHammingDistCountOcc(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError);
unsigned long long BWTHammingDistMatch(const BWT *bwt, const unsigned char *convertedKey, const HitCombination *hitCombination,
                        const unsigned long long *cachedSaIndex, const unsigned long long cachedSaIndexNumOfChar,
                        SaIndexGroupNew* __restrict saIndexGroup, const unsigned int maxSaIndexGroup);
int BWTHammingDistMatchOld(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
                           SaIndexGroupNew* __restrict saIndexGroup, const unsigned long long maxSaIndexGroup,
                           const unsigned long long posQuery, const unsigned long long info);
unsigned long long BWTEditDistMaxSaIndexGroup(const unsigned int keyLength, const unsigned int maxError);
// Does not insert characters on pattern boundary
unsigned long long BWTEditDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
                     SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned int maxSaIndexGroup);
unsigned long long BWTEditDistMatchOld(const unsigned char *convertedKey, const unsigned int keyLength, const BWT *bwt, const unsigned int maxError,
                     SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned int maxSaIndexGroup);


unsigned long long BWTEliminateDupSaIndexGroup(SaIndexGroupWithLengthError* __restrict saIndexGroup, const unsigned long long numOfSaGroup);

unsigned long long BWTSubPatternHammingDistCountOcc(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
                                     const unsigned int maxError, const unsigned int skip);
int BWTSubPatternHammingDistSaIndex(const BWT *bwt, const unsigned char *convertedKey, const int keyLength, const int skip, 
                                    const HitCombination *hitCombination, 
                                    const unsigned long long *cachedSaIndex, const unsigned long long cachedSaIndexNumOfChar,
                                    SaIndexGroupNew* __restrict saIndexGroup, const int maxnumOfSaIndexGroup,
                                    int* __restrict firstSaIndexGroupForSubPattern);
int BWTSubPatternHammingDistSaIndexOld(const unsigned char *convertedKey, const int keyLength, const int subPatternLength, const BWT *bwt, 
                                    const int maxError, const int skip, const int lengthProcessed, int* __restrict lengthInCurrentRound,
                                    SaIndexGroupNew* __restrict saIndexGroup, const int maxnumOfSaIndexGroup);
int BWTDPHit(const BWT *bwt, SaIndexGroupNew* __restrict saIndexGroup, const long long numOfSaIndexGroup, 
             const int firstSaIndexGrouptoProcess, int* __restrict saIndexGroupProcessed,
             const int discardDiagonalHit,
             char* workingMemory, const int workingMemorySize,
             BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics);
//unsigned int BWTSubPatternHammingDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
//                                  const unsigned int maxError, const unsigned int skip, const unsigned int matchBitVector,
//                                  const SaIndexRange *saIndexRange, const unsigned int saIndexRangeNumOfChar,
//                                  char *workingMemory, const unsigned int workingMemorySize, const unsigned int sortOption,
//                                  BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics);
unsigned long long BWTSubPatternEditDistMatch(const unsigned char *convertedKey, const unsigned int keyLength, const unsigned int subPatternLength, const BWT *bwt, 
                                  const unsigned int maxError, const unsigned int skip, const unsigned int maxnumOfHit, 
                                  HitListWithPosQueryLengthError* __restrict hitList, BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics, 
                                  const unsigned int eliminateDuplicateStartingPos);
int BWTGappedDPDBvsQuery(BWT *bwt, const unsigned char *convertedKey, const int queryPatternLength, 
                         char* __restrict workingMemory, const int workingMemorySize, int* __restrict totalNumOfQueryPos,
                         const int matchScore, const int mismatchScore,
                         const int gapOpenScore, const int gapExtendScore,
                         const int cutoffScore,
                         BWTDPStatistics* __restrict bwtDPStatistics,
                         const int printProgressDepth);

// Text retrieval functions
// Position in text will be placed at the first word of hitListSizeInWord

// startSaIndex + resultInfo must be sorted in increasing order; there must be no overlapping groups except that one group can completely enclose another
unsigned long long BWTTextPosition(const BWT *bwt, const SaIndexGroupNew *saIndexGroup, const unsigned long long numOfSaIndexGroups, 
                    HitList* __restrict hitList, 
                    SaIndexList* __restrict tempSaIndexList1, SaIndexList* __restrict tempSaIndexList2, 
                    BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics, const unsigned long long eliminateDuplicateStartingPos);
                    
unsigned long long BWTDecodeTextPosition(const BWT *bwt, SaIndexList* __restrict evenIterationSaIndexList, SaIndexList* __restrict oddIterationSaIndexList,
                          const unsigned long long numOfSaIndex, const unsigned long long numOfIterationProcessed, const unsigned long long maxnumOfIteration, 
                          HitList* __restrict hitList, BWTSaRetrievalStatistics* __restrict bwtSaRetrievalStatistics);


unsigned long long BWTDecompressText(const BWT *bwt, const unsigned long long endSaIndex, const unsigned long long length, const unsigned char *reverseCharMap,
                      unsigned char *decompressedText);
unsigned long long BWTDecompressTextAsWordPacked(const BWT *bwt, const unsigned long long endSaIndex, const unsigned long long length, unsigned long long *decompressedText);

void BWTPrintDPStatistics(FILE * outFile, const BWTDPStatistics* bwtDPStatistics);
void BWTInitializeSaRetrievalStatistics(BWTSaRetrievalStatistics *bwtSaRetrievalStatistics);
void BWTAllocateDPStatistics(BWTDPStatistics *bwtDPStatistics);
void BWTInitializeDPStatistics(BWTDPStatistics *bwtDPStatistics);
void BWTFreeDPStatistics(BWTDPStatistics *bwtDPStatistics);

// QSort comparison functions
int SaIndexGroupStartSaIndexOrder(const void *saIndexGroup, const long long index1, const long long index2);
int SaIndexGroupStartSaIndexLengthErrorOrder(const void *saIndexGroup, const long long index1, const long long index2);
int HitListPosTextErrorLengthOrder(const void *hitList, const long long index1, const long long index2);
int HitListPosText16BitOrder(const void *hitList, const long long index1, const long long index2);
int HitListPosTextOrder(const void *hitList, const long long index1, const long long index2);
int GappedHitListScorePosTextOrder(const void *gappedHitList, const long long index1, const long long index2);
int GappedHitListDbSeqIndexScorePosTextOrder(const void *gappedHitList, const long long index1, const long long index2);


#endif
