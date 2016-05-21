/*

#    Copyright (C) 2015, The University of Hong Kong.
#
#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 3
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  03110-1301, USA.

   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Packaging 2BWT library as a separate product.
            Enhancing the 2BWT library with interface functions
            for basic BWT search.

*/

#include "2BWT-Interface.h"

Idx2BWT * BWTLoad2BWT(const char * indexFilePrefix, const char * saFileNameExtension) {

    Idx2BWT * idx2BWT = (Idx2BWT*) malloc(sizeof(Idx2BWT));
    BWT * bwt;
    BWT * rev_bwt;
    HSP * hsp;
    
    char bwtFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char bwtOccFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char saFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char rev_bwtFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char rev_bwtOccFilename[MAX_INDEX_FILENAME_LENGTH];
    char packedDnaFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char annotationFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char ambiguityFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char translateFilename[MAX_INDEX_FILENAME_LENGTH]; 
    
    strcpy(bwtFilename,indexFilePrefix);
    strcpy(bwtOccFilename,indexFilePrefix);
    strcpy(saFilename,indexFilePrefix);
    strcpy(rev_bwtFilename,indexFilePrefix);
    strcpy(rev_bwtOccFilename,indexFilePrefix);
    strcpy(packedDnaFilename,indexFilePrefix);
    strcpy(annotationFilename,indexFilePrefix);
    strcpy(ambiguityFilename,indexFilePrefix);
    strcpy(translateFilename,indexFilePrefix);
    
    strcat(bwtFilename,".bwt");
    strcat(bwtOccFilename,".fmv");
    strcat(saFilename,saFileNameExtension);
    strcat(rev_bwtFilename,".rev.bwt");
    strcat(rev_bwtOccFilename,".rev.fmv");
    strcat(packedDnaFilename,".pac");
    strcat(annotationFilename,".ann");
    strcat(ambiguityFilename,".amb");
    strcat(translateFilename,".tra");
    
    MMMasterInitialize(3, 0, FALSE, NULL);
    MMPool * mmPool = MMPoolCreate(2097152);
    
    bwt = BWTLoad(mmPool, bwtFilename, bwtOccFilename, saFilename, NULL, NULL, NULL);
    rev_bwt = BWTLoad(mmPool, rev_bwtFilename, rev_bwtOccFilename, NULL, NULL, NULL, NULL);
    hsp = HSPLoad(mmPool, packedDnaFilename, annotationFilename, ambiguityFilename,translateFilename, 1);
    
    HSPFillCharMap(idx2BWT->charMap);
    HSPFillComplementMap(idx2BWT->complementMap);
     
    idx2BWT->bwt = bwt;
    idx2BWT->rev_bwt = rev_bwt;
    idx2BWT->hsp = hsp;
    idx2BWT->mmPool = mmPool;

    idx2BWT->numReads = 0;
    idx2BWT->readIDtable = NULL;
    
    return idx2BWT;
}

Idx2BWT *BWTLoad2BWTLite(const char * indexFilePrefix) {
    Idx2BWT * idx2BWT = (Idx2BWT*) malloc(sizeof(Idx2BWT));
    BWT * bwt;
    BWT * rev_bwt;

    char bwtFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char bwtOccFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char rev_bwtFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char rev_bwtOccFilename[MAX_INDEX_FILENAME_LENGTH];
    char bwtReadIDtableFilename[MAX_INDEX_FILENAME_LENGTH];

    strcpy(bwtFilename,indexFilePrefix);
    strcpy(bwtOccFilename,indexFilePrefix);
    strcpy(rev_bwtFilename,indexFilePrefix);
    strcpy(rev_bwtOccFilename,indexFilePrefix);
    strcpy(bwtReadIDtableFilename, indexFilePrefix);

    strcat(bwtFilename,".bwt");
    strcat(bwtOccFilename,".fmv");
    strcat(rev_bwtFilename,".rev.bwt");
    strcat(rev_bwtOccFilename,".rev.fmv");
    strcat(bwtReadIDtableFilename, ".ridt");

    static int is_init = 0;
    if (!is_init)
        MMMasterInitialize(6, 0, FALSE, NULL);
    is_init = 1;
    MMPool * mmPool = MMPoolCreate(2097152);

    bwt = BWTLoad(mmPool, bwtFilename, bwtOccFilename, NULL, NULL, NULL, NULL);
    rev_bwt = BWTLoad(mmPool, rev_bwtFilename, rev_bwtOccFilename, NULL, NULL, NULL, NULL);

    HSPFillCharMap(idx2BWT->charMap);
    HSPFillComplementMap(idx2BWT->complementMap);

    idx2BWT->bwt = bwt;
    idx2BWT->rev_bwt = rev_bwt;
    idx2BWT->hsp = NULL;
    idx2BWT->mmPool = mmPool;

    // read read id table
    FILE *readIDtableFile = fopen64(bwtReadIDtableFilename, "rb");
    if (readIDtableFile == NULL) {
        fprintf(stderr, "%s: cannot open read_id_table!\n", __func__);
        exit(1);
    }
    fread(&(idx2BWT->numReads), sizeof(idx2BWT->numReads), 1, readIDtableFile);
    if (bwt->cumulativeFreq[4] % idx2BWT->numReads != 0) {
        fprintf(stderr, "Number of reads and read length do not match!\n");
        exit(1);
    }
    idx2BWT->readLength = bwt->cumulativeFreq[4] / idx2BWT->numReads;
    fprintf(stderr, "Read length %d\n", idx2BWT->readLength);
    idx2BWT->readIDtable = (unsigned int*) malloc(sizeof(unsigned int) * idx2BWT->numReads);
    if (idx2BWT->readIDtable == NULL) {
        fprintf(stderr, "%s:%s:%d: memory alloc failed!\n", __FILE__, __func__, __LINE__);
        exit(1);
    }
    fread(idx2BWT->readIDtable, sizeof(unsigned int), idx2BWT->numReads, readIDtableFile);

    fclose(readIDtableFile);
    
    return idx2BWT;
}

void BWTFree2BWT(Idx2BWT * idx2BWT) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    HSP * hsp = idx2BWT->hsp;
    MMPool * mmPool = idx2BWT->mmPool;
    
    if (hsp != NULL) {
        HSPFree(mmPool, hsp, 1);
    }
    BWTFree(mmPool, bwt);
    BWTFree(mmPool, rev_bwt);
    MMPoolFree(mmPool);

    if (idx2BWT->readIDtable != NULL) {
        free(idx2BWT->readIDtable);
    }
    
    free(idx2BWT);
}

int BWTGetQuality(Idx2BWT * idx2BWT, unsigned long long readID, int offset) {

    BWT * bwt = idx2BWT->bwt;
    unsigned long long absOffset = readID * idx2BWT->readLength + offset;
    return (bwt->bwtCode[absOffset / CHAR_PER_WORD] >> (absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1))) & 1;
}

int BWTCheckUsed(Idx2BWT * idx2BWT, unsigned long long readID) {

    BWT * bwt = idx2BWT->bwt;
    unsigned long long absOffset = bwt->cumulativeFreq[4] + readID;
    return (bwt->bwtCode[absOffset / CHAR_PER_WORD] >> (absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1))) & 1;
}

void BWTSetUsed(Idx2BWT * idx2BWT, unsigned long long readID, int used) {
    BWT * bwt = idx2BWT->bwt;
    unsigned long long absOffset = bwt->cumulativeFreq[4] + readID;
    unsigned int *wordToUpdate = bwt->bwtCode + absOffset / CHAR_PER_WORD;
    int offsetInWord = absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1);

    if (used) {
        *wordToUpdate |= 1 << offsetInWord;
    } else {
        *wordToUpdate &= ~(1 << offsetInWord);
    }
}

void BWTCheckFlag(Idx2BWT * idx2BWT, unsigned long long readID, int *used, int *cons) {

    BWT * bwt = idx2BWT->bwt;
    unsigned long long absOffset = bwt->cumulativeFreq[4] + readID;
    *used = (bwt->bwtCode[absOffset / CHAR_PER_WORD] >> (absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1))) & 1;
    
    bwt = idx2BWT->rev_bwt;
    absOffset = bwt->cumulativeFreq[4] + readID;
    *cons = (bwt->bwtCode[absOffset / CHAR_PER_WORD] >> (absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1))) & 1;
}

void BWTSetFlag(Idx2BWT * idx2BWT, unsigned long long readID, int used, int cons) {
    BWT * bwt = idx2BWT->bwt;
    unsigned long long absOffset = bwt->cumulativeFreq[4] + readID;
    unsigned int *wordToUpdate = bwt->bwtCode + absOffset / CHAR_PER_WORD;
    int offsetInWord = absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1);

    if (used) {
        *wordToUpdate |= 1 << offsetInWord;
    } else {
        *wordToUpdate &= ~(1 << offsetInWord);
    }

    bwt = idx2BWT->rev_bwt;
    absOffset = bwt->cumulativeFreq[4] + readID;
    wordToUpdate = bwt->bwtCode + absOffset / CHAR_PER_WORD;
    offsetInWord = absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1);

    if(cons)
        *wordToUpdate |= 1 << offsetInWord;
    else
        *wordToUpdate &= ~(1 << offsetInWord);

}
void BWTConvertPattern(Idx2BWT * idx2BWT, const char * patternSource, int patternLength, unsigned char * patternDestination) {

    int i;
    for (i=0;i<patternLength;i++) {
        patternDestination[i] = idx2BWT->charMap[patternSource[i]];
    }
    patternDestination[i]='\0';
}


void BWTSARangeInitial(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned long long *saIndexLeft, unsigned long long *saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    (*saIndexLeft) = bwt->cumulativeFreq[c]; //+1;
    (*saIndexRight) = bwt->cumulativeFreq[c+1]-1;

}


void BWTSARangeBackward(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned long long *saIndexLeft, unsigned long long *saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;

    
    unsigned long long l = (*saIndexLeft);
    unsigned long long r = (*saIndexRight);
    (*saIndexLeft)  = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c);// + 1;
    (*saIndexRight) = bwt->cumulativeFreq[c] + BWTOccValue(bwt, r + 1, c) - 1;
}

void BWTSARangeBackward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned long long *saIndexLeft, unsigned long long *saIndexRight,
                        unsigned long long *rev_saIndexLeft, unsigned long long *rev_saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned long long oL[ALPHABET_SIZE];
    unsigned long long oR[ALPHABET_SIZE];
    unsigned long long oCount[ALPHABET_SIZE];
    unsigned long long l = (*saIndexLeft);
    unsigned long long r = (*saIndexRight);
    unsigned long long rev_l = (*rev_saIndexLeft);
    unsigned long long rev_r = (*rev_saIndexRight);
    int k;
    
    BWTAllOccValue(bwt,l,oL);
    BWTAllOccValue(bwt,r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }

    l = bwt->cumulativeFreq[c] + oL[c] + 1;
    r = bwt->cumulativeFreq[c] + oR[c];
    rev_r = rev_r - oCount[c];
    rev_l = rev_r - (r-l);

    l--;
    r--;

    (*saIndexLeft) = l;
    (*saIndexRight) = r;
    (*rev_saIndexLeft) = rev_l;
    (*rev_saIndexRight) = rev_r;
    
}
void BWTSARangeForward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned long long *saIndexLeft, unsigned long long *saIndexRight,
                        unsigned long long *rev_saIndexLeft, unsigned long long *rev_saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned long long oL[ALPHABET_SIZE];
    unsigned long long oR[ALPHABET_SIZE];
    unsigned long long oCount[ALPHABET_SIZE];
    unsigned long long l = (*saIndexLeft);
    unsigned long long r = (*saIndexRight);
    unsigned long long rev_l = (*rev_saIndexLeft);
    unsigned long long rev_r = (*rev_saIndexRight);
    int k;
    
    BWTAllOccValue(rev_bwt,rev_l,oL);
    BWTAllOccValue(rev_bwt,rev_r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }
    
    rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
    rev_r = bwt->cumulativeFreq[c] + oR[c];
    r = r - oCount[c];
    l = r - (rev_r-rev_l);

    rev_r--;
    rev_l--;
    
    (*saIndexLeft) = l;
    (*saIndexRight) = r;
    (*rev_saIndexLeft) = rev_l;
    (*rev_saIndexRight) = rev_r;
}

void BWTAllSARangesBackward(Idx2BWT * idx2BWT, 
                        const unsigned long long saIndexLeft, const unsigned long long saIndexRight,
                        unsigned long long *resultSaIndexesLeft, unsigned long long *resultSaIndexesRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned long long oL[ALPHABET_SIZE];
    unsigned long long oR[ALPHABET_SIZE];
    unsigned long long oCount[ALPHABET_SIZE];
    unsigned long long l = saIndexLeft;
    unsigned long long r = saIndexRight;
    int k;
    
    BWTAllOccValue(bwt,l,oL);
    BWTAllOccValue(bwt,r + 1,oR);
    
    for (k=0;k<ALPHABET_SIZE;k++) {
        resultSaIndexesLeft[k]  = bwt->cumulativeFreq[k] + oL[k] + 1 - 1;
        resultSaIndexesRight[k] = bwt->cumulativeFreq[k] + oR[k] - 1;
    }

}


void BWTAllSARangesBackward_Bidirection(Idx2BWT * idx2BWT, 
                        const unsigned long long saIndexLeft, const unsigned long long saIndexRight,
                        const unsigned long long rev_saIndexLeft, const unsigned long long rev_saIndexRight,
                        unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight,
                        unsigned long long *rev_resultSaIndexLeft, unsigned long long *rev_resultSaIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned long long oL[ALPHABET_SIZE];
    unsigned long long oR[ALPHABET_SIZE];
    unsigned long long oCount[ALPHABET_SIZE];
    unsigned long long l = saIndexLeft;
    unsigned long long r = saIndexRight;
    unsigned long long rev_l = rev_saIndexLeft;
    unsigned long long rev_r = rev_saIndexRight;
    int k;
    
    BWTAllOccValue(bwt,l,oL);
    BWTAllOccValue(bwt,r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }
    
    for (k=0;k<ALPHABET_SIZE;k++) {

        resultSaIndexLeft[k] = bwt->cumulativeFreq[k] + oL[k] + 1 - 1;
        resultSaIndexRight[k] = bwt->cumulativeFreq[k] + oR[k] - 1;
        rev_resultSaIndexRight[k] = rev_r - oCount[k] - 1;
        rev_resultSaIndexLeft[k] = rev_resultSaIndexRight[k] - 
        (resultSaIndexRight[k]-resultSaIndexLeft[k]) - 1;

    }

}

void BWTAllSARangesForward_Bidirection(Idx2BWT * idx2BWT, 
                        const unsigned long long saIndexLeft, const unsigned long long saIndexRight,
                        const unsigned long long rev_saIndexLeft, const unsigned long long rev_saIndexRight,
                        unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight,
                        unsigned long long *rev_resultSaIndexLeft, unsigned long long *rev_resultSaIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned long long oL[ALPHABET_SIZE];
    unsigned long long oR[ALPHABET_SIZE];
    unsigned long long oCount[ALPHABET_SIZE];
    unsigned long long l = saIndexLeft;
    unsigned long long r = saIndexRight;
    unsigned long long rev_l = rev_saIndexLeft;
    unsigned long long rev_r = rev_saIndexRight;
    int k;
    
    BWTAllOccValue(rev_bwt,rev_l,oL);
    BWTAllOccValue(rev_bwt,rev_r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }
    
    for (k=0;k<ALPHABET_SIZE;k++) {

        rev_resultSaIndexLeft[k] = bwt->cumulativeFreq[k] + oL[k] + 1 - 1;
        rev_resultSaIndexRight[k] = bwt->cumulativeFreq[k] + oR[k] - 1;
        resultSaIndexRight[k] = r - oCount[k] - 1;
        resultSaIndexLeft[k] = resultSaIndexRight[k] - 
        (rev_resultSaIndexRight[k]-rev_resultSaIndexLeft[k]) - 1;

    }
}

void BWTRetrievePositionFromSAIndex(Idx2BWT * idx2BWT, 
                                    unsigned long long saIndex, 
                                    unsigned int * sequenceId, unsigned long long * offset) {
    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    HSP * hsp = idx2BWT->hsp;
    unsigned short * ambiguityMap = hsp->ambiguityMap;
    Translate * translate = hsp->translate;
    
    unsigned long long ambPosition = BWTSaValue(bwt,saIndex);
    unsigned long long approxIndex = ambPosition>>GRID_SAMPLING_FACTOR_2_POWER;
    unsigned long long approxValue = ambiguityMap[approxIndex];
    while (translate[approxValue].startPos>ambPosition) {
        approxValue--;
    }
    ambPosition-=translate[approxValue].correction;
    
    (*sequenceId) = translate[approxValue].chrID;
    (*offset) = ambPosition;
}

