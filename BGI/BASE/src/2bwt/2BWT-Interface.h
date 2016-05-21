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

#ifndef __2BWT_INTERFACE_H__
#define __2BWT_INTERFACE_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "BWT.h"
#include "HSP.h"
#include "MemManager.h"

#define MAX_INDEX_FILENAME_LENGTH 1024
#define POOLSIZE                  2097152

typedef struct Idx2BWT {
    MMPool * mmPool;
    BWT * bwt;
    BWT * rev_bwt;
    HSP * hsp;
    unsigned char charMap[256];
    unsigned char complementMap[256];

    // for short read index
    unsigned int numReads;
    unsigned int readLength;
    unsigned int *readIDtable;
} Idx2BWT;



//=============================================================================
//  BWTLoad2BWT
//  This function loads up the 2BWT index from memory and return the placeholder
//  Input :     indexFilePrefix - the file name prefix of the 2BWT index files
//                                e.g. "ncbi.genome.fa.index"
//              saFileNameExtension - the file name extension of suffix array in
//                                    the 2BWT index. This string should correspond
//                                    to the last part of the SaValueFileName value 
//                                    in 2bwt-builder.ini when the index was built.
//                                    e.g. ".sa"
//  Output :    function returns the placeholder of 2BWT index
//=============================================================================
Idx2BWT * BWTLoad2BWT(const char * indexFilePrefix, const char * saFileNameExtension);
Idx2BWT * BWTLoad2BWTLite(const char * indexFilePrefix);

//=============================================================================
//  BWTFree2BWT
//  This function frees the 2BWT index from memory and destory the placeholder
//  Input :     idx2BWT - placeholder of 2BWT index
//  Output :    none
//=============================================================================
void BWTFree2BWT(Idx2BWT * idx2BWT);

//=============================================================================
//  BWTGetQuality
//  For short read BWT only
//  Return the 0/1 quality of a base specified by readID and offset(0 based)
//  Input :     idx2BWT, readID, offset
//  Output :    quality (0 or 1)
//=============================================================================
int BWTGetQuality(Idx2BWT * idx2BWT, unsigned long long readID, int offset);

//=============================================================================
//  BWTCheckUsed
//  For short read BWT only
//  Check whether a read has been used
//  Input :     idx2BWT, readID
//  Output :    0 for not unused and 1 for used
//=============================================================================
int BWTCheckUsed(Idx2BWT * idx2BWT, unsigned long long readID);

//=============================================================================
//  BWTSetUsed
//  For short read BWT only
//  Check whether a read has been used
//  Input :     idx2BWT, readID, used (1 for used and 0 for unused)
//  Output :    none
//=============================================================================
void BWTSetUsed(Idx2BWT * idx2BWT, unsigned long long readID, int used);


//=============================================================================
//  BWTCheckFlag
//  For short read BWT only
//  Check whether a read has been used
//  Input :     idx2BWT, readID
//  Output :    0 for not unused and 1 for used
//=============================================================================
void BWTCheckFlag(Idx2BWT * idx2BWT, unsigned long long readID, int *used, int *cons);
//=============================================================================
//  BWTSetFlag
//  For short read BWT only
//  Check whether a read has been used
//  Input :     idx2BWT, readID, used (1 for used and 0 for unused)
//  Output :    none
//=============================================================================
void BWTSetFlag(Idx2BWT * idx2BWT, unsigned long long readID, int used, int cons);

//=============================================================================
//  BWTConvertPattern
//  This function converts a human readable nucleotide sequence (e.g. ACCAACAT) 
//  into a coding scheme recognised by 2BWT index (e.g. 01100103)
//  Input :     idx2BWT - placeholder of 2BWT index
//              patternSource - the human readable nucleotide sequence
//              patternLength - the length of the pattern, obviously
//  Output :    patternDestination - the converted pattern
//=============================================================================
void BWTConvertPattern(Idx2BWT * idx2BWT, const char * patternSource, int patternLength, unsigned char * patternDestination);

//=============================================================================
//  BWTSARangeInitial
//  This function performs the 2BWT search of the first character 
//  of the pattern and returns the resulting SA range to the input parameter 
//  resultSaIndexLeft and resultSaIndexRight.
//  Input :     idx2BWT - placeholder of 2BWT index
//              c       - the character to search for
//  Output :    resultSaIndexLeft  - lower bound of the resulting SA range
//              resultSaIndexRight - upper bound of the resulting SA range
//=============================================================================
void BWTSARangeInitial(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight);
     
//=============================================================================
//  BWTSARangeBackward
//  This function performs the 2BWT backward search for 1 character
//  and returns the resulting SA range to the input parameter 
//  resultSaIndexLeft and resultSaIndexRight.
//  Input :     idx2BWT - placeholder of 2BWT index
//              c   - the character to search for
//  Output :    saIndexLeft  - lower bound of the resulting SA range
//              saIndexRight - upper bound of the resulting SA range
//=============================================================================
void BWTSARangeBackward(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned long long *saIndexLeft, unsigned long long *saIndexRight);

//=============================================================================
//  BWTSARangeBackward_Bidirection
//  This function performs the 2BWT backward search for 1 character
//  and returns the resulting SA ranges to the input parameter 
//  saIndexLeft and saIndexRight.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              c       - the character to search for
//
//              saIndexLeft  - lower bound of the intermediate SA ranges
//              saIndexRight - upper bound of the intermediate SA ranges
//              rev_saIndexLeft  - lower bound of the intermediate Rev SA ranges
//              rev_saIndexRight - upper bound of the intermediate Rev SA ranges
//
//  Output :    saIndexLeft  - lower bound of the resulting SA range
//              saIndexRight - upper bound of the resulting SA range
//              rev_saIndexLeft  - lower bound of the resulting SA range
//              rev_saIndexRight - upper bound of the resulting SA range
//=============================================================================
void BWTSARangeBackward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned long long *saIndexLeft, unsigned long long *saIndexRight,
                        unsigned long long *rev_saIndexLeft, unsigned long long *rev_saIndexRight);


//=============================================================================
//  BWTSARangeForward_Bidirection
//  This function performs the 2BWT forward search for 1 character
//  and returns the resulting SA ranges to the input parameter 
//  saIndexLeft and saIndexRight.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              c       - the character to search for
//
//              saIndexLeft  - lower bound of the intermediate SA ranges
//              saIndexRight - upper bound of the intermediate SA ranges
//              rev_saIndexLeft  - lower bound of the intermediate Rev SA ranges
//              rev_saIndexRight - upper bound of the intermediate Rev SA ranges
//
//  Output :    saIndexLeft  - lower bound of the resulting SA range
//              saIndexRight - upper bound of the resulting SA range
//              rev_saIndexLeft  - lower bound of the resulting SA range
//              rev_saIndexRight - upper bound of the resulting SA range
//=============================================================================

void BWTSARangeForward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned long long *saIndexLeft, unsigned long long *saIndexRight,
                        unsigned long long *rev_saIndexLeft, unsigned long long *rev_saIndexRight);
                        
//=============================================================================
//  BWTAllSARangesBackward
//  This function performs the 2BWT backward search for all characters in the 
//  alphabet set defined and supported, and returns the resulting SA ranges 
//  to the input parameter resultSaIndexLeft and resultSaIndexRight.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              saIndexLeft  - lower bound of the intermediate SA ranges
//              saIndexRight - upper bound of the intermediate SA ranges
//
//  Output :    resultSaIndexLefts  - lower bound of the resulting SA ranges.
//                                    Array in size of ALPHABET_SIZE.
//              resultSaIndexRights - upper bound of the resulting SA ranges.
//                                    Array in size of ALPHABET_SIZE.
//=============================================================================

void BWTAllSARangesBackward(Idx2BWT * idx2BWT, 
                        const unsigned long long saIndexLeft, const unsigned long long saIndexRight,
                        unsigned long long *resultSaIndexesLeft, unsigned long long *resultSaIndexesRight);
//=============================================================================
//  BWTAllSARangesBackward_Bidirection
//  This function performs the 2BWT backward search for all characters in the 
//  alphabet set defined and supported and returns the resulting SA ranges 
//  to the input parameter resultSaIndexLeft and resultSaIndexRight.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              c       - the character to search for
//
//              saIndexLeft  - lower bound of the intermediate SA ranges
//              saIndexRight - upper bound of the intermediate SA ranges
//              rev_saIndexLeft  - lower bound of the intermediate Rev SA ranges
//              rev_saIndexRight - upper bound of the intermediate Rev SA ranges
//
//  Output :    resultSaIndexLeft  - lower bound of the resulting SA range
//                                   Array in size of ALPHABET_SIZE.
//              resultSaIndexRight - upper bound of the resulting SA range
//                                   Array in size of ALPHABET_SIZE.
//              rev_resultSaIndexLeft  - lower bound of the resulting SA range
//                                       Array in size of ALPHABET_SIZE.
//              rev_resultSaIndexRight - upper bound of the resulting SA range
//                                       Array in size of ALPHABET_SIZE.
//=============================================================================
void BWTAllSARangesBackward_Bidirection(Idx2BWT * idx2BWT, 
                        const unsigned long long saIndexLeft, const unsigned long long saIndexRight,
                        const unsigned long long rev_saIndexLeft, const unsigned long long rev_saIndexRight,
                        unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight,
                        unsigned long long *rev_resultSaIndexLeft, unsigned long long *rev_resultSaIndexRight);

//=============================================================================
//  BWTAllSARangesForward_Bidirection
//  This function performs the 2BWT forward search for all characters in the 
//  alphabet set defined and supported and returns the resulting SA ranges 
//  to the input parameter resultSaIndexLeft and resultSaIndexRight.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              c       - the character to search for
//
//              saIndexLeft  - lower bound of the intermediate SA ranges
//              saIndexRight - upper bound of the intermediate SA ranges
//              rev_saIndexLeft  - lower bound of the intermediate Rev SA ranges
//              rev_saIndexRight - upper bound of the intermediate Rev SA ranges
//
//  Output :    resultSaIndexLeft  - lower bound of the resulting SA range
//                                   Array in size of ALPHABET_SIZE.
//              resultSaIndexRight - upper bound of the resulting SA range
//                                   Array in size of ALPHABET_SIZE.
//              rev_resultSaIndexLeft  - lower bound of the resulting SA range
//                                       Array in size of ALPHABET_SIZE.
//              rev_resultSaIndexRight - upper bound of the resulting SA range
//                                       Array in size of ALPHABET_SIZE.
//=============================================================================
void BWTAllSARangesForward_Bidirection(Idx2BWT * idx2BWT, 
                        const unsigned long long saIndexLeft, const unsigned long long saIndexRight,
                        const unsigned long long rev_saIndexLeft, const unsigned long long rev_saIndexRight,
                        unsigned long long *resultSaIndexLeft, unsigned long long *resultSaIndexRight,
                        unsigned long long *rev_resultSaIndexLeft, unsigned long long *rev_resultSaIndexRight);
                        

//=============================================================================
//  BWTRetrievePositionFromSAIndex
//  This function translates an SA index (an index within a reported SA range)
//  into a text position in the original reference sequence in the format of
//  1-based sequenceId and offset.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              saIndex - the SA index within an SA range
//  Output :    sequenceId  - 1-based sequence id in the original reference
//                            sequence. Value can be 1,2,3....
//              offset      - 1-based offset in the above sequence. Value 
//                            can be 1,2,3...
//=============================================================================
void BWTRetrievePositionFromSAIndex(Idx2BWT * idx2BWT, 
                                    unsigned long long saIndex, 
                                    unsigned int * sequenceId, unsigned long long * offset);
                                    
#endif

