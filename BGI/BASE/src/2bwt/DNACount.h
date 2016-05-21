/*

   DNACount.h        DNA Count

   This module contains DNA occurrence counting functions. The DNA must be
   in word-packed format.

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

#ifndef __DNA_COUNT_H__
#define __DNA_COUNT_H__

#include "TypeNLimit.h"

// DNA
#define DNA_ALPHABET_SIZE            8
#define DNA_CHAR_PER_WORD            8
#define DNA_BIT_PER_CHAR            4

// DNA occurrence count table
#define DNA_OCC_CNT_TABLE_SIZE_IN_WORD    65536
#define DNA_OCC_SUM_EXCEPTION(sum)            ((sum & 0xeeeeeeeeeeeeeeef) == 0) //fefefeff

// DNA with 'n'
#define DNA_N_ALPHABET_SIZE            5
#define DNA_N_CHAR_PER_WORD            10
#define DNA_N_BIT_PER_CHAR            3

// DNA with 'n' occurrence count table
#define DNA_N_OCC_CNT_TABLE_SIZE_IN_WORD    32786


void GenerateDNAOccCountTable(unsigned long long *dnaDecodeTable);

// The following functions can only count up to 255 characters
unsigned long long ForwardDNAOccCount(const unsigned int* dna, const unsigned long long index, const unsigned int character, const unsigned long long* dnaDecodeTable);
unsigned long long BackwardDNAOccCount(const unsigned int* dna, const unsigned long long index, const unsigned int character, const unsigned long long* dnaDecodeTable);
void ForwardDNAAllOccCount(const unsigned int* dna, const unsigned long long index, unsigned long long* __restrict occCount, const unsigned long long* dnaDecodeTable);
void BackwardDNAAllOccCount(const unsigned int* dna, const unsigned long long index, unsigned long long* __restrict  occCount, const unsigned long long* dnaDecodeTable);
unsigned long long Forward1OccCount(const unsigned int* bitVector, const unsigned long long index, const unsigned long long* dnaDecodeTable);    // Count number of 1 bit
unsigned long long Backward1OccCount(const unsigned int* bitVector, const unsigned long long index, const unsigned long long* dnaDecodeTable); // Count number of 1 bit

// The following functions have no limit on the number of characters
unsigned long long ForwardDNAOccCountNoLimit(const unsigned int* dna, const unsigned long long index, const unsigned int character, const unsigned long long* dnaDecodeTable);
unsigned long long BackwardDNAOccCountNoLimit(const unsigned int* dna, const unsigned long long index, const unsigned int character, const unsigned long long* dnaDecodeTable);
void ForwardDNAAllOccCountNoLimit(const unsigned int* dna, const unsigned long long index, unsigned long long* __restrict  occCount, const unsigned long long* dnaDecodeTable);
void BackwardDNAAllOccCountNoLimit(const unsigned int* dna, const unsigned long long index, unsigned long long* __restrict  occCount, const unsigned long long* dnaDecodeTable);


void GenerateDNA_NOccCountTable(unsigned int *dnaDecodeTable);

// The following functions have no limit on the number of characters
unsigned long long ForwardDNA_NOccCount(const unsigned int* dna, const unsigned long long index, const unsigned int character, const unsigned long long* dnaDecodeTable);
unsigned long long BackwardDNA_NOccCount(const unsigned int* dna, const unsigned long long index, const unsigned int character, const unsigned long long* dnaDecodeTable);
void ForwardDNA_NAllOccCount(const unsigned int* dna, const unsigned long long index, unsigned long long* __restrict  occCount, const unsigned long long* dnaDecodeTable);
void BackwardDNA_NAllOccCount(const unsigned int* dna, const unsigned long long index, unsigned long long* __restrict  occCount, const unsigned long long* dnaDecodeTable);

// The following functions have no limit on the number of characters
unsigned long long ForwardDNAnOccCountNoLimit(const unsigned int* dna, const unsigned long long index, const unsigned int character, const unsigned long long* dnaDecodeTable);
unsigned long long BackwardDNA_NOccCountNoLimit(const unsigned int* dna, const unsigned long long index, const unsigned int character, const unsigned long long* dnaDecodeTable);
void ForwardDNA_NAllOccCountNoLimit(const unsigned int* dna, const unsigned long long index, unsigned long long* __restrict  occCount, const unsigned long long* dnaDecodeTable);
void BackwardDNA_NAllOccCountNoLimit(const unsigned int* dna, const unsigned long long index, unsigned long long* __restrict  occCount, const unsigned long long* dnaDecodeTable);

// The first character from startAddr is indexed as 1
// DNA_NAllOccCount only count occurrence from character 0 to 3

// The following functions work for any word packed text
unsigned long long ForwardOccCount(const unsigned int* packed, const unsigned long long index, const unsigned int character, const unsigned int alphabetSize);
unsigned long long BackwardOccCount(const unsigned int* packed, const unsigned long long index, const unsigned int character, const unsigned int alphabetSize);
void ForwardAllOccCount(const unsigned int* packed, const unsigned long long index, const unsigned int alphabetSize, unsigned long long* occCount);
void BackwardAllOccCount(const unsigned int* packed, const unsigned long long index, const unsigned int alphabetSize, unsigned long long* occCount);


#endif

