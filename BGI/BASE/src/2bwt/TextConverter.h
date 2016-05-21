/*

   TextConverter.h		Text Converter

   This module contains miscellaneous text conversion functions.

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

#ifndef __TEXTCONVERTOR_H__
#define __TEXTCONVERTOR_H__

#include "TypeNLimit.h"
#include "MemManager.h"

#define INVALID_CHAR 0xFF
#define CHAR_MAP_SIZE 256
#define PACKED_BUFFER_SIZE			(PACKED_BUFFER_SIZE_IN_WORD * BYTES_IN_WORD)
#define PACKED_BUFFER_SIZE_IN_WORD	65536
#define MAX_SEQ_NAME_LENGTH			256
#define RANDOM_SUBSTITUTE			'R'

// charMap is a char array of size 256. The index of the array is the input text value
// and the content of the array is the output text value. e.g. A -> 0, C -> 1
// If the value of an entry = INVALID_CHAR, the indexed text value is an invalid input

// Retrieve word packed text
unsigned long long GetWordPackedText(const unsigned int *packedText, const unsigned long long index, const unsigned long long shift, const unsigned long long numberOfBit, const unsigned long long vacantBit);

// Character map functions
unsigned long long ReadCharMap(unsigned char *charMap, const char *inputFileName, const unsigned char defaultMapping);
void GenerateReverseCharMap(const unsigned char *charMap, unsigned char *reverseCharMap);

// Word packed text functions
unsigned int BitPerWordPackedChar(const unsigned int alphabetSize);
unsigned long long TextLengthFromWordPacked(unsigned long long wordPackedLength, unsigned int bitPerChar, unsigned int lastWordLength);
unsigned long long WordPackedLengthFromText(unsigned long long textLength, unsigned int bitPerChar);
unsigned int LastWordLength(unsigned long long textLength, unsigned int bitPerChar);

// Byte packed text functions
unsigned int BitPerBytePackedChar(const unsigned int alphabetSize);
unsigned long long TextLengthFromBytePacked(unsigned long long bytePackedLength, unsigned int bitPerChar, unsigned int lastByteLength);
unsigned long long BytePackedLengthFromText(unsigned long long textLength, unsigned int bitPerChar);
unsigned char LastByteLength(unsigned long long textLength, unsigned int bitPerChar);

// Conversion functions
void ConvertTextToWordPacked(const unsigned char *input, unsigned int *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertTextToBytePacked(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertWordPackedToText(const unsigned int *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertBytePackedToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertBytePackedToCode(const unsigned char *input, unsigned char *output, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertWordPackedToBytePacked(const unsigned int *input, unsigned char *output, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertBytePackedToWordPacked(const unsigned char *input, unsigned int *output, const unsigned int alphabetSize, const unsigned long long textLength);
void ConvertTextToCode(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned long long textLength);
void ConvertCodeToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned long long textLength);

// Pack text with all shift
void PackTextWithAllShift(const unsigned char *input, unsigned int **output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned long long textLength);

// Full load function
unsigned long long ReadTextAsWordPacked(const char *inputFileName, const unsigned char *charMap, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned long long maxTextLength);
unsigned long long ReadBytePackedAsWordPacked(const char *inputFileName, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned long long maxTextLength);
void *DNALoadPacked(const char *inputFileName, unsigned long long *textLength, const unsigned int convertToWordPacked, const unsigned int trailerBufferInWord);
void DNAFreePacked(void* packedDna, const unsigned long long textLength, const unsigned int trailerBufferInWord);

// Save functions
void SaveText(const char *outputFileName, const unsigned char *text, const unsigned long long textLength);
void SaveBytePacked(const char *outputFileName, const unsigned char *wordPacked, const unsigned long long textLength, const unsigned int alphabetSize);
void SaveWordPacked(const char *outputFileName, const unsigned int *wordPacked, const unsigned long long textLength, const unsigned int alphabetSize);

// Incremental load functions (start from end of text)
FILE *InitialLoadPackedIncFromEnd(const char* inputFileName, unsigned char *packedOutput, const unsigned int alphabetSize, const unsigned long long packedLengthPerLoad, unsigned long long *textLength, unsigned long long *textLengthForThisLoad);
void LoadPackedIncFromEnd(FILE *packedFile, unsigned char *packedOutput, const unsigned long long packedLengthPerLoad);
FILE *InitialLoadTextIncFromEnd(const char* inputFileName, unsigned char *textOutput, const unsigned long long textLengthPerLoad, unsigned long long *textLength, unsigned long long *textLengthForThisLoad);
void LoadTextIncFromEnd(FILE *textFile, unsigned char *textOutput, const unsigned long long textLengthPerLoad);


#endif
