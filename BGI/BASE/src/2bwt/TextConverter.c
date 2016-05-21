/*

   TextConverter.c		Text Converter

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TextConverter.h"
#include "MiscUtilities.h"
#include "r250.h"


unsigned long long GetWordPackedText(const unsigned int *packedText, const unsigned long long index, const unsigned long long shift, const unsigned long long numberOfBit, const unsigned long long vacantBit) {

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
#ifdef DNA_ONLY
		text = (packedText[index] << shift) | (packedText[index + 1] >> (BITS_IN_WORD - shift));
#else
		text = (packedText[index] << shift) | (packedText[index + 1] >> (BITS_IN_WORD - shift) << vacantBit);
#endif
	} else {
		text = packedText[index];
	}

	if (numberOfBit < BITS_IN_WORD) {
		// Fill unused bit with zero
		text &= mask[numberOfBit];
	}

	return text;
}


unsigned long long ReadCharMap(unsigned char *charMap, const char *inputFileName, const unsigned char defaultMapping) {

	FILE *inputFile;
	char c;
	unsigned long long v, alphabetSize;

	inputFile = (FILE*)fopen64(inputFileName, "r");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadCharMap() : Cannot open character map!\n");
		exit(1);
	}

	for (v=0; v<CHAR_MAP_SIZE; v++) {
		charMap[v] = defaultMapping;
	}

	alphabetSize = 0;

	while (!feof(inputFile)) {
		fscanf(inputFile, " %c %llu \n", &c, &v);
		if (v > CHAR_MAP_SIZE) {
			fprintf(stderr, "ReadCharMap() : Invalid charMap!\n");
			return 0;
		}
		charMap[(unsigned int)c] = (unsigned char)v;
		if (v > alphabetSize) {
			alphabetSize = v;
		}
	}

	fclose(inputFile);

	alphabetSize++;

	return alphabetSize;

}

void GenerateReverseCharMap(const unsigned char *charMap, unsigned char *reverseCharMap) {

	unsigned long long i, j;

	for (i=0; i<CHAR_MAP_SIZE; i++) {
		reverseCharMap[i] = INVALID_CHAR;
		for (j=0; j<CHAR_MAP_SIZE; j++) {
			if (charMap[j] == i) {
				reverseCharMap[i] = (unsigned char)j;
				break;
			}
		}
	}

}

unsigned int BitPerWordPackedChar(const unsigned int alphabetSize) {

	#ifdef DEBUG
	if (alphabetSize < 2) {
		fprintf(stderr, "BitPerWordPackedChar() : alphabetSize < 2!\n");
		exit(1);
	}
	#endif
	
	//return ceilLog2(alphabetSize);
	return BitPerBytePackedChar(alphabetSize);

}

unsigned long long TextLengthFromWordPacked(unsigned long long wordPackedLength, unsigned long long bitPerChar, unsigned long long lastWordLength) {

	return (wordPackedLength - 1) * (BITS_IN_WORD / bitPerChar) + lastWordLength;

}

unsigned long long WordPackedLengthFromText(unsigned long long textLength, unsigned int bitPerChar) {

	return (textLength + (BITS_IN_WORD / bitPerChar) - 1) / (BITS_IN_WORD / bitPerChar);

}

unsigned int LastWordLength(unsigned long long textLength, unsigned int bitPerChar) {

	return textLength % (BITS_IN_WORD / bitPerChar);

}

unsigned int BitPerBytePackedChar(const unsigned int alphabetSize) {

	unsigned int bitPerChar;

	#ifdef DEBUG
	if (alphabetSize < 2) {
		fprintf(stderr, "BitPerBytePackedChar() : alphabetSize < 2!\n");
		exit(1);
	}
	#endif

	bitPerChar = ceilLog2(alphabetSize);

	#ifdef DEBUG
	if (bitPerChar > BITS_IN_BYTE) {
		fprintf(stderr, "BitPerBytePackedChar() : bitPerChar > BITS_IN_BYTE!\n");
		exit(1);
	}
	#endif

	// Return the largest number of bit that does not affect packing efficiency
	if (BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar) > bitPerChar) {
		bitPerChar = BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar);
	}
	if (bitPerChar<4) return 4;
	return bitPerChar;
}

unsigned long long TextLengthFromBytePacked(unsigned long long bytePackedLength, unsigned int bitPerChar, unsigned int lastByteLength) {

	if (bytePackedLength > ALL_ONE_MASK_64 / (BITS_IN_BYTE / bitPerChar)) {
		fprintf(stderr, "TextLengthFromBytePacked(): text length > 2^32!\n");
		exit(1);
	}
	return (bytePackedLength - 1) * (BITS_IN_BYTE / bitPerChar) + lastByteLength;

}

unsigned long long BytePackedLengthFromText(unsigned long long textLength, unsigned int bitPerChar) {

	return (textLength + (BITS_IN_BYTE / bitPerChar) - 1) / (BITS_IN_BYTE / bitPerChar);

}

unsigned char LastByteLength(unsigned long long textLength, unsigned int bitPerChar) {

	return (unsigned char)(textLength % (BITS_IN_BYTE / bitPerChar));

}

void ConvertTextToWordPacked(const unsigned char *input, unsigned int *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned long long textLength) {

	unsigned int bitPerChar, charPerWord;
	unsigned long long i, j, k;
	unsigned int c;
	unsigned long long charValue;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	for (i=0; i<textLength/charPerWord; i++) {
		c = 0;
		j = i * charPerWord;
		for (k=0; k<charPerWord; k++) {
			charValue = charMap[input[j+k]];
			if (charValue >= alphabetSize) {
				charValue = 0;
			}
			c = c | (charValue << (BITS_IN_WORD - (k+1) * bitPerChar));
		}
		output[i] = c;
	}
	if (i * charPerWord < textLength) {
		c = 0;
		j = i * charPerWord;
		for (k=0; j+k < textLength; k++) {
			charValue = charMap[input[j+k]];
			if (charValue >= alphabetSize) {
				charValue = 0;
			}
			c = c | (charValue << (BITS_IN_WORD - (k+1) * bitPerChar));
		}
		output[i] = c;
	}

}

void ConvertTextToBytePacked(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned long long textLength) {

	unsigned int bitPerChar, charPerByte;
	unsigned long long i, j, k;
	unsigned char c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = 0;
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			c = c | (unsigned char)(charMap[input[j+k]] << (BITS_IN_BYTE - (k+1) * bitPerChar));
		}
		output[i] = c;
	}
	if (i * charPerByte < textLength) {
		c = 0;
		j = i * charPerByte;
		for (k=0; j+k < textLength; k++) {
			c = c | (unsigned char)(charMap[input[j+k]] << (BITS_IN_BYTE - (k+1) * bitPerChar));
		}
		output[i] = c;
	}

}

void ConvertWordPackedToText(const unsigned int *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned long long textLength) {

	unsigned int bitPerChar, charPerWord;
	unsigned long long i, j, k;
	unsigned int c;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	for (i=0; i<textLength/charPerWord; i++) {
		c = input[i];
		j = i * charPerWord;
		for (k=0; k<charPerWord; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_WORD - bitPerChar)];
			c <<= bitPerChar;
		}
	}
	if (i * charPerWord < textLength) {
		c = input[i];
		j = i * charPerWord;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_WORD - bitPerChar)];
			c <<= bitPerChar;
		}
	}

}

void ConvertBytePackedToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned int alphabetSize, const unsigned long long textLength) {

	unsigned int bitPerChar, charPerByte;
	unsigned long long i, j, k;
	unsigned char c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_BYTE - bitPerChar)];
			c <<= bitPerChar;
		}
	}
	if (i * charPerByte < textLength) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_BYTE - bitPerChar)];
			c <<= bitPerChar;
		}
	}

}

void ConvertBytePackedToCode(const unsigned char *input, unsigned char *output, const unsigned int alphabetSize, const unsigned long long textLength) {

	unsigned int bitPerChar, charPerByte;
	unsigned long long i, j, k;
	unsigned char c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			output[j+k] = c >> (unsigned char)(BITS_IN_BYTE - bitPerChar);
			c <<= bitPerChar;
		}
	}
	if (i * charPerByte < textLength) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = c >> (unsigned char)(BITS_IN_BYTE - bitPerChar);
			c <<= bitPerChar;
		}
	}

}

void ConvertWordPackedToBytePacked(const unsigned int *input, unsigned char *output, const unsigned int alphabetSize, const unsigned long long textLength) {

	unsigned long long i, j, k;
	unsigned int c;
	unsigned long long bitPerBytePackedChar;
	unsigned long long bitPerWordPackedChar;
	unsigned int charPerWord;
	unsigned int charPerByte;
	unsigned long long bytePerIteration;
	unsigned long long byteProcessed = 0;
	unsigned long long wordProcessed = 0;
	unsigned long long mask, shift;
	
	unsigned long long buffer[BITS_IN_WORD];

	bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
	bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerBytePackedChar;
	charPerByte = BITS_IN_BYTE / bitPerWordPackedChar;

	bytePerIteration = charPerWord / charPerByte;
	mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
	shift = BITS_IN_WORD - bitPerWordPackedChar;

	while ((wordProcessed + 1) * charPerWord < textLength) {

		c = input[wordProcessed];
		for (i=0; i<charPerWord; i++) {
			buffer[i] = c >> shift;
			c <<= bitPerWordPackedChar;
		}
		wordProcessed++;

		k = 0;
		for (i=0; i<bytePerIteration; i++) {
			c = 0;
			for (j=0; j<charPerByte; j++) {
				c |= buffer[k] << (BITS_IN_BYTE - (j+1) * bitPerBytePackedChar);
				k++;
			}
			output[byteProcessed] = (unsigned char)c;
			byteProcessed++;
		}

	}

	c = input[wordProcessed];
	for (i=0; i < textLength - wordProcessed * charPerWord; i++) {
		buffer[i] = c >> shift;
		c <<= bitPerWordPackedChar;
	}

	k = 0;
	while (byteProcessed * charPerByte < textLength) {
		c = 0;
		for (j=0; j < textLength - wordProcessed * charPerWord; j++) {
			c |= buffer[k] << (BITS_IN_BYTE - (j+1) * bitPerBytePackedChar);
			k++;
		}
		output[byteProcessed] = (unsigned char)c;
		byteProcessed++;
	}

}

void ConvertBytePackedToWordPacked(const unsigned char *input, unsigned int *output, const unsigned int alphabetSize, const unsigned long long textLength) {

	unsigned long long i, j, k;
	unsigned int c;
	unsigned long long bitPerBytePackedChar;
	unsigned long long bitPerWordPackedChar;
	unsigned int charPerWord;
	unsigned int charPerByte;
	unsigned long long bytePerIteration;
	unsigned long long byteProcessed = 0;
	unsigned long long wordProcessed = 0;
	unsigned long long mask, shift;
	
	unsigned long long buffer[BITS_IN_WORD];

	bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
	bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerBytePackedChar;
	charPerWord = BITS_IN_WORD / bitPerWordPackedChar;

	bytePerIteration = charPerWord / charPerByte;
	mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
	shift = BITS_IN_WORD - BITS_IN_BYTE + bitPerBytePackedChar - bitPerWordPackedChar;

	while ((wordProcessed + 1) * charPerWord < textLength) {

		k = 0;
		for (i=0; i<bytePerIteration; i++) {
			c = (unsigned int)input[byteProcessed] << shift;
			for (j=0; j<charPerByte; j++) {
				buffer[k] = c & mask;
				c <<= bitPerBytePackedChar;
				k++;
			}
			byteProcessed++;
		}

		c = 0;
		for (i=0; i<charPerWord; i++) {
			c |= buffer[i] >> bitPerWordPackedChar * i;
		}
		output[wordProcessed] = c;
		wordProcessed++;

	}

	k = 0;
	for (i=0; i < (textLength - wordProcessed * charPerWord - 1) / charPerByte + 1; i++) {
		c = (unsigned int)input[byteProcessed] << shift;
		for (j=0; j<charPerByte; j++) {
			buffer[k] = c & mask;
			c <<= bitPerBytePackedChar;
			k++;
		}
		byteProcessed++;
	}

	c = 0;
	for (i=0; i<textLength - wordProcessed * charPerWord; i++) {
		c |= buffer[i] >> bitPerWordPackedChar * i;
	}
	output[wordProcessed] = c;

}

void ConvertTextToCode(const unsigned char *input, unsigned char *output, const unsigned char *charMap, const unsigned long long textLength) {

	unsigned long long i;

	for (i=0; i< textLength; i++) {
		output[i] = charMap[input[i]];
	}

}

void ConvertCodeToText(const unsigned char *input, unsigned char *output, const unsigned char *reverseCharMap, const unsigned long long textLength) {

	unsigned long long i;

	for (i=0; i< textLength; i++) {
		output[i] = reverseCharMap[input[i]];
	}

}

void PackTextWithAllShift(const unsigned char *input, unsigned int **output, const unsigned char *charMap, const unsigned int alphabetSize, const unsigned long long textLength) {

	unsigned int bitPerChar;
	unsigned long long numberOfShift;
	unsigned long long numberOfWord;
	unsigned long long shift;

	unsigned long long i, j;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	numberOfShift = BITS_IN_WORD / bitPerChar;
	numberOfWord = WordPackedLengthFromText(textLength, bitPerChar);

	ConvertTextToWordPacked(input, output[0], charMap, alphabetSize, textLength);

	for (i=1; i<numberOfShift; i++) {
		shift = i * bitPerChar;
		output[i][0] = output[0][0] >> shift;
		for (j=1; j<=numberOfWord; j++) {
			output[i][j] = (output[0][j] >> shift) | (output[0][j-1] << (BITS_IN_WORD - shift));
		}
	}

}


unsigned long long ReadTextAsWordPacked(const char *inputFileName, const unsigned char *charMap, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned long long maxTextLength) {

	FILE *inputFile;
	unsigned char *buffer;
	unsigned int charPerWord;
	unsigned long long charRead;
	unsigned long long charProcessed = 0, wordProcessed = 0;
	unsigned long long charPerBuffer;

	inputFile = (FILE*)fopen64(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadTextAsWordPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	charPerWord = BITS_IN_WORD / BitPerWordPackedChar(alphabetSize);
	charPerBuffer = PACKED_BUFFER_SIZE / charPerWord * charPerWord;

	buffer = (unsigned char*) MMUnitAllocate(charPerBuffer);

	charRead = (unsigned int)fread(buffer, 1, charPerBuffer, inputFile);
	while (charRead > 0 && charProcessed + charRead < maxTextLength) {
		ConvertTextToWordPacked(buffer, targetAddress + wordProcessed, charMap, alphabetSize, charRead);
		wordProcessed += charRead / charPerWord;
		charProcessed += charRead;
		charRead = (unsigned int)fread(buffer, 1, charPerBuffer, inputFile);
	}

	if (charRead > 0 && charProcessed < maxTextLength) {
		ConvertTextToWordPacked(buffer, targetAddress + wordProcessed, charMap, alphabetSize, min(charRead, maxTextLength - charProcessed));
		charProcessed += charRead;
	}

	MMUnitFree(buffer, charPerBuffer);

	fclose(inputFile);

	return charProcessed;

}

unsigned long long ReadBytePackedAsWordPacked(const char *inputFileName, const unsigned int alphabetSize, unsigned int *targetAddress, const unsigned long long maxTextLength) {

	FILE *inputFile;
	unsigned char *buffer1, *buffer2;
	unsigned int charPerByte, charPerWord;
	unsigned long long charPerBuffer, wordPerBuffer;
	unsigned long long charProcessed = 0, wordProcessed = 0;
	unsigned long long byteRead, tempByteRead;
	unsigned long long charInLastBuffer;
	unsigned long long bufferSize;

	inputFile = (FILE*)fopen64(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadBytePackedAsWordPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	charPerByte = BITS_IN_BYTE / BitPerBytePackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / BitPerWordPackedChar(alphabetSize);
	bufferSize = PACKED_BUFFER_SIZE / charPerByte / charPerWord * charPerByte * charPerWord;

	charPerBuffer = bufferSize * charPerByte;
	wordPerBuffer = charPerBuffer / charPerWord;

	buffer1 = (unsigned char*) MMUnitAllocate(bufferSize);
	buffer2 = (unsigned char*) MMUnitAllocate(bufferSize);

	byteRead = (unsigned int)fread(buffer1, 1, bufferSize, inputFile);
	tempByteRead = (unsigned int)fread(buffer2, 1, bufferSize, inputFile);

	while (tempByteRead > 1 && charProcessed + charPerBuffer < maxTextLength) {
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, charPerBuffer);
		charProcessed += charPerBuffer;
		wordProcessed += wordPerBuffer;
		memcpy(buffer1, buffer2, bufferSize);
		byteRead = tempByteRead;
		tempByteRead = (unsigned int)fread(buffer2, 1, bufferSize, inputFile);
	}

	if (tempByteRead > 1) {
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, maxTextLength - charProcessed);
		charProcessed += charPerBuffer;
	} else {
		if (tempByteRead == 1) {
			charInLastBuffer = charPerBuffer - charPerByte + buffer2[0];
		} else {
			charInLastBuffer = (byteRead - 2) * charPerByte + buffer1[byteRead - 1];
		}
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, min(maxTextLength - charProcessed, charInLastBuffer));
		charProcessed += charInLastBuffer;
	}

	MMUnitFree(buffer1, bufferSize);
	MMUnitFree(buffer2, bufferSize);

	fclose(inputFile);

	return charProcessed;

}

// Alphabet size of DNA must be 8
void *DNALoadPacked(const char *inputFileName, unsigned long long *textLength, const unsigned int convertToWordPacked, const unsigned int trailerBufferInWord) {

	FILE *inputFile;
	unsigned char tempChar[4];
	unsigned int *packedText;
	off64_t packedFileLen;
	unsigned char lastByteLength;
	unsigned long long wordToProcess;
	unsigned long long i;
	unsigned long long trailerBufferIn128;

	trailerBufferIn128 = (trailerBufferInWord + 3) / 4 * 4;

	inputFile = (FILE*)(FILE*)fopen64(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "DNALoadPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(inputFile, -1, SEEK_END);
	packedFileLen = ftello64(inputFile);
	if (packedFileLen == -1) {
		fprintf(stderr, "DNALoadPacked(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, inputFile);

	*textLength = (packedFileLen - 1) * 2 + lastByteLength;

	wordToProcess = (*textLength + 32 - 1) / 32 * 4 + trailerBufferIn128 * 4;		// allocate multiple of 128 bit + trailer buffer

	packedText = (unsigned int*) MMUnitAllocate(wordToProcess * sizeof(unsigned int));
	for (i=(*textLength)/8; i<wordToProcess; i++) {
		packedText[i] = 0;
	}

	fseek(inputFile, 0, SEEK_SET);
	fread(packedText, 1, packedFileLen, inputFile);
	fclose(inputFile);

	if (convertToWordPacked) {

		for (i=0; i<wordToProcess; i++) {
	
			*(unsigned int*)tempChar = packedText[i];
			packedText[i] = (tempChar[0] << 24) | (tempChar[1] << 16) | (tempChar[2] << 8) | tempChar[3];

		}

	}

	return (void*)packedText;

}

void DNAFreePacked(void* packedDNA, const unsigned long long textLength, const unsigned int trailerBufferInWord) {

	unsigned int trailerBufferIn128;

	trailerBufferIn128 = (trailerBufferInWord + 3) / 4 * 4;

	MMUnitFree(packedDNA, ((textLength + 64 - 1) / 64 * 4 + trailerBufferIn128 * 4) * sizeof(unsigned int));

}

void SaveText(const char *outputFileName, const unsigned char *text, const unsigned long long textLength) {

	FILE *outputFile;

	outputFile = (FILE*)fopen64(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveText() : Cannot open output file!\n");
		exit(1);
	}

	fwrite(text, sizeof(unsigned char), textLength, outputFile);
	fclose(outputFile);

}

void SaveBytePacked(const char *outputFileName, const unsigned char *bytePacked, const unsigned long long textLength, const unsigned int alphabetSize) {

	FILE *outputFile;
	unsigned long long bitPerChar, charPerByte, bytePackedLen;
	unsigned char lastByteLen;
	unsigned char zero = 0;

	outputFile = (FILE*)fopen64(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveBytePacked() : Cannot open output file!\n");
		exit(1);
	}

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	bytePackedLen = BytePackedLengthFromText(textLength, bitPerChar);
	lastByteLen = LastByteLength(textLength, bitPerChar);

	fwrite(bytePacked, sizeof(unsigned char), bytePackedLen, outputFile);
	if (lastByteLen == 0) {
		fwrite(&zero, sizeof(unsigned char), 1, outputFile);
	}
	fwrite(&lastByteLen, sizeof(unsigned char), 1, outputFile);
	fclose(outputFile);

}

void SaveWordPacked(const char *outputFileName, const unsigned int *wordPacked, const unsigned long long textLength, const unsigned int alphabetSize) {

	FILE *outputFile;
	unsigned long long bitPerChar, charPerWord, wordPackedLen;
	unsigned long long lastWordLen;
	unsigned int zero = 0;

	outputFile = (FILE*)fopen64(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveWordPacked() : Cannot open output file!\n");
		exit(1);
	}

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	wordPackedLen = WordPackedLengthFromText(textLength, bitPerChar);
	lastWordLen = LastWordLength(textLength, bitPerChar);

	fwrite(wordPacked, sizeof(unsigned int), wordPackedLen, outputFile);
	if (lastWordLen == 0) {
		fwrite(&zero, sizeof(unsigned int), 1, outputFile);
	}
	fwrite(&lastWordLen, sizeof(unsigned int), 1, outputFile);
	fclose(outputFile);

}

FILE *InitialLoadPackedIncFromEnd(const char* inputFileName, unsigned char *packedOutput, const unsigned int alphabetSize, 
								  const unsigned long long packedLengthPerLoad, unsigned int *textLength, unsigned int *textLengthForThisLoad) {

	FILE *packedFile;
	off64_t packedFileLen;
	unsigned long long len, packedFileLenForThisLoad;
	unsigned char lastByteLength;
	unsigned int bitPerChar, charPerWord;

	packedFile = (FILE*)fopen64(inputFileName, "rb");

	if (packedFile == NULL) {
		fprintf(stderr, "InitialLoadPackedIncFromEnd() : Cannot open inputFileName!\n");
		exit(1);
	}

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	fseek(packedFile, -1, SEEK_END);
	packedFileLen = ftello64(packedFile);
	if (packedFileLen == -1) {
		fprintf(stderr, "InitialLoadPackedIncFromEnd(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, packedFile);

	len = TextLengthFromBytePacked(packedFileLen, bitPerChar, lastByteLength);

	if (lastByteLength == 0 && (packedFileLen - 1) % packedLengthPerLoad == 0) {
		packedFileLenForThisLoad = 0;
		fseek(packedFile, -((int)(2+packedLengthPerLoad)), SEEK_END);
		*textLength = len;
		*textLengthForThisLoad = 0;
		return packedFile;
	}

	if (packedFileLen % packedLengthPerLoad == 0) {
		packedFileLenForThisLoad = packedLengthPerLoad;
	} else {
		packedFileLenForThisLoad = packedFileLen % packedLengthPerLoad;
	}
	fseek(packedFile, -1, SEEK_END);

	fseek(packedFile, -((int)packedFileLenForThisLoad), SEEK_CUR);
	fread(packedOutput, sizeof(unsigned char), packedFileLenForThisLoad, packedFile);
	fseek(packedFile, -((int)packedFileLenForThisLoad), SEEK_CUR);
	if (packedFileLen > packedFileLenForThisLoad) {
		fseek(packedFile, -((int)packedLengthPerLoad), SEEK_CUR);
	}

	*textLength = len;
	*textLengthForThisLoad = TextLengthFromBytePacked(packedFileLenForThisLoad, bitPerChar, lastByteLength);

	return packedFile;

}

void LoadPackedIncFromEnd(FILE *packedFile, unsigned char *packedOutput, const unsigned long long packedLengthPerLoad) {
	
	fread(packedOutput, sizeof(unsigned char), packedLengthPerLoad, packedFile);
	fseek(packedFile, -(2*(int)packedLengthPerLoad), SEEK_CUR);

}


FILE *InitialLoadTextIncFromEnd(const char* inputFileName, unsigned char *textOutput, const unsigned long long textLengthPerLoad, unsigned int *textLength, unsigned int *textLengthForThisLoad) {

	FILE *textFile;
	unsigned long long textLenForThisLoad;
	off64_t len;

	textFile = (FILE*)fopen64(inputFileName, "rb");

	if (textFile == NULL) {
		fprintf(stderr, "InitialLoadTextIncFromEnd() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(textFile, 0, SEEK_END);
	len = ftello64(textFile);
	if (len == -1) {
		fprintf(stderr, "InitialLoadTextIncFromEnd(): Cannot determine file length!\n");
		exit(1);
	}

	textLenForThisLoad = len % textLengthPerLoad;

	if (textLenForThisLoad > 0) {
		fseek(textFile, -((int)textLenForThisLoad), SEEK_END);
		fread(textOutput, sizeof(unsigned char), textLenForThisLoad, textFile);
		fseek(textFile, -((int)textLenForThisLoad), SEEK_END);
	}

	*textLength = len;
	*textLengthForThisLoad = textLenForThisLoad;

	return textFile;
}

void LoadTextIncFromEnd(FILE *textFile, unsigned char *textOutput, const unsigned long long textLengthPerLoad) {

	if (ftello64(textFile) < (off64_t) textLengthPerLoad) {
		fprintf(stderr, "LoadTextIncFromEnd(): file pointer is not correctly placed!\n");
		exit(1);
	}

	fseek(textFile, -((int)textLengthPerLoad), SEEK_CUR);
	fread(textOutput, sizeof(unsigned char), textLengthPerLoad, textFile);
	fseek(textFile, -((int)textLengthPerLoad), SEEK_CUR);

}

