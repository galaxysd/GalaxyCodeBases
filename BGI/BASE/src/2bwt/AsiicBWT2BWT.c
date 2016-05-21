/*

   AsiicBWT2BWT.c        Prepare 2BWT files.

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
   Author : Edward MK Wu, CHi Man Liu, Dinghua, Li

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include "DNACount.h"
#include "HSP.h"
#include "BWT.h"
#include <zlib.h>

using namespace std;

#define BUFFER_SIZE 10485760

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

void GenerateDNAOccCountTable(unsigned long long *dnaDecodeTable) {

    unsigned long long i, j, c, t;

    for (i=0; i<DNA_OCC_CNT_TABLE_SIZE_IN_WORD; i++) {
        dnaDecodeTable[i] = 0;
        c = i;
        for (j=0; j<4; j++) {
            t = c & 0x0000000F;
            dnaDecodeTable[i] += (unsigned long long)1 << (t * 4); //Each alphabet gets 4 bits
            c >>= 4;
        }
    }

}

void BWTGenerateOccValueFromBwt(const unsigned int*  bwt, unsigned int* __restrict occValue, unsigned long long* __restrict occValueMajor,
                                const unsigned long long textLength, const unsigned long long*  decodeTable) {

    unsigned long long numberOfOccValueMajor, numberOfOccValue;
    unsigned long long wordBetweenOccValue;
    unsigned long long numberOfOccIntervalPerMajor;
    unsigned long long c;
    unsigned long long i, j, k;
    unsigned long long occMajorIndex;
    unsigned long long occIndex, bwtIndex;
    unsigned long long sum;
    unsigned int tempOccValue0[ALPHABET_SIZE], tempOccValue1[ALPHABET_SIZE];
    unsigned int tempOccValue2[ALPHABET_SIZE], tempOccValue3[ALPHABET_SIZE];

    wordBetweenOccValue = OCC_INTERVAL / CHAR_PER_WORD;

    // Calculate occValue
    numberOfOccValue = (textLength + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;                // Value at both end for bi-directional encoding
    numberOfOccIntervalPerMajor = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
    numberOfOccValueMajor = (numberOfOccValue + numberOfOccIntervalPerMajor - 1) / numberOfOccIntervalPerMajor;

    for (i=0;i<ALPHABET_SIZE;i++) {
        tempOccValue0[i] = 0;
        occValueMajor[i] = 0;
        tempOccValue1[i] = 0;
        tempOccValue2[i] = 0;
        tempOccValue3[i] = 0;
    }

    occIndex = 0;
    bwtIndex = 0;
    for (occMajorIndex=1; occMajorIndex<numberOfOccValueMajor; occMajorIndex++) {

        for (i=0; i<numberOfOccIntervalPerMajor/OCC_VALUE_PER_LONG; i++) {

            sum = 0;
            for (j=0;j<ALPHABET_SIZE;j++) {
                tempOccValue1[j] = tempOccValue0[j];
            }

            for (j=0; j<wordBetweenOccValue; j++) {
                c = bwt[bwtIndex];
                sum = decodeTable[c >> 16];
                for (k=0;k<ALPHABET_SIZE;k++) {
                    tempOccValue1[k] += (sum & 0x0000000F);    sum >>= 4;
                }
                sum = decodeTable[c & 0x0000FFFF];
                for (k=0;k<ALPHABET_SIZE;k++) {
                    tempOccValue1[k] += (sum & 0x0000000F);    sum >>= 4;
                }
                bwtIndex++;
            }
            
            sum = 0;
            for (j=0;j<ALPHABET_SIZE;j++) {
                tempOccValue2[j] = tempOccValue1[j];
            }

            for (j=0; j<wordBetweenOccValue; j++) {
                c = bwt[bwtIndex];
                sum = decodeTable[c >> 16];
                for (k=0;k<ALPHABET_SIZE;k++) {
                    tempOccValue2[k] += (sum & 0x0000000F);    sum >>= 4;
                }
                sum = decodeTable[c & 0x0000FFFF];
                for (k=0;k<ALPHABET_SIZE;k++) {
                    tempOccValue2[k] += (sum & 0x0000000F);    sum >>= 4;
                }
                bwtIndex++;
            }
            
            sum = 0;
            for (j=0;j<ALPHABET_SIZE;j++) {
                tempOccValue3[j] = tempOccValue2[j];
            }

            for (j=0; j<wordBetweenOccValue; j++) {
                c = bwt[bwtIndex];
                sum = decodeTable[c >> 16];
                for (k=0;k<ALPHABET_SIZE;k++) {
                    tempOccValue3[k] += (sum & 0x0000000F);    sum >>= 4;
                }
                sum = decodeTable[c & 0x0000FFFF];
                for (k=0;k<ALPHABET_SIZE;k++) {
                    tempOccValue3[k] += (sum & 0x0000000F);    sum >>= 4;
                }
                bwtIndex++;
            }
            
            for (k=0;k<ALPHABET_SIZE;k++) {
                occValue[occIndex * ALPHABET_SIZE * 2 + k * 2] = (tempOccValue0[k] << 16) | tempOccValue1[k];
                occValue[occIndex * ALPHABET_SIZE * 2 + k * 2 + 1] = (tempOccValue2[k] << 16) | tempOccValue3[k];
                tempOccValue0[k] = tempOccValue3[k];
            }
            sum = 0;

            occIndex++;

            for (j=0; j<wordBetweenOccValue; j++) {
                c = bwt[bwtIndex];
                sum = decodeTable[c >> 16];
                for (k=0;k<ALPHABET_SIZE;k++) {
                    tempOccValue0[k] += (sum & 0x0000000F);    sum >>= 4;
                }
                sum = decodeTable[c & 0x0000FFFF];
                for (k=0;k<ALPHABET_SIZE;k++) {
                    tempOccValue0[k] += (sum & 0x0000000F);    sum >>= 4;
                }
                bwtIndex++;
            }
        }

        for (k=0;k<ALPHABET_SIZE;k++) {
            occValueMajor[occMajorIndex * ALPHABET_SIZE + k] = occValueMajor[(occMajorIndex - 1) * ALPHABET_SIZE + k] + tempOccValue0[k];
            tempOccValue0[k] = 0;
            
        }

    }

    while (occIndex < (numberOfOccValue-1)/OCC_VALUE_PER_LONG) {
        sum = 0;
        
        for (k=0;k<ALPHABET_SIZE;k++) {
            tempOccValue1[k] = tempOccValue0[k];
        }

        for (j=0; j<wordBetweenOccValue; j++) {
            c = bwt[bwtIndex];
            sum = decodeTable[c >> 16];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue1[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            sum = decodeTable[c & 0x0000FFFF];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue1[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            bwtIndex++;
        }

        for (k=0;k<ALPHABET_SIZE;k++) {
            tempOccValue2[k] = tempOccValue1[k];
        }

        for (j=0; j<wordBetweenOccValue; j++) {
            c = bwt[bwtIndex];
            sum = decodeTable[c >> 16];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue2[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            sum = decodeTable[c & 0x0000FFFF];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue2[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            bwtIndex++;
        }
        

        for (k=0;k<ALPHABET_SIZE;k++) {
            tempOccValue3[k] = tempOccValue2[k];
        }

        for (j=0; j<wordBetweenOccValue; j++) {
            c = bwt[bwtIndex];
            sum = decodeTable[c >> 16];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue3[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            sum = decodeTable[c & 0x0000FFFF];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue3[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            bwtIndex++;
        }

        for (k=0;k<ALPHABET_SIZE;k++) {
            occValue[occIndex * ALPHABET_SIZE * 2 + k * 2] = (tempOccValue0[k] << 16) | tempOccValue1[k];
            occValue[occIndex * ALPHABET_SIZE * 2 + k * 2 + 1] = (tempOccValue2[k] << 16) | tempOccValue3[k];
            tempOccValue0[k] = tempOccValue3[k];
        }

        sum = 0;
        occIndex++;

        for (j=0; j<wordBetweenOccValue; j++) {
            c = bwt[bwtIndex];
            sum = decodeTable[c >> 16];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue0[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            sum = decodeTable[c & 0x0000FFFF];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue0[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            bwtIndex++;
        }
    }

    sum = 0;
    for (k=0;k<ALPHABET_SIZE;k++) {
        tempOccValue1[k] = tempOccValue0[k];
    }

    if (occIndex * OCC_VALUE_PER_LONG < numberOfOccValue - 1) {
        for (j=0; j<wordBetweenOccValue; j++) {
            c = bwt[bwtIndex];
            sum = decodeTable[c >> 16];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue1[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            sum = decodeTable[c & 0x0000FFFF];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue1[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            bwtIndex++;
        }

        for (k=0;k<ALPHABET_SIZE;k++) {
            tempOccValue2[k] = tempOccValue1[k];
        }

        for (j=0; j<wordBetweenOccValue; j++) {
            c = bwt[bwtIndex];
            sum = decodeTable[c >> 16];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue2[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            sum = decodeTable[c & 0x0000FFFF];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue2[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            bwtIndex++;
        }
        

        for (k=0;k<ALPHABET_SIZE;k++) {
            tempOccValue3[k] = tempOccValue2[k];
        }

        for (j=0; j<wordBetweenOccValue; j++) {
            c = bwt[bwtIndex];
            sum = decodeTable[c >> 16];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue3[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            sum = decodeTable[c & 0x0000FFFF];
            for (k=0;k<ALPHABET_SIZE;k++) {
                tempOccValue3[k] += (sum & 0x0000000F);    sum >>= 4;
            }
            bwtIndex++;
        }
    }
    
    for (k=0;k<ALPHABET_SIZE;k++) {
        occValue[occIndex * ALPHABET_SIZE * 2 + k * 2] = (tempOccValue0[k] << 16) | tempOccValue1[k];
        occValue[occIndex * ALPHABET_SIZE * 2 + k * 2 + 1] = (tempOccValue2[k] << 16) | tempOccValue3[k];
    }

}

unsigned long long BWTGenerateOccValueToFileFromBwt(unsigned int *bwt, unsigned long long *cumulativeFreq, unsigned long long inverseSa0, 
                                                    const char *occValueFileName, unsigned long long*  decodeTable) 
{
    FILE *occValueFile;
    //unsigned long long bwtFileSizeInWord, bwtResidentSizeInWord;
    unsigned long long textLength;
    unsigned int *occValue;
    unsigned long long *occValueMajor;
    unsigned long long occSizeInWord, occMajorSizeInWord;
    //unsigned long long i;

    textLength = cumulativeFreq[ALPHABET_SIZE - 1];

    // occValue File
    occValueFile = (FILE*)fopen64(occValueFileName, "wb");
    if (occValueFile == NULL) {
        fprintf(stderr, "BWTGenerateOccValueToFileFromBwt(): Cannot open occ value file!\n");
        exit(1);
    }

    fwrite(&inverseSa0, sizeof(unsigned long long), 1, occValueFile);
    fwrite(cumulativeFreq, sizeof(unsigned long long), ALPHABET_SIZE, occValueFile);

    occSizeInWord = BWTOccValueMinorSizeInWord(textLength);
    occMajorSizeInWord = BWTOccValueMajorSizeInWord(textLength);
    occValue = (unsigned int*) malloc(occSizeInWord * sizeof(unsigned int));
    occValueMajor = (unsigned long long*) malloc(occMajorSizeInWord * sizeof(unsigned long long));

    if (decodeTable == NULL) {
        decodeTable = (unsigned long long*) malloc(DNA_OCC_CNT_TABLE_SIZE_IN_WORD * sizeof(unsigned long long));
        GenerateDNAOccCountTable(decodeTable);
    }

    BWTGenerateOccValueFromBwt(bwt, occValue, occValueMajor, textLength, decodeTable);

    fwrite(occValue, sizeof(unsigned int), occSizeInWord, occValueFile);
    fwrite(occValueMajor, sizeof(unsigned long long), occMajorSizeInWord, occValueFile);
    fclose(occValueFile);

    free(occValue);
    free(occValueMajor);
    free(decodeTable);

    return textLength;

}

void setQuality(unsigned int *bwtCode, unsigned long long offset) {
    unsigned long long whichWord = offset / CHAR_PER_WORD;
    int posAtWord = offset % CHAR_PER_WORD;
    bwtCode[whichWord] |= (1 << (posAtWord * BIT_PER_CHAR + BIT_PER_CHAR - 1));
}

void print_usage(char **argv) {
    fprintf(stderr, "Usage: %s output_prefix [quality_file] [read_len] < *.bwt.ascii \n", argv[0]);
    fprintf(stderr, "Input:\n 1) ASIIC BWT should be input by stdin\n 2) If quality file and read length are provided, quality will be saved to BWT.\n");
    fprintf(stderr, "Output:\n output_prefix.bwt and output_prefix.fmv.\n");
}

int main(int argc, char **argv) {
	if (argc < 2) {
		print_usage(argv);
		exit(1);
	}

	int charMap[256];
	char *output_prefix = argv[1];
	char *bwt_name = (char*) malloc((strlen(output_prefix) + 5) * sizeof(char));
	char *fmv_name = (char*) malloc((strlen(output_prefix) + 5) * sizeof(char));
	char *buffer = (char*) malloc(sizeof(char) * BUFFER_SIZE);
	unsigned int packedShift[CHAR_PER_WORD];
	unsigned int curr_word = 0;
	int char_index = 0;
	unsigned long long bwt_count[ALPHABET_SIZE];
	unsigned long long inverseSa0 = 0x3FFFFFFFFFFFFFFFULL;
	vector<unsigned int> bwtCode;

	charMap['A'] = charMap['a'] = 0;
	charMap['C'] = charMap['c'] = 1;
	charMap['G'] = charMap['g'] = 2;
	charMap['T'] = charMap['t'] = 3;
	charMap['$'] = 4;

	memset(bwt_count, 0, sizeof(bwt_count));

	for (int i = 0; i < CHAR_PER_WORD; i++) {
        packedShift[i] = BITS_IN_WORD - (i+1) * BIT_PER_CHAR;
    }

    strcpy(bwt_name, output_prefix);
    strcat(bwt_name, ".bwt");
    strcpy(fmv_name, output_prefix);
    strcat(fmv_name, ".fmv");

	// read asiic bwt
	while (true) {
		size_t read_size = fread(buffer, sizeof(char), BUFFER_SIZE, stdin);

		if (read_size == 0) {
			break;
		}

		for (size_t i = 0; i < read_size; ++i) {
			bwt_count[charMap[(unsigned char)(buffer[i])]]++;
			curr_word |= charMap[(unsigned char)(buffer[i])] << packedShift[char_index];
			++char_index;
			if (char_index >= CHAR_PER_WORD) {
				bwtCode.push_back(curr_word);
				curr_word = 0;
				char_index = 0;
			}
		}
	}
	if (char_index != 0) {
		bwtCode.push_back(curr_word);
	}

	for (int i = 1; i < ALPHABET_SIZE; ++i) {
		bwt_count[i] += bwt_count[i - 1];
	}

    unsigned long long bwtResidentSizeInWord = BWTResidentSizeInWord(bwt_count[ALPHABET_SIZE - 1]);
    unsigned long long bwtTotalWordCount = bwtCode.size();
    for (unsigned long long i=bwtTotalWordCount; i<bwtResidentSizeInWord; i++) {
        bwtCode.push_back(0);
    }
    bwtCode.reserve(bwtCode.size());
    puts("Generating Occ Value...");
    BWTGenerateOccValueToFileFromBwt(&bwtCode[0], &bwt_count[0], 0x3FFFFFFFFFFFFFFFULL, fmv_name, NULL); 
    puts("Generate Occ Value...DONE!");

    if (argc >= 4) {
        unsigned long long idx = 0;
        FILE *qual_file = fopen64(argv[2], "rb");
        int read_len = atoi(argv[3]);
        //long long total_num_reads = bwt_count[ALPHABET_SIZE - 1] / (read_len + 1);
        assert(bwt_count[ALPHABET_SIZE - 1] % (read_len + 1) == 0);

        int words_per_read_qual = (read_len + BITS_IN_WORD - 1) / BITS_IN_WORD;
        const int kReadsPerBatch = (1 << 20); // 1M
        unsigned int *buffer = (unsigned int*) malloc(sizeof(unsigned int) * words_per_read_qual * kReadsPerBatch);

        puts("Reading quality file and update quality score...");
        while (true) {
            size_t num_words = fread(buffer, sizeof(unsigned int), kReadsPerBatch * words_per_read_qual, qual_file);
            size_t num_reads = num_words / words_per_read_qual;
            assert(num_words % words_per_read_qual == 0);

            if (num_words == 0) { break; }

            unsigned int *qual_p = buffer;
            for (size_t i = 0; i < num_reads; ++i) {
                int word_index = 0;
                for (int j = 0; j < read_len; ++j) {
                    if (*qual_p & (1 << (BITS_IN_WORD - 1 - word_index))) {
                        setQuality(&bwtCode[0], idx);
                    }
                    ++idx;
                    ++word_index;
                    if (word_index % BITS_IN_WORD == 0) {
                        word_index = 0;
                        ++qual_p;
                    }
                }
                if (word_index) {
                    ++qual_p;
                }
            }
        }
        puts("DONE!");
        free(buffer);
        fclose(qual_file);
    }

    FILE *bwt_file = fopen64(bwt_name, "wb");
    puts("Writing BWT Code...");
    fwrite(&inverseSa0, sizeof(unsigned long long), 1, bwt_file);
    for (int i = 0; i < ALPHABET_SIZE; ++i) {
        fwrite(&bwt_count[i], sizeof(unsigned long long), 1, bwt_file);
    }
    fwrite(&bwtCode[0], sizeof(unsigned int), bwtTotalWordCount, bwt_file);
    fclose(bwt_file);
    puts("DONE!");

	return 0;
}
