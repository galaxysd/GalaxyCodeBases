/*

   overlap-api.cpp        Basic BWT operations in BASE.

   This program contains basic functions for BWT operations.

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

   Date   : 1st Jan 2015
   Author : Dinghua LI
   Change : developped this progame for BASE.
*/
#include "overlap-api.h"

inline unsigned int _getBWTvalue(BWT *bwt, unsigned long long i) {
    static unsigned int bwtCharMask = ((1 << BIT_PER_CHAR) - 1) & 7;
    unsigned int c = bwt->bwtCode[i / CHAR_PER_WORD] >> (BITS_IN_WORD - BIT_PER_CHAR * (i % CHAR_PER_WORD + 1));
    return c & bwtCharMask;
}

int _dfsExtractReadInf(Idx2BWT * idx2BWT, ULL saL, ULL saR, ReadInf* &ri, int outputLimit, int strand, int _pos)
{
    ULL l, r;
    unsigned int c;
    ULL numSAs = 0;
    int totalCount = 0;
    l = saL, r = saR;

    if (saL == saR) { // only one route, backward search to the end
        do {
            c = _getBWTvalue(idx2BWT->bwt, saL);
            saL = idx2BWT->bwt->cumulativeFreq[c] + BWTOccValue(idx2BWT->bwt, saL, c);
            ++_pos;
        } while (c != 4);

        ri->read_id = idx2BWT->readIDtable[saL - idx2BWT->bwt->cumulativeFreq[4]];
        ri->strand = strand;
        ri->pos = _pos - 1;
        ++ri;
        return 1;
    } 

    for (int i = 0; i < 4; ++i) { // try A C G T
        BWTSARangeBackward(idx2BWT, i ,&saL, &saR);
        if (saL <= saR && saR-saL <= r-l) {
            int outputCount = _dfsExtractReadInf(idx2BWT, saL, saR, ri, outputLimit, strand, _pos + 1);
            outputLimit -= outputCount;
            totalCount += outputCount;
            if (outputLimit <= 0) return totalCount;
            
            numSAs += saR - saL + 1;
            if (numSAs == r - l + 1) { // no branch at all
                return totalCount;
            }
        }
        saL = l; saR = r;  // reset ranges
    }

    // try add $ 
    BWTSARangeBackward(idx2BWT, 4 ,&saL, &saR); // 4 is $
    if (saL <= saR) {
        int limit = saR - saL + 1;
        if (limit > outputLimit) {
            limit = outputLimit;
        }

        for (int i = 0; i < limit; ++i) {
            ri->read_id = idx2BWT->readIDtable[saL + i - idx2BWT->bwt->cumulativeFreq[4]];
            ri->strand = strand;
            ri->pos = _pos;
            ++ri;
        }
        totalCount += limit;
    }

    return totalCount;
}

int extractReadInf(Idx2BWT * idx2BWT, ExactRange range, ReadInf* ri, int outputLimit)
{
    return _dfsExtractReadInf(idx2BWT, range.saL, range.saR, ri, outputLimit, range.strand, 0);
}

/*
    Function: Find the sequences and quality of a given read id.
    Input: read id, sequence, quality.
    Output: updated sequence and quality.
    Comment: This is a function learn from SGA:src/SuffixTools/BWTAlgorithms.cpp:extractString, and it will be used in the step to generate contig sequences.
*/

int extractRead(Idx2BWT * idx2BWT, unsigned long long readID, char* seq, char* qual) 
{
    static char dnaChars[] = {'A','C','G','T'};
    int readLen = 0;
    BWT *bwt = idx2BWT->bwt;
    ULL l = idx2BWT->bwt->cumulativeFreq[4] + readID;
    for (;;) {
        unsigned int c = _getBWTvalue(idx2BWT->bwt, l);
        if (c == 4) {
          for (int i = 0, j = readLen - 1; i < j; ++i, --j) {
              seq[i] ^= seq[j];
              seq[j] ^= seq[i];
              seq[i] ^= seq[j];
          }
          seq[readLen] = 0;
          l = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c);
          fprintf(stderr, "ReadID verify: %d. Just for debug, comment this at file %s, line: %d.\n", idx2BWT->readIDtable[l - bwt->cumulativeFreq[4]], __FILE__, __LINE__);
          return readLen;
        } else {
          l = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c);
          seq[readLen] = dnaChars[c];
          ++readLen;
        }
    }
}
