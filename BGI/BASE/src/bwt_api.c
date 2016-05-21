/*

   bwt_api.c       BWT operations in BASE.

   This program contains basic functions for BWT operation and CAS operation.

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
    Author : Binghang Liu, bhliu@cs.hku.hk
*/

#include "bwt_api.h"

/*Sparse bits operation functions*/
int check_and_set_bwt_sparse_thread_and_used(ULL id, int size, int threadId) 
{
    BWT * bwt = idx2BWT[size]->rev_bwt;
    unsigned int oldvalue, newvalue;
    unsigned int new_word = ((unsigned int)(1) << 31) | ((unsigned int) (1) << 27), blank_word=0;
    new_word |= id_to_sparse_word[threadId];
    blank_word = ~((unsigned int) (id_to_sparse_word[63]));

    while (1) {
        oldvalue = bwt-> bwtCode[id];
        if(bwt->bwtCode[id] >> 31)
            return 0;
        newvalue = oldvalue & blank_word | new_word;
        if(__sync_bool_compare_and_swap(&(bwt->bwtCode[id]), oldvalue, newvalue))
            return 1;
        fprintf(stderr, "check and set conf\n");
    }
    return 0;
}

int is_used_bwt_sparse(ULL id, int size)
{
    BWT * bwt = idx2BWT[size]->rev_bwt;
    return bwt->bwtCode[id] >> 31;
}

int get_threadId_bwt_sparse(ULL id, int size)
{
    BWT * bwt = idx2BWT[size]->rev_bwt;
    if(((bwt->bwtCode[id] >> 27) & 1) == 0)
        return -1;

    int j;
    unsigned int word = 0, value = bwt->bwtCode[id];
    for(j=5; j>=0; j--)
    {
       word = word << 1 | ((value >> (3 + 4*j)) & 1 );
    }
    return word ;
}

int clean_current_bwt_sparse(ULL id, int size)
{
    BWT * bwt = idx2BWT[size]->rev_bwt;
    unsigned int new_word = ~((unsigned int)(1) << 27), oldvalue, newvalue;

    while (1) {
        oldvalue = bwt-> bwtCode[id];
        if(((bwt->bwtCode[id] >> 27) & 1) == 0)
            return 1;
        newvalue = oldvalue & new_word;
        if(__sync_bool_compare_and_swap(&(bwt->bwtCode[id]), oldvalue, newvalue))
            return 1;
        fprintf(stderr, "clean current conf\n");
    }
    return 0;
}

int clean_use_bwt_sparse(ULL id, int size)
{
    BWT * bwt = idx2BWT[size]->rev_bwt;
    unsigned int new_word = ~((unsigned int)(1) << 31), oldvalue, newvalue;

    while (1) {
        oldvalue = bwt-> bwtCode[id];
        if(bwt->bwtCode[id] >> 31 == 0)
            return 1;
        newvalue = oldvalue & new_word;
        if(__sync_bool_compare_and_swap(&(bwt->bwtCode[id]), oldvalue, newvalue))
            return 1;
        fprintf(stderr, "clean use conf\n");
    }
    return 0;
}

/*One bit operation functions*/
int is_used_bwt(ULL id, int size)
{
    return BWTCheckUsed(idx2BWT[size], id);
}

void set_used_bwt(ULL id, int size, int flag)
{
    BWTSetUsed(idx2BWT[size], id, flag);
}

void check_flag_bwt(ULL id, int size, int* used, int* cons)
{
    BWTCheckFlag(idx2BWT[size], id, used, cons);
}

void set_flag_bwt(ULL id, int size, int flag, int cons)
{
    BWTSetFlag(idx2BWT[size], id, flag, cons);
}

/*Extract reads sequence and quality by id*/
inline unsigned int _getBWTvalue(BWT *bwt, unsigned long long i) {
    static unsigned int bwtCharMask = ((1 << BIT_PER_CHAR) - 1) & 7;
    unsigned int c = bwt->bwtCode[i / CHAR_PER_WORD] >> (BITS_IN_WORD - BIT_PER_CHAR * (i % CHAR_PER_WORD + 1));
    return c & bwtCharMask;
}


int extractReadQual(Idx2BWT * idx2BWT, unsigned long long readID, char* seq, int* qual)
{
    static char dnaChars[] = {'A','C','G','T'};
    int readLen = 0;
    BWT *bwt = idx2BWT->bwt;
    ULL l = idx2BWT->bwt->cumulativeFreq[4] + readID;
    ULL absOffset = readID * idx2BWT->readLength + readLen;

    for (;;) {
        unsigned int c = _getBWTvalue(idx2BWT->bwt, l);
        if (c == 4) {
          for (int i = 0, j = readLen - 1; i < j; ++i, --j) {
              seq[i] ^= seq[j];
              seq[j] ^= seq[i];
              seq[i] ^= seq[j];
          }
          seq[readLen] = 0;
          qual[readLen] = '\0';
          //l = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c);
          return readLen;
        } else {
          l = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c);
          seq[readLen] = dnaChars[c];
          qual[readLen] = (bwt->bwtCode[absOffset / CHAR_PER_WORD] >> (absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1))) & 1;
          ++readLen;
          ++absOffset;
        }
    }
}

//Compare the pos-seed part of reads to the assembled consensus.
int plus_extract_and_compare(Idx2BWT * idx2BWT, unsigned long long readID, int start, int end, char *check_seq, int check_len, int startp)
{
    static char dnaChars[] = {'A','C','G','T'};

    BWT *bwt = idx2BWT->bwt;
    ULL l = idx2BWT->bwt->cumulativeFreq[4] + readID;
    ULL absOffset = readID * idx2BWT->readLength + idx2BWT->readLength-1;

    char b;
    unsigned int c, q, i, mismatch=0;
//if(DEBUG) fprintf(stdout, "plus read id %lld, start %d, end %d, startp %d\n", readID, start, end, startp);
    for (i=0;i<Read_length;i++) {
        unsigned int c = _getBWTvalue(idx2BWT->bwt, l);
        if (c == 4 || i > start) 
          return 1;
        q = (bwt->bwtCode[absOffset / CHAR_PER_WORD] >> (absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1))) & 1;
        b = dnaChars[c];

        if(i >= end && startp - i + end < check_len)
        {
//if(DEBUG) fprintf(stdout, "pos %d, cons %c seq %c, qual %d\n", i, check_seq[startp - i + end], b, q);
            if(b != check_seq[startp - i + end] && q == 1)
            {
                if(q == 1 || mismatch >=5)
                {
if(DEBUG)   fprintf(stdout, "Plus i: %d,%c->%c,%d; with startp %d, end %d\n", i, check_seq[startp - i + end], b, q, startp, end);
                    return 0;
                }
                mismatch++;
            }
        }

        l = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c);
        absOffset--;
    }
    return 1;
}


int minus_extract_and_compare(Idx2BWT * idx2BWT, unsigned long long readID, int start, int end, char *check_seq, int check_len, int startp)
{
    static char dnaChars[] = {'A','C','G','T'};

    BWT *bwt = idx2BWT->bwt;
    ULL l = idx2BWT->bwt->cumulativeFreq[4] + readID;
    ULL absOffset = readID * idx2BWT->readLength;

    char b;
    unsigned int c, q, i, mismatch=0;
//if(DEBUG) fprintf(stdout, "minus read id %lld, start %d, end %d, startp %d\n", readID, start, end, startp);fflush(stdout);
    for (i=0;i<Read_length;i++) {
        unsigned int c = _getBWTvalue(idx2BWT->bwt, l);
        if (c == 4) 
          return 1;
        q = (bwt->bwtCode[absOffset / CHAR_PER_WORD] >> (absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1))) & 1;
        b = dnaChars[3-c];

        if(i >= start && startp + i - start < check_len)
        {
//if(DEBUG) fprintf(stdout, "pos %d, cons %c seq %c, qual %d\n", i, check_seq[startp + i - start], b, q);fflush(stdout);
            if(b != check_seq[startp + i - start])
            {
                if(q == 1 || mismatch >= 5)
                {
if(DEBUG)   fprintf(stdout, "Minus i: %d,%c->%c,%d\n", i, check_seq[startp + i -start], b, q);
                    return 0;
                }
                mismatch++;
            }
        }

        l = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c);
        absOffset++;
    }
    return 1;
}


void extractQual(Idx2BWT * idx2bwt, ULL id, int* qual)
{
    BWT * bwt = idx2bwt->bwt;
    int i, q;
    ULL absOffset;

    for(i=0; i< idx2bwt->readLength; i++)
    {
        absOffset = id * idx2bwt->readLength + i;
        qual[i] = (bwt->bwtCode[absOffset / CHAR_PER_WORD] >> (absOffset % CHAR_PER_WORD * BIT_PER_CHAR + (BIT_PER_CHAR - 1))) & 1;
    }
    qual[i]='\0';
}

void get_id_seq(char* seq, int* qual, ULL id, int size)
{
    extractReadQual(idx2BWT[size], id, seq, qual);
}

/*
Attention: index from new 2bwt need to minus 1 when computing the l,r,ll, rr.
*/

void BWTSARangeBackwardFor(Idx2BWT * idx2bwt, const unsigned char c,
                        unsigned long long *saIndexLeft, unsigned long long *saIndexRight) {

    BWT * l_bwt = idx2bwt->bwt;
    BWT * l_rev_bwt = idx2bwt->rev_bwt;

    unsigned long long l = (*saIndexLeft);
    unsigned long long r = (*saIndexRight);
    (*saIndexLeft)  = l_bwt->cumulativeFreq[c] + BWTOccValue(l_bwt, l, c);// + 1;
    (*saIndexRight) = l_bwt->cumulativeFreq[c] + BWTOccValue(l_bwt, r + 1, c)-1;
}

void BWTSARangeBackwardRev(Idx2BWT* idx2bwt, unsigned char rcc, ULL* l, ULL* r, ULL* rev_l, ULL* rev_r)
{
      ULL rr = (*r);
      ULL ll = (*l);
      ULL rev_ll= (*rev_l);
      ULL rev_rr= (*rev_r);

      BWT *l_rev_bwt = idx2bwt->rev_bwt;
      BWT *l_bwt = idx2bwt->bwt;

      ULL occL[ALPHABET_SIZE];
      ULL occR[ALPHABET_SIZE];
      ULL oCount[ALPHABET_SIZE];

      BWTAllOccValue(l_rev_bwt,rev_ll, occL);
      BWTAllOccValue(l_rev_bwt,rev_rr + 1,occR);

      oCount[ALPHABET_SIZE-1]=0;
      int k;
      for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+occR[k+1]-occL[k+1];
      }

      rev_ll = l_bwt->cumulativeFreq[rcc] + occL[rcc] + 1;
      rev_rr = l_bwt->cumulativeFreq[rcc] + occR[rcc];
      rr = rr - oCount[rcc];
      ll = rr - (rev_rr - rev_ll);

      rev_ll--; //new.
      rev_rr--; //new.
      (*r) = rr; 
      (*l) = ll; 
      (*rev_l) = rev_ll; 
      (*rev_r) = rev_rr;
}

void initial_bases(Base* seeds, BWT **bwts, unsigned char c, unsigned char rc, int count)
{
    int i;
    for(i=0; i< count; i++)
    {
        seeds[i].saL = bwts[i]->cumulativeFreq[c];
        seeds[i].saR = bwts[i]->cumulativeFreq[c+1]-1;
        seeds[i].saLL = bwts[i]->cumulativeFreq[rc]; 
        seeds[i].saRR = bwts[i]->cumulativeFreq[rc+1]-1;
        seeds[i].saLLrev = seeds[i].saLL;
        seeds[i].saRRrev = seeds[i].saRR;
        seeds[i].plus = 0; seeds[i].minus = 0;
        seeds[i].plus = seeds[i].saR >= seeds[i].saL ? seeds[i].saR - seeds[i].saL + 1 : 0;
        seeds[i].minus = seeds[i].saRR >= seeds[i].saLL ? seeds[i].saRR - seeds[i].saLL + 1 : 0;
    }
}

unsigned long long update_bases(Base* seeds, BWT **bwts, BWT **rev_bwts, unsigned char c, unsigned char rc, int count, unsigned long long *occL, unsigned long long *occR, unsigned long long *oCount)
{
    int i;
    unsigned long long depth=0;
    for(i=0; i< count; i++)
    {
        if(seeds[i].plus > 0)
        {
            seeds[i].saL = bwts[i]->cumulativeFreq[c] + BWTOccValue(bwts[i], seeds[i].saL, c) ; //remove +1.
            seeds[i].saR = bwts[i]->cumulativeFreq[c] + BWTOccValue(bwts[i], seeds[i].saR + 1, c)-1; //add -1.
            seeds[i].plus = seeds[i].saR >= seeds[i].saL ? seeds[i].saR - seeds[i].saL + 1 : 0;
        }
        if(seeds[i].minus > 0)
        {
            BWTAllOccValue(rev_bwts[i],seeds[i].saLLrev,occL);
            BWTAllOccValue(rev_bwts[i],seeds[i].saRRrev + 1,occR);

            oCount[ALPHABET_SIZE-1]=0;
            int k;
            for (k=ALPHABET_SIZE-2;k>=0;k--) {
                oCount[k]=oCount[k+1]+occR[k+1]-occL[k+1];
            }

            seeds[i].saLLrev = bwts[i]->cumulativeFreq[rc] + occL[rc] + 1;
            seeds[i].saRRrev = bwts[i]->cumulativeFreq[rc] + occR[rc];
            seeds[i].saRR = seeds[i].saRR - oCount[rc];
            seeds[i].saLL = seeds[i].saRR - (seeds[i].saRRrev - seeds[i].saLLrev);

            seeds[i].saLLrev--; //new.
            seeds[i].saRRrev--; //new.

            seeds[i].minus = seeds[i].saRR >= seeds[i].saLL ? seeds[i].saRR - seeds[i].saLL + 1 : 0;
        }
        seeds[i].depth = seeds[i].plus + seeds[i].minus;

        if(FLAGS[i] == 1)
            depth += seeds[i].depth;
    }
    return depth;
}

/*Extract SA range for substring of one sequence.*/
void set_SAranges(const char* seq, int start, int end, Base* seed_tmp, unsigned long long* depth)
{
    unsigned long long occL[ALPHABET_SIZE];
    unsigned long long occR[ALPHABET_SIZE];
    unsigned long long oCount[ALPHABET_SIZE];
    unsigned char c, rcc;
    int p = end; 
    unsigned long long curSACount=0;
    c = charMap[ seq[p] ]; rcc = 3 - c; //in seq, there is only atgc.

    initial_bases(seed_tmp, bwts, c, rcc, bwtCount);
    p--;
    while(p >= start)
    {
        c = charMap[ seq[p] ]; rcc = 3 - c;
        curSACount = update_bases(seed_tmp, bwts, rev_bwts, c, rcc, bwtCount, occL, occR, oCount);
        p--;
    }
    *depth = curSACount;
}

/*Check bubbles structures */
int dual_bubble_check(const char* seq, int initial_depth, int threadId)
{
    unsigned long long occL[ALPHABET_SIZE];
    unsigned long long occR[ALPHABET_SIZE];
    unsigned long long oCount[ALPHABET_SIZE];
    unsigned char c1, rcc1;
    int p = strlen(seq)-1, end = strlen(seq)-1, dis = 1;
    c1 = charMap[ seq[p] ]; rcc1 = 3 - c1;

    initial_bases(sbase[threadId][0], bwts, c1, rcc1, bwtCount);
    p--; dis++;

    double higher_bound = initial_depth * Upper_bound, current_base_depth;

if(DEBUG)   fprintf(stdout, "seq: %s\t initial_depth: %d, higher_bound: %.1f\n", seq, initial_depth, higher_bound);

    unsigned long long curSACount1=0, curSACount2=0;
    while(p >= 0)
    {
        c1 = charMap[ seq[p] ]; rcc1 = 3 - c1;

        //set_base(sbase[threadId][2], sbase[threadId][0]);
        curSACount1 = update_bases(sbase[threadId][0], bwts, rev_bwts, c1, rcc1, bwtCount, occL, occR, oCount);
        current_base_depth = double(curSACount1 * Read_length)/double(Read_length-dis+1);

if(DEBUG) fprintf(stdout, "Pos%d,dis%d:%lld, %.1f, %.1f\n", p, dis, curSACount1, current_base_depth, current_base_depth*Upper_bound);

        if(current_base_depth < Expect_depth * Upper_bound && current_base_depth < higher_bound && initial_depth < current_base_depth * Upper_bound)
            return 1;
        p--; dis++;
    }
    return 0;
}

int get_other_depth(int c, int p, int threadId, unsigned long long *occL, unsigned long long *occR, unsigned long long *oCount)
{
    int i, depth, max_depth=0;
    for(i=0; i<4; i++)
    {
        if(i==c)
            continue;
        set_base(sbase[threadId][p], sbase[threadId][p+2]);
        depth = update_bases(sbase[threadId][p], bwts, rev_bwts, i, 3-i, bwtCount, occL, occR, oCount);
        if(max_depth < depth)
            max_depth = depth;
    }
    return max_depth;
}

/* Repeat and sequencing error analysis*/
int dual_seq_check(FILE* logs, const char* seq1, const char* seq2, int tlen, double base_depth, int d2, double diff, int threadId)
{
    unsigned long long occL[ALPHABET_SIZE];
    unsigned long long occR[ALPHABET_SIZE];
    unsigned long long oCount[ALPHABET_SIZE];
    unsigned char c1, rcc1, c2, rcc2;
    int p = strlen(seq1)-1, end = strlen(seq1)-1, dis = 1, hete = 0, has_sharp = 0, has_sharp2= 0;
    c1 = charMap[ seq1[p] ]; rcc1 = 3 - c1;
    c2 = charMap[ seq2[p] ]; rcc2 = 3 - c2;

if(DEBUG)   fprintf(logs, "seq1: %s\nseq2: %s tlen: %d\n", seq1, seq2, tlen);

    initial_bases(sbase[threadId][0], bwts, c1, rcc1, bwtCount);
    initial_bases(sbase[threadId][1], bwts, c2, rcc2, bwtCount);
    p--; dis++;

    unsigned long long curSACount1=0, curSACount2=0, first_s1_depth=0, first_s2_depth=0, prevSACount1 = 0, prevSACount2 = 0;
    double current_expect_depth;
    int depth_conse_length = Read_length - int(double(Read_length * Consensus_depth) / double(Expect_depth)) + 1;
    int depth_error_length = Read_length - int(double(Read_length * Error_depth) / double(Expect_depth)) + 1;

    double current_base_depth = base_depth;
    if(current_base_depth > Expect_depth)
        current_base_depth = Expect_depth;

    while(p >=0)
    {
        c1 = charMap[ seq1[p] ]; rcc1 = 3 - c1;
        c2 = charMap[ seq2[p] ]; rcc2 = 3 - c2;

        current_expect_depth = double(Read_length - dis + 1)/double(Read_length)* current_base_depth - 1;
        if(current_expect_depth < Consensus_depth)
            current_expect_depth = Consensus_depth;

        set_base(sbase[threadId][2], sbase[threadId][0]);
        set_base(sbase[threadId][3], sbase[threadId][1]);
        curSACount1 = update_bases(sbase[threadId][0], bwts, rev_bwts, c1, rcc1, bwtCount, occL, occR, oCount);
        curSACount2 = update_bases(sbase[threadId][1], bwts, rev_bwts, c2, rcc2, bwtCount, occL, occR, oCount);

if(DEBUG)  fprintf(logs, "Pos%d,dis%d:%lld&%lld,%.1f\n", p, dis, curSACount1, curSACount2, current_expect_depth*Upper_bound);

        //to too long check.
        if(dis >= depth_conse_length || dis >= depth_error_length)
            return 0;

        if(dis > tlen && current_expect_depth >= 4 * Error_depth && curSACount1 + curSACount2 >= 4 * Error_depth)
        {
            if(curSACount1 + curSACount2 < current_expect_depth*Upper_bound )
                hete = 1;
        }

        //check for low depth for path1.
        if(curSACount1 <= 3*Error_depth)
        {
            //for sequencing error.
            if(first_s1_depth == 0 && curSACount1 == 1
               && curSACount2 > Error_depth
               && current_expect_depth > 2 * Error_depth
               && (diff < MAXDIFRATIO || diff * d2 <= MAXDIFNUM|| curSACount1 < current_expect_depth * 0.01)
              ) 
            {
if(DEBUG)       fprintf(logs, "Sequencing error tips 1: depth1 %d depth2 %d expect depth %.1f diff %.1f\n", curSACount1, curSACount2, current_expect_depth, diff);
                return 2;
            }

            //for sequencing error2.
            if(curSACount1 <= 3*Error_depth && curSACount2 > 2* Error_depth && current_expect_depth > 2*Error_depth
              && (diff < MAXDIFRATIO || diff * d2 < MAXDIFNUM || (curSACount1 != 0 && curSACount1 < current_expect_depth * 0.05))
              && (curSACount1 < double(curSACount2) / double(Consensus_depth) && curSACount1 < current_expect_depth /double(2 * Consensus_depth))
              && (curSACount1 >= Consensus_depth  
                  || (curSACount1 > 1 && curSACount2 >= 4* Consensus_depth && current_expect_depth >= 5 * Consensus_depth && !(curSACount2 > 2 * current_expect_depth && current_expect_depth < 6 * Consensus_depth)) 
                  || (curSACount1 <= 1 && curSACount2 > 3* Consensus_depth && current_expect_depth >= 4 * Consensus_depth)
                 )
              )
            {
if(DEBUG)       fprintf(logs, "Sequencing error tips 2: depth1 %d depth2 %d expect depth %.1f diff %.1f\n", curSACount1, curSACount2, current_expect_depth, diff);
                return 2;
            }

            //sharp decrease.
            if(prevSACount1 - curSACount1 >= Error_depth && curSACount1 < Error_depth
               && (((curSACount1 == 0 && curSACount2 > 1 && current_expect_depth > Error_depth )//for low depth.
                    && (current_expect_depth >= 2* Consensus_depth 
                        || (curSACount2 <= current_expect_depth && current_expect_depth >= Consensus_depth) 
                        || (prevSACount1 >= 2*Consensus_depth && current_expect_depth*Upper_bound >= 2*Consensus_depth)
                        || (prevSACount1 >= 2*Consensus_depth && curSACount2 <= current_expect_depth*Upper_bound && curSACount2 >= Consensus_depth)
                       ) //for high depth.
                   )
                  || ((curSACount1 != 0 && curSACount2 > 2* Error_depth && current_expect_depth > 2*Error_depth) //left depth is caused by sequencing error.
                    && (diff < MAXDIFRATIO || diff * d2 < MAXDIFNUM|| curSACount1 < current_expect_depth * 0.05 || (curSACount1 < current_expect_depth * 0.1 && prevSACount1 - curSACount1 >= 3*Consensus_depth))
                    && (curSACount1 < double(curSACount2) / double(Consensus_depth) && curSACount1 < current_expect_depth /double(2*Consensus_depth))
                    && (curSACount1 >= Consensus_depth
                        || (curSACount1 == 1 && curSACount2 >= 2* Consensus_depth && current_expect_depth*Upper_bound > 3* Consensus_depth)
                        || (curSACount1 > 1 && prevSACount2 - curSACount2 < Error_depth)
                        || (curSACount1 > 1 && curSACount2 >=4 * Consensus_depth && current_expect_depth >= 4 * Consensus_depth && (diff < MAXDIFRATIO/2 || diff * d2 < MAXDIFNUM/2+1))
                       )
                   )
                 )
              )
            {
                int other_depth = get_other_depth(c1, 0, threadId, occL, occR, oCount);
                if(other_depth >= 2*Error_depth 
                   || (other_depth >= Error_depth && curSACount2 > Error_depth) 
                   || (other_depth >= Error_depth && curSACount1 == 0 && current_expect_depth >= 2* Consensus_depth)
                   || (other_depth >= curSACount2 && curSACount1 == 0 && current_expect_depth > curSACount2)
                )
                {
if(DEBUG)           fprintf(logs, "Butterfly structure 1: prev depth1 %d depth1 %d depth2 %d otherdepth %d expect depth %.1f diff %.1f\n", prevSACount1, curSACount1, curSACount2, other_depth, current_expect_depth, diff);
		    if(dis <= tlen)
                        return 2;
                    return 1;
                }
if(DEBUG)       fprintf(logs, "other depth %d\n", other_depth);
            }

            if(prevSACount1 - curSACount1 >= Error_depth && curSACount1 <= 2*Error_depth )//&& curSACount1 < current_expect_depth /double(2*Consensus_depth))
                has_sharp = 1;

            if(curSACount1 <= 1 && prevSACount1 <= Error_depth && has_sharp == 1 && curSACount2 >= 3* Consensus_depth && current_expect_depth > 4 * Consensus_depth
                //&& ((curSACount1 == 0) || (diff < 0.5 || diff * d2 < 10))
              )
            {
if(DEBUG)       fprintf(logs, "Butterfly structure 3: prev depth1 %d depth1 %d depth2 %d expect depth %.1f diff %.1f\n", prevSACount1, curSACount1, curSACount2, current_expect_depth, diff);
                return 1;
            }

        }

        //check errors in cons path.
        if(curSACount2 <= Error_depth)
        {
           if(prevSACount2 - curSACount2 >= Error_depth && curSACount2 < Error_depth //&& (current_expect_depth > 3 * Consensus_depth || (curSACount2 + curSACount1) < current_expect_depth * Upper_bound)
               && ( ((curSACount2 == 0 && curSACount1 > 1 && current_expect_depth > Error_depth )//for high depth.
                     && (current_expect_depth >= 2* Consensus_depth 
                         || (curSACount1 <= current_expect_depth && current_expect_depth >= Consensus_depth) 
                         || (prevSACount2 >= 2*Consensus_depth && current_expect_depth*Upper_bound >= 2*Consensus_depth)
                         || (prevSACount2 >= 2*Consensus_depth && curSACount1 <= current_expect_depth*Upper_bound && curSACount1 >= Consensus_depth)
                        ) //for low depth.
                    )
                  || ((curSACount2 != 0 && curSACount1 > 2* Error_depth && current_expect_depth > 2*Error_depth) //left depth is caused by sequencing error
                     && (diff < MAXDIFRATIO || diff * d2 < MAXDIFNUM || curSACount2 < current_expect_depth * 0.05 || (curSACount2 < current_expect_depth * 0.1 && prevSACount1 - curSACount1 >= 3*Consensus_depth))
                     && (curSACount2 < double(curSACount1) / double(Consensus_depth) && curSACount2 < current_expect_depth /double(2*Consensus_depth))
                     && (curSACount2 >= Consensus_depth 
                         || (curSACount2 == 1 && curSACount1 >= 2* Consensus_depth && current_expect_depth*Upper_bound > 3* Consensus_depth)
                         || (curSACount2 > 1 && prevSACount1 - curSACount1 < Error_depth)
                         || (curSACount2 > 1 && curSACount1 >= 4* Consensus_depth && current_expect_depth >= 4 * Consensus_depth && (diff < MAXDIFRATIO/2 || diff * d2 < MAXDIFNUM/2+1))
                        )
                     )
                  )
             )
           {
               int other_depth = get_other_depth(c2, 1, threadId, occL, occR, oCount);
               if(other_depth >= 2*Error_depth 
                  || (other_depth >= Error_depth && curSACount1 >Error_depth) 
                  || (other_depth >= Error_depth && curSACount2 == 0 && current_expect_depth >= 2* Consensus_depth)
                  || (other_depth >= curSACount1 && curSACount2 == 0 && current_expect_depth > curSACount1)
                )
               {
if(DEBUG)          fprintf(logs, "Butterfly structure 2: prev depth2 %d depth1 %d depth2 %d otherdepth %d expect depth %.1f diff %.1f\n", prevSACount2, curSACount1, curSACount2, other_depth, current_expect_depth, diff);
                   return -1;
               }
if(DEBUG)      fprintf(logs, "other depth %d\n", other_depth);
           }
        }

        if(prevSACount2 - curSACount2 >= Error_depth && curSACount2 <= 2*Error_depth )//&& curSACount2 < current_expect_depth /double(2*Consensus_depth))
                has_sharp2 = 1;

        if(curSACount2 <= 1 && prevSACount2 <= Error_depth && has_sharp2 == 1 && curSACount1 >= 3* Consensus_depth && current_expect_depth > 4 * Consensus_depth 
           //&& ((curSACount2 == 0) || (diff < 0.5 || diff * d2 < 10))
          )
        {
if(DEBUG)       fprintf(logs, "Butterfly structure 4: prev depth1 %d depth1 %d depth2 %d expect depth %.1f diff %.1f\n", prevSACount2, curSACount2, curSACount1, current_expect_depth, diff);
                return -1;
        }

        if(curSACount2 * curSACount1 == 0)
        {
            if(has_sharp2 == 1 && has_sharp == 1 && curSACount2 ==0 && curSACount1 == 0)
                return 5;
            if(hete == 1 && curSACount2 + curSACount1 < 2 * Error_depth)
            {
                if(curSACount2 == 0
                   && (curSACount1 >= Error_depth
                       || (curSACount1 > 1 && prevSACount2 >= Error_depth )
                      )
                   )
                    return 4;
                return 3;
            }
            return 0;
        }
        if(curSACount1 <= Error_depth && curSACount2 <= Error_depth)
        {
            if(hete == 1)
                return 3;
            return 0;
        }
        if(dis == MINOVERLAP )
        {
            first_s1_depth = curSACount1;
            first_s2_depth = curSACount2;
            //remove tips.
        }
        prevSACount1 = curSACount1;
        prevSACount2 = curSACount2;
        p--; dis++;
    }
    if(p == -1 && hete == 1)
        return 3;
    return 0;
}


//Obtain seeds within extension.
int extend_backwards_rep(const char* seq, int start, int end, Base* seed_tmp, unsigned long long* bdepth, int *rep_start)
{
    unsigned long long occL[ALPHABET_SIZE];
    unsigned long long occR[ALPHABET_SIZE];
    unsigned long long oCount[ALPHABET_SIZE];
    unsigned char c, rcc;
    int p = end; //begin and 0 is 0 oriented.
    unsigned long long curSACount=0, prevSACount=0;
    c = charMap[ seq[p] ]; rcc = 3 - c; //in seq, there is only atgc.

    initial_bases(seed_tmp, bwts, c, rcc, bwtCount);
    p--;

    float edepth = Upper_bound * (float)Expect_depth, nedepth;
    int is_unique = 0, use_rep_start = 0, is_sharp_decrease=0, last_sharp_p=0;
    int sharp_start = -1, high_start = -1, sharp_depth = 0, high_depth = 0;
    while(p >= start)
    {
        c = charMap[ seq[p] ]; rcc = 3 - c;
        curSACount = update_bases(seed_tmp, bwts, rev_bwts, c, rcc, bwtCount, occL, occR, oCount);
        int dis = end-p+1;
   
        if(curSACount <= Error_depth)
        {
            break;
        }
        nedepth = edepth*(float)((Read_length-dis+1))/double(Read_length) ;

        if(dis >= MINOVERLAP/2
           && ((dis >= MINOVERLAP && curSACount < prevSACount/10 && curSACount <= nedepth/3 && curSACount > 4 * Error_depth)
              )
          )
        {
            break;
        }

        if( ((prevSACount - curSACount > 3*Error_depth)
             ||(prevSACount - curSACount >= nedepth/2 && prevSACount - curSACount > Error_depth)
            ) 
            &&((nedepth > 3 * Error_depth)  
              || (last_sharp_p - p > 3 * Consensus_depth && curSACount > 2 * Error_depth)
              )
          )
        {
            is_sharp_decrease = 1;
            last_sharp_p = p;
        }else
            is_sharp_decrease = 0;

if(DEBUG )  fprintf(stdout, "dis %d depth %lld nedepth %.1f sharp %d higher %d\n", dis, curSACount,nedepth, is_sharp_decrease, curSACount >= nedepth/Upper_bound*1.8);	

        if(dis < depth_6_length && curSACount > 3*Error_depth && curSACount < MAXDEPTH 
           && ((nedepth > 4 * Error_depth && curSACount > 4 * Error_depth && curSACount >= nedepth/Upper_bound*1.8)
               || (nedepth > 2 * Error_depth && curSACount >= 2 * Error_depth && is_sharp_decrease ==1)
              )
          )
        {    *rep_start = p;
        }

        if(dis >= MINOVERLAP && curSACount <= nedepth)
        {
           if(curSACount <= 3*Error_depth && (last_sharp_p - p> 3 * Consensus_depth || (prevSACount > 4 * curSACount && prevSACount < 4 *nedepth) || (prevSACount > 2 * curSACount) && prevSACount < 2 * nedepth))
               break;

           if(prevSACount > curSACount * 10 && curSACount < 4 * Error_depth && nedepth < 4 * Error_depth)
               break;

           if(curSACount <= 2 * Error_depth && prevSACount - curSACount >= Error_depth && prevSACount - curSACount >= curSACount) 
               break;

           is_unique=1;
           break;
        }

        if(nedepth < 2 * Error_depth || curSACount < 2*Error_depth)
        {
            break;
        }
        prevSACount = curSACount;
        p--;

    }

    *bdepth = curSACount;
    if(is_unique == 0)
    {
        return -2;
    }
    return p;
}

//Obtain seeds to initial extension.
int extend_backwards(const char* seq, int start, int end, Base* seed_tmp, unsigned long long* bdepth, int *break_start)
{

    unsigned long long occL[ALPHABET_SIZE];
    unsigned long long occR[ALPHABET_SIZE];
    unsigned long long oCount[ALPHABET_SIZE];
    unsigned char c, rcc;
    int p = end; //begin and 0 is 0 oriented.
    unsigned long long curSACount=0, prevSACount=0;
    c = charMap[ seq[p] ]; rcc = 3 - c; //in seq, there is only atgc.

    initial_bases(seed_tmp, bwts, c, rcc, bwtCount);
    p--;

    float edepth = Upper_bound * (float)Expect_depth, nedepth;
    int is_unique = 0;
    while(p >= start)
    {
        c = charMap[ seq[p] ]; rcc = 3 - c;
        curSACount = update_bases(seed_tmp, bwts, rev_bwts, c, rcc, bwtCount, occL, occR, oCount);
        int dis = end-p+1;

        if(curSACount <= 2*Error_depth)
        {
            *break_start = p;
            break;
        }

        nedepth = edepth*(double)((Read_length-dis+1))/double(Read_length) ;
if(DEBUG ) fprintf(stdout, "dis %d depth %lld nedepth %.1f\n", dis, curSACount,nedepth);

        if(dis >= MINOVERLAP && curSACount <= nedepth/3 && curSACount > 4 * Error_depth)
        {
            *break_start = -1;
            break;
        }

        if(dis >= MINOVERLAP &&  (curSACount <= int(nedepth) && curSACount >= int(nedepth)/2) )
        {
           is_unique=1;
           *break_start = p;
           break;
        }

        if(dis > depth_6_length)
        {
            *break_start = p;
            break;
        }
        prevSACount = curSACount;
        p--;
    }

    if(p < start)
        *break_start = 0;

    *bdepth = curSACount;
    if(is_unique == 0)
        return -1;
    return p;
}

int is_allhigh_qual(int* qual, int s, int e)
{
    int i;
    for(i=s; i<= e; i++)
    {
        if(qual[i] == 0)
            return 0;
    }
    return 1;
}

void set_base(Base* to, Base* from)
{
    int i;
    for(i=0; i< bwtCount; i++)
    {
        to[i] = from[i];
    }
}

void print_b(Base* b)
{
    int i;
    for(i=0; i< bwtCount; i++)
    {
        fprintf(stderr, "i: %d l %lld r %lld ll %lld rr %lld plus %d minus %d\n", i, b->saL, b->saR, b->saLL, b->saRR, b->plus, b->minus);
    }
}


int initial_contig_seed(FILE *logs, ULL id, int size, int threadId)
{
    char seq[Read_length+1];
    int qual[Read_length+1];
    
    get_id_seq(seq, qual, id, size);

    int start=0, end = strlen(seq)-1;
    int flag = -1;
    unsigned long long depth=0, first_repeat_depth=0;
    int first_repeat_start=-1, first_repeat_end, break_start = 0;

    while(flag == -1 && end >= start + MINOVERLAP-1)
    {
        while(qual[end] == 0 && end > start + MINOVERLAP-1)
            end --;
        if(end < start + MINOVERLAP-1)
        {
            flag == -1;
            break;
        }
if(DEBUG)    fprintf(stdout, "Ending at %d\n", end);
        flag = extend_backwards(seq, start, end, g_iters[threadId][0].b, &(depth), &break_start);

	   //pure repetitive.
        if(flag == -1 && break_start != -1 && depth > 2*Error_depth && first_repeat_start == -1 && depth < MAXDEPTH)
    	{
    		first_repeat_start = break_start;
    		first_repeat_end = end;
    		if(is_allhigh_qual(qual, first_repeat_start, first_repeat_end) == 1)
    		{
    			first_repeat_depth = depth;
    			set_base(tmpbase[threadId] , g_iters[threadId][0].b);
    		}else{
    			first_repeat_start = -1;
    		}
    	}
if(DEBUG)    fprintf(stdout, "flag %d first_repeat_start %d first_repeat_end %d\n", flag, first_repeat_start, first_repeat_end);
        if(flag != -1 && is_allhigh_qual(qual, flag, end) == 1)
        {
            break;
        }
        end -= MINOVERLAP;
        flag = -1;
    }

    //pure repetitive.
    if(flag == -1 && first_repeat_start != -1)
    {
        flag = first_repeat_start;
    	end = first_repeat_end;
    	depth = first_repeat_depth;
    	set_base(g_iters[threadId][0].b, tmpbase[threadId]);
    }

    if(flag == -1)
    {
        g_iters[threadId][0].depth = 0;
        return 0;
    }

    int i=0, len = end-flag + 1;
    for(i=0; i< len; i++)
    {
        ctgs[threadId][ ctg_num[threadId] ].seq[i] = seq[flag + (len -1 - i)]; //the sequence of contig is only reversed.
        g_iters[threadId][0].seq[i]  = seq[flag + i];
    }

    ctgs[threadId][ ctg_num[threadId] ].seq[i]='\0';
    g_iters[threadId][0].seq[i]='\0';

    ctgs[threadId][ ctg_num[threadId] ].len = i;

    g_iters[threadId][0].slen = i;
    if(depth > MAXDEPTH)
        g_iters[threadId][0].depth = MAXDEPTH;
    else
        g_iters[threadId][0].depth = depth;

    return 1;
}

int count_lib(char* lib_file)
{
    FILE* fp=fopen(lib_file, "r");
    char line[256];
    int c=0;

    while(fgets ( line, sizeof ( line ), fp ) != NULL)
    {
        if(line[0] == '#')
            continue;
        c++;
    }
    fclose(fp);
    return c;
}

int load_bwts(char* lib_file)
{
    char prefix[256], line[256];
    int insert, sd, flag, i=0;

    int libCount=count_lib(lib_file);
    fprintf(stderr, "lib number is: %d\n", libCount);
    idx2BWT = (Idx2BWT**)calloc(libCount, sizeof(Idx2BWT*));

    InsertSize = (int*) calloc(libCount, sizeof(int));
    SD = (int*) calloc(libCount, sizeof(int));
    FLAGS = (int*) calloc(libCount, sizeof(int));

    FILE* fp=fopen(lib_file, "r");
    while(fgets ( line, sizeof ( line ), fp ) != NULL)
    {
        if(line[0] == '#')
            continue;
        sscanf(line, "%s %d %d %d", prefix, &insert, &sd, &flag);
        fprintf(stderr, "Loading index: %s...\n", prefix);
        idx2BWT[i] = BWTLoad2BWTLite(prefix);
        InsertSize[i] = insert;
        SD[i] = sd;
        FLAGS[i] = flag;
        i++;
    }
    fprintf(stderr, "Loading %d libraries.\n", i);
    fclose(fp);
    return i;
} 

void free_bwts(int num)
{
    int i;
    for(i=0; i< num; i++)
        BWTFree2BWT(idx2BWT[i]);
    free(idx2BWT);

    free(InsertSize);
    free(SD);
    free(FLAGS);
}

