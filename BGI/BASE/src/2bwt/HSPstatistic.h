/*

   HSPstatistic.h     

   Copyright (C) 2014 The University of Hong Kong.

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
/* ==========================================================================
 *
 * HSPstatistic.h
 *
 * author: Wong Swee Seong
 * date: 30-03-2006
 *
 * gapped and ungapped extension for string matching
 * 
 * ==========================================================================*/


#ifndef _HSP_statistic_
#define _HSP_statistic_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define int4 int32_t
#define uint4 uint32_t
#define int8 int64_t
#define uint8 uint64_t
#define int2 int16_t
#define uint2 uint16_t
#define DNA_CHAR_SIZE 18

typedef struct stat__score_block
{
  int match, mismatch, gapOpen, gapExt;
} STAT_scoreBlock;

typedef struct stat__Xdropoffs
{
  double ungapXdropoff, gapXdropoff, gapXdropoffFinal;
} STAT_Xdropoffs;

void initializeHSPstatistic(uint8 databaseSize, uint8 numOfSequences, uint8 minSeqLength, double *dbCharProb,
                            uint8 querySize, uint8 numOfContext, double *queryCharProb, int scoringMatrix[DNA_CHAR_SIZE][DNA_CHAR_SIZE]);
void printHSPstatistic(FILE * outFile);
int4 calcUngapCutoffScore();
int4 calcGapCutoffScore();
int4 getUngapXdropoff();
int4 getGapXdropoff();
int4 getGapXdropoffFinal();

int4 stat_ungapNormalized2nominal(double normalizedScore);
double stat_ungapNominal2normalized(int4 nominalScore);
double stat_ungapCalcEvalue(double normalizedScore);
int4 stat_ungapEvalue2nominal(double evalue);
double stat_gapNominal2normalized(int4 nominalScore);
int4 stat_gapNormalized2nominal(double normalizedScore);
double stat_gapCalcEvalue(double normalizedScore);
int4 stat_gapEvalue2nominal(double evalue);
 
#ifdef __cplusplus
}
#endif

#endif
