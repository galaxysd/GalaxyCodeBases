/*

   karlin.h      

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

#ifndef _karlin_
#define _karlin_

/*
wss:
The following definitions are added to remove dependency on blast.h
*/
#include <stdint.h>
#define int4 int32_t
#define uint4 uint32_t
#define int8 int64_t
#define uint8 uint64_t
#define int2 int16_t
#define uint2 uint16_t

extern double BlastKarlin_lambda;
extern double BlastKarlin_K;
extern double BlastKarlin_H;

void BlastKarlinBlkCalc(double* scoreProbabilities, int4 min, int4 max);

int4 BlastComputeLengthAdjustment(double K,
                             double logK,
                             double alpha_d_lambda,
                             double beta,
                             int4 query_length,
                             uint4 db_length,
                             int4 db_num_seqs,
                             int4 *length_adjustment);
#endif

