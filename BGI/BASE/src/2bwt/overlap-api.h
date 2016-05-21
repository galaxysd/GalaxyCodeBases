#ifndef __OVERLAP_API_H
#define __OVERLAP_API_H

/*

   overlap-api.h        Basic BWT operations in BASE.

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

#include "2BWT-Interface.h"

typedef unsigned long long ULL;

/*
 * ExactRange represents an SA range (= a bunch of reads) found during an exact match
 * overlap search.
 */

struct ReadInf {
  int read_id; 
  int strand; // 1 means positive, 2 means negative.
  int pos; //position in read.
};

struct ExactRange {
  int overlapLength;
  int strand; // 1 means positive, 2 means negative
  ULL saL, saR; // left and right boundaries of SA range
  ULL rsaL, rsaR; // left and right boundaries of reverse SA range
  ExactRange(int overlapLength, int strand, ULL saL, ULL saR, ULL rsaL, ULL rsaR) : // constructor
    overlapLength(overlapLength), strand(strand), saL(saL), saR(saR), rsaL(rsaL), rsaR(rsaR) {}
  ExactRange() : overlapLength(0), strand(0), saL(0), saR(0), rsaL(0), rsaR(0) {} // default constructor
};

int extractReadInf(Idx2BWT * idx2BWT, ExactRange range, ReadInf* ri, int outputLimit);
int extractRead(Idx2BWT * idx2BWT, unsigned long long readID, char* seq, char* qual = NULL);

#endif

