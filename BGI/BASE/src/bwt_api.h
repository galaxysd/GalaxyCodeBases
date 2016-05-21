#ifndef __BWT_API__H_
#define __BWT_API__H_

/*

   bwt_api.cpp        BWT operations in BASE.

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
   Author : Binghang LIU
   Change : write the codes in the process to develop BASE.
*/

#include<linux/types.h>
#include <math.h>
#include "./2bwt/overlap-api.h"
#include "global.h"

int extract_and_compare(Idx2BWT * idx2BWT, unsigned long long readID, int cut, char *check_seq, int check_len, int startp);
int plus_extract_and_compare(Idx2BWT * idx2BWT, unsigned long long readID, int start, int end, char *check_seq, int check_len, int startp);
int minus_extract_and_compare(Idx2BWT * idx2BWT, unsigned long long readID, int start, int end, char *check_seq, int check_len, int startp);

int extend_backwards(const char* seq, int start, int end, Base* seed_tmp, unsigned long long* depth, int *break_start);
int extend_backwards_rep(const char* seq, int start, int end, Base* seed_tmp, unsigned long long* bdepth, int *rep_start);
int initial_contig_seed(FILE* logs,ULL id, int size, int threadId);
void get_id_seq(char* seq, int* qual, ULL id, int size);

void set_SAranges(const char* seq, int start, int end, Base* seed_tmp, unsigned long long* depth);
void set_SArange(int size, const char* seq, int start, int end, Base* seed_tmp);

int dual_seq_check(FILE* logs, const char* seq1, const char* seq2, int tlen, double base_depth, int d2, double diff, int threadId);
int dual_bubble_check(const char* seq, int initial_depth, int threadId);

void BWTSARangeBackwardRev(Idx2BWT* idx2bwt, unsigned char rcc, ULL* l, ULL* r, ULL* rev_l, ULL* rev_r);
void BWTSARangeBackwardFor(Idx2BWT * idx2bwt, const unsigned char c,
                        unsigned long long *saIndexLeft, unsigned long long *saIndexRight);

int check_and_set_bwt_sparse_thread_and_used(ULL id, int size, int threadId);
int clean_use_bwt_sparse(ULL id, int size);
int is_used_bwt_sparse(ULL id, int size);
int get_threadId_bwt_sparse(ULL id, int size);
int clean_current_bwt_sparse(ULL id, int size);

int load_bwts(char* lib_file);
void free_bwts(int num);

void set_base(Base* to, Base* from);

#endif

