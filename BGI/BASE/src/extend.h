#ifndef __EXTEND__H_
#define __EXTEND__H_
/*

   extend.h        extend in one direction of BASE.

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

    Date   : 28th Jan 2016
    Author : Binghang Liu, bhliu@cs.hku.hk
*/

#include "bwt_api.h"

int backward_extending_bwtpe(FILE* logs, int is_first, int threadId);

void free_thread_tmpRead(int threadId, int flag);
void clean_trimed_reads(int ctg_len, int read_num, int threadId);

ULL get_pair_id(ULL id);
int is_unique(int depth, int len);
int is_repeat(int depth, int len);

void reverse_com(char* from, char* to, int len);
void str_copy(char* seq, char* from, int start, int len);

#endif
