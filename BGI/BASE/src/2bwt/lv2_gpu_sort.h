/*

   lv2_gpu_sort.h		sort lev2 with GPU

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
    Change : Generate this file by fishing from Chi Man LIU's code.

*/
#ifndef LV2_GPU_SORT_H
#define LV2_GPU_SORT_H

#include <stdint.h>
#include "misc.h"
// single thread
#define THREADS_PER_BLOCK 256 
void lv2_gpu_sort(uint32_t *lv2_substrings, uint32_t *permutation, int words_per_substring, int64_t width, int64_t lv2_num_items);
size_t get_gpu_mem();

#endif

