/*
     lv2_cpu_sort.h
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
    Author : Binghang LIU, Dinghua Li
    Change : Generate this file by fishing from Chi Man LIU's code.
*/

#include <parallel/algorithm>
#include <assert.h>
#include <stdint.h>

struct CompareHigh32Bits {
    bool operator() (uint64_t a, uint64_t b) {
        return (a >> 32) < (b >> 32);
    }
};

void lv2_cpu_sort(uint32_t *lv2_substrings, uint32_t *permutation, uint64_t *cpu_sort_space, int words_per_substring, int64_t width, int64_t lv2_num_items) {
#pragma omp parallel for
    for (uint32_t i = 0; i < lv2_num_items; ++i) {
        cpu_sort_space[i] = permutation[i];
	//fprintf(stderr, "i %d, per %d, mper %d\n", i, permutation[i], permutation[i] << 3 >> 3);
    }

    for (int64_t iteration = words_per_substring - 1; iteration >= 0; --iteration) {
        uint32_t *lv2_substr_p = lv2_substrings + width * iteration;
#pragma omp parallel for
        for (uint32_t i = 0; i < lv2_num_items; ++i) {
            cpu_sort_space[i] &= 0xFFFFFFFF;
            cpu_sort_space[i] |= uint64_t(*(lv2_substr_p + (cpu_sort_space[i] & 0x1FFFFFFF))) << 32;
        }
        // pss::parallel_stable_sort(cpu_sort_space, cpu_sort_space + lv2_num_items, CompareHigh32Bits());
        __gnu_parallel::stable_sort(cpu_sort_space, cpu_sort_space + lv2_num_items, CompareHigh32Bits());
    }

    // copy answer back to host
#pragma omp parallel for
    for (uint32_t i = 0; i < lv2_num_items; ++i) {
        permutation[i] = cpu_sort_space[i] & 0xFFFFFFFF;
    }
}

