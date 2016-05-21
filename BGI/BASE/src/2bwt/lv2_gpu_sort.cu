/*

   lv2_gpu_sort.cu       sort lev2 with GPU

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

#include <stdio.h>
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>

#include <omp.h>
#include <pthread.h>
#include <assert.h>
#include <stdarg.h>
#include <vector>
#include <algorithm>


#include "lv2_gpu_sort.h"

#include "b40c/radix_sort/enactor.cuh"
#include "b40c/util/multiple_buffering.cuh"

size_t get_gpu_mem()
{
    size_t free_gpu_mem, total_gpu_mem;
    assert(cudaMemGetInfo(&free_gpu_mem, &total_gpu_mem) == cudaSuccess);
    fprintf(stderr, "Free GPU mem: %lld\n", free_gpu_mem);
    return free_gpu_mem;
}

__global__ void permutation_kernel( uint32_t* index, uint32_t* val, uint32_t* new_val, int num_elements ) {
    int tid = blockIdx.x * THREADS_PER_BLOCK + threadIdx.x;

    if (tid < num_elements)
      new_val[tid] = val[index[tid] & 0x1FFFFFFF];
}


void lv2_gpu_sort(uint32_t *lv2_substrings, uint32_t *permutation, int uint32_ts_per_substring, int64_t width, int64_t lv2_num_items){

	uint32_t *d_keys, *d_values;
        assert(cudaMalloc( (void**) &d_keys, sizeof(uint32_t) * lv2_num_items ) == cudaSuccess);
        assert(cudaMalloc( (void**) &d_values, sizeof(uint32_t) * lv2_num_items ) == cudaSuccess);
        b40c::radix_sort::Enactor enactor;
        b40c::util::DoubleBuffer<uint32_t, uint32_t> ss(d_keys, d_values);

        assert(cudaMemcpy( ss.d_values[ss.selector], permutation, sizeof(uint32_t) * lv2_num_items, cudaMemcpyHostToDevice ) == cudaSuccess);
        for (int iteration = uint32_ts_per_substring - 1; iteration >= 0; --iteration ) { // TODO uint32_ts_per_suffix ?
        	if ( iteration == uint32_ts_per_substring - 1 ) { // TODO uint32_ts_per_suffix
                	assert(cudaMemcpy( ss.d_keys[ss.selector], lv2_substrings + ( iteration * width ),
                        sizeof(uint32_t) * lv2_num_items, cudaMemcpyHostToDevice ) == cudaSuccess);
                } else {
                	assert(cudaMemcpy( ss.d_keys[1-ss.selector], lv2_substrings + ( iteration * width ),
                        sizeof(uint32_t) * lv2_num_items, cudaMemcpyHostToDevice ) == cudaSuccess);
                        int num_blocks = ( lv2_num_items + THREADS_PER_BLOCK - 1 ) / THREADS_PER_BLOCK;

                        permutation_kernel<<<num_blocks, THREADS_PER_BLOCK>>>( ss.d_values[ss.selector], ss.d_keys[1-ss.selector],
                                                                               ss.d_keys[ss.selector], lv2_num_items );
                }
                assert(enactor.Sort<b40c::radix_sort::LARGE_SIZE>( ss, lv2_num_items ) == cudaSuccess);
        }

        // free device memory EXCEPT sort_indexes
        if (ss.d_keys[ss.selector]) cudaFree(ss.d_keys[ss.selector]);
        if (ss.d_keys[1-ss.selector]) cudaFree(ss.d_keys[1-ss.selector]);
        if (ss.d_values[1-ss.selector]) cudaFree(ss.d_values[1-ss.selector]);
        ///////////////////// END GPU SORT ////////////////////////////
        assert(cudaMemcpy( permutation, ss.d_values[ss.selector], sizeof(int) * lv2_num_items, cudaMemcpyDeviceToHost ) == cudaSuccess);
	if (ss.d_values[ss.selector]) cudaFree(ss.d_values[ss.selector]);
}

