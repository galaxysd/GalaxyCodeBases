/*
 * 
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#ifndef __HEAP_RJ_H
#define __HEAP_RJ_H

#include "list.h"

define_list(rjheapv, void*);

typedef int (*heap_comp_func)(const void *e1, const void *e2, void *ref);

typedef struct {
	rjheapv *ptrs;
	void *ref;
	heap_comp_func cmp;
} Heap;

static inline Heap* init_heap(heap_comp_func cmp, void *ref){
	Heap *heap;
	heap = malloc(sizeof(Heap));
	heap->ptrs = init_rjheapv(8);
	heap->cmp  = cmp;
	heap->ref  = ref;
	return heap;
}

static inline void free_heap(Heap *heap){ free_rjheapv(heap->ptrs); free(heap); }

static inline void clear_heap(Heap *heap){ clear_rjheapv(heap->ptrs); }

static inline void push_heap(Heap *heap, void *p){
	void *pp;
	size_t i;
	i = count_rjheapv(heap->ptrs);
	push_rjheapv(heap->ptrs, p);
	while(i && heap->cmp(get_rjheapv(heap->ptrs, i), get_rjheapv(heap->ptrs, (i - 1) >> 1), heap->ref) < 0){
		pp = get_rjheapv(heap->ptrs, i);
		set_rjheapv(heap->ptrs, i, get_rjheapv(heap->ptrs, (i - 1) >> 1));
		set_rjheapv(heap->ptrs, (i - 1) >> 1, pp);
		i = (i - 1) >> 1;
	}
}

static inline size_t count_heap(Heap *heap){ return count_rjheapv(heap->ptrs); }

static inline void* peer_heap(Heap *heap){ return (count_rjheapv(heap->ptrs)? get_rjheapv(heap->ptrs, 0) : NULL );}

static inline void remove_heap(Heap *heap, size_t idx){
	void *pp;
	size_t swap;
	set_rjheapv(heap->ptrs, idx, get_rjheapv(heap->ptrs, count_rjheapv(heap->ptrs) - 1));
	trunc_rjheapv(heap->ptrs, 1);
	while((idx << 1) + 1 < count_rjheapv(heap->ptrs)){
		swap = idx;
		if(heap->cmp((const void*)get_rjheapv(heap->ptrs, swap), (const void*)get_rjheapv(heap->ptrs, (idx << 1) + 1), heap->ref) > 0){
			swap = (idx << 1) + 1;
		}
		if((idx << 1) + 2 < count_rjheapv(heap->ptrs) && heap->cmp((const void*)get_rjheapv(heap->ptrs, swap), (const void*)get_rjheapv(heap->ptrs, (idx << 1) + 2), heap->ref) > 0){
			swap = (idx << 1) + 2;
		}
		if(swap == idx) break;
		pp = get_rjheapv(heap->ptrs, idx);
		set_rjheapv(heap->ptrs,  idx, get_rjheapv(heap->ptrs, swap));
		set_rjheapv(heap->ptrs, swap, pp);
		idx = swap;
	}
}

static inline void* pop_heap(Heap *heap){
	void *p;
	if(count_rjheapv(heap->ptrs)){
		p = get_rjheapv(heap->ptrs, 0);
		remove_heap(heap, 0);
		return p;
	} else return NULL;
}

#endif
