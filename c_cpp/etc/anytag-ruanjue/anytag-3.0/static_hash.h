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
 
#ifndef __STATIC_HASH_H_RJ
#define __STATIC_HASH_H_RJ

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "hashset.h"
#include "bitvec.h"

#define define_static_hashset(hash_type, primitive_hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
typedef struct {	\
	BitVec *flags;	\
	hash_key_type *array;	\
	size_t size, count, iter;	\
} hash_type;	\
static inline hash_type* primitive_hash_type##2##hash_type(primitive_hash_type *old){	\
	hash_type *set;	\
	hash_key_type *e;	\
	size_t i;	\
	assert(old->ocp == old->count);	\
	set = malloc(sizeof(hash_type));	\
	set->size = old->size;	\
	set->count = old->count;	\
	set->flags = init_bitvec(set->size);	\
	set->array = malloc(sizeof(hash_key_type) * set->count);	\
	if(set->array == NULL){	\
		fprintf(stderr, " -- Out of memory --\n"); fflush(stderr); abort();	\
	}	\
	reset_iter_##primitive_hash_type(old);	\
	i = 0;	\
	while((e = ref_iter_##primitive_hash_type(old))){	\
		set->array[i++] = *e;	\
		one_bitvec(set->flags, offset_##primitive_hash_type(old, e));	\
	}	\
	assert(i == set->count);	\
	index_bitvec(set->flags);	\
	return set;	\
}	\
static inline int64_t count_##hash_type(hash_type *set){ return set->count; }	\
static inline hash_key_type* get_##hash_type(hash_type *set, hash_key_type key){	\
	size_t hc, ec;	\
	hc = hash_code_macro(key) % set->size;	\
	ec = set->size;	\
	while(1){	\
		if(get_bitvec(set->flags, hc)){	\
			if(ec == set->size) ec = rank_bitvec(set->flags, hc);	\
			if(hash_equal_macro(set->array[ec], key)) return set->array + ec;	\
		} else {	\
			return NULL;	\
		}	\
		if(hc + 1 == set->size){	\
			hc = ec = 0;	\
		} else { hc ++; ec ++; }	\
	}	\
}	\
static inline void reset_iter_##hash_type(hash_type *set){ set->iter = 0; } \
static inline hash_key_type* ref_iter_##hash_type(hash_type *set){	\
	if(set->iter < set->count) return set->array + set->iter ++;	\
	else return NULL;	\
}	\
static inline size_t offset_##hash_type(hash_type *set, hash_key_type *ptr){ return ptr - set->array; }	\
static inline size_t dump_##hash_type(hash_type *set, FILE *out){              \
	size_t n;                                                          \
	n = 0;	\
	n += ffwrite(&set->size, sizeof(size_t), 1, out);                         \
	n += ffwrite(&set->count, sizeof(size_t), 1, out);                         \
	n += ffwrite(set->array, sizeof(hash_key_type), set->count, out);	\
	n += dump_bitvec(set->flags, out);	\
	return n;                                                                  \
}	\
static inline hash_type* load_##hash_type(FILE *in){	\
	hash_type *set;	\
	set = malloc(sizeof(hash_type));	\
	set->iter = 0;	\
	fread(&set->size, sizeof(size_t), 1, in);	\
	fread(&set->count, sizeof(size_t), 1, in);	\
	set->array = malloc(sizeof(hash_key_type) * set->count);	\
	fread(set->array, sizeof(hash_key_type), set->count, in);	\
	set->flags = load_bitvec(in);	\
	index_bitvec(set->flags);	\
	return set;	\
}	\
static inline void free_##hash_type(hash_type *set){	\
	free_bitvec(set->flags);	\
	free(set->array);	\
	free(set);	\
}

// ---- basic static hash --- //
define_static_hashset(U32hash, u32hash, uint32_t, u32hash_code, uxxhash_equals);
define_static_hashset(U64hash, u64hash, uint64_t, u64hash_code, uxxhash_equals);
define_static_hashset(I32hash, i32hash, int, i32hash_code, i32hash_equals);
define_static_hashset(Chash, chash, char*, chash_code, chash_equals);
define_static_hashset(UUhash, uuhash, uuhash_t, uuhash_code, uuhash_equals);
define_static_hashset(CUhash, cuhash, cuhash_t, cuhash_code, cuhash_equals);

#endif
