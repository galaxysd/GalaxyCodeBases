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
 
#ifndef __BITS_VEC_RJ_H
#define __BITS_VEC_RJ_H

#include <stdint.h>
#include <string.h>
#include <stdlib.h>

/* Useful functions when n_bit > 8 */

static inline void u8byte2bits(uint64_t val, uint8_t *dat, uint64_t offset, uint8_t size){
	uint8_t i, *v;
	v = (uint8_t*)&val;
#if __BYTE_ORDER == 1234
	for(i=0;i<size;i++){
		if(v[i >> 3] & (1U << (i & 0x7U))) dat[offset >> 3] |= 1U << (offset & 0x7U);
		else dat[offset >> 3] &= ~(1U << (offset & 0x7U));
		offset ++;
	}
#else
	for(i=0;i<size;i++){
		if(v[7 - (i >> 3)] & (1U << (i & 0x7U))) dat[offset >> 3] |= 1U << (offset & 0x7U);
		else dat[offset >> 3] &= ~(1U << (offset & 0x7U));
		offset ++;
	}
#endif
}

static inline uint64_t bits2u8byte(uint8_t *dat, uint64_t offset, uint8_t size){
	uint64_t ret;
	uint8_t i, *v;
	ret = 0;
	v = (uint8_t*)&ret;
#if __BYTE_ORDER == 1234
	for(i=0;i<size;i++){
		if(dat[offset >> 3] & (1U << (offset & 0x7U))) v[i >> 3] |= 1U << (i & 0x7U);
		offset ++;
	}
#else
	for(i=0;i<size;i++){
		if(dat[offset >> 3] & (1U << (offset & 0x7U))) v[7 - (i >> 3)] |= 1U << (i & 0x7U);
		offset ++;
	}
#endif
	return ret;
}

typedef struct {
	uint8_t *bits;
	uint64_t size;
	uint64_t cap;
	uint32_t n_bit, mask;
} BitsVec;

static inline BitsVec* init_bitsvec(uint64_t size, uint32_t n_bit){
	BitsVec *vec;
	if(n_bit == 0) n_bit = 1;
	else if(n_bit > 8) n_bit = 8;
	if(size < 8) size = 8;
	vec = malloc(sizeof(BitsVec));
	vec->n_bit = n_bit;
	vec->mask  = (1U << n_bit) - 1U;
	vec->size  = 0;
	vec->cap   = size;
	vec->bits  = calloc(1, (size * vec->n_bit + 7) / 8);
	return vec;
}

static inline void clear_bitsvec(BitsVec *vec){
	memset(vec->bits, 0, (vec->cap * vec->n_bit + 7) / 8);
	vec->size = 0;
}

static inline void free_bitsvec(BitsVec *vec){
	free(vec->bits);
	free(vec);
}

static inline int encap_bitsvec(BitsVec *vec, uint32_t n){
	uint64_t cap;
	if(vec->size + n <= vec->cap) return 0;
	cap = vec->cap;
	while(vec->size + n > vec->cap){
		if(vec->cap < 1024 * 1024){
			vec->cap <<= 1;
		} else {
			vec->cap += 1024 * 1024;
		}
	}
	vec->bits = realloc(vec->bits, (vec->cap * vec->n_bit + 7) / 8);
	memset(vec->bits + (cap * vec->n_bit + 7) / 8, 0, (vec->cap * vec->n_bit + 7) / 8 - (cap * vec->n_bit + 7) / 8);
	return 1;
}

static inline void set_bitsvec(BitsVec *vec, uint64_t idx, uint8_t dat){
	uint64_t off;
	off = (idx * vec->n_bit);
	dat &= vec->mask;
	vec->bits[off >> 3] = (vec->bits[off >> 3] & (~(vec->mask << (off & 0x07U)))) | (dat << (off & 0x07U));
	if((off & 0x07U) + vec->n_bit > 8){
		vec->bits[(off >> 3) + 1] = (vec->bits[(off >> 3) + 1] & (~(vec->mask >> (8 - (off & 0x07U))))) | (dat >> (8 - (off & 0x07U)));
		//vec->bits[(off >> 3) + 1] |= dat >> (8 - (off & 0x07U));
	}
}

static inline void push_bitsvec(BitsVec *vec, uint8_t dat){
	encap_bitsvec(vec, 1);
	set_bitsvec(vec, vec->size, dat);
	vec->size ++;
}

static inline uint8_t get_bitsvec(BitsVec *vec, uint64_t idx){
	uint64_t off;
	uint8_t dat;
	off = (idx * vec->n_bit);
	dat = (vec->bits[off >> 3] >> (off & 0x07U));
	if((off & 0x07U) + vec->n_bit > 8){
		dat |= vec->bits[(off >> 3) + 1] << (8 - (off & 0x07U));
	}
	return dat & vec->mask;
}

static inline int pop_bitvec(BitsVec *vec){
	if(vec->size == 0) return -1;
	vec->size --;
	return get_bitsvec(vec, vec->size);
}

#endif
