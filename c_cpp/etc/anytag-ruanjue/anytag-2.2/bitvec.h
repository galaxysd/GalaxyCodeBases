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
 
#ifndef __BIT_VEC_RJ_H
#define __BIT_VEC_RJ_H

#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define get_bit8(bits, idx) ((((bits)[(idx) >> 3]) >> ((idx) & 0x07)) & 0x01)
#define get_bit16(bits, idx) ((((bits)[(idx) >> 4]) >> ((idx) & 0x0F)) & 0x01)
#define get_bit32(bits, idx) ((((bits)[(idx) >> 5]) >> ((idx) & 0x1F)) & 0x01)
#define get_bit64(bits, idx) ((((bits)[(idx) >> 6]) >> ((idx) & 0x3F)) & 0x01)

#define get_2bit8(bits, idx) ((((bits)[(idx) >> 2]) >> (((idx) & 0x03) << 1)) & 0x03)
#define get_2bit16(bits, idx) ((((bits)[(idx) >> 3]) >> (((idx) & 0x07) << 1)) & 0x03)
#define get_2bit32(bits, idx) ((((bits)[(idx) >> 4]) >> (((idx) & 0x0F) << 1)) & 0x03)
#define get_2bit64(bits, idx) ((((bits)[(idx) >> 5]) >> (((idx) & 0x1F) << 1)) & 0x03)

static const uint8_t byte_ones_table[256] = {
	0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

typedef struct {
	uint64_t *bits;
	uint64_t n_bit;
	uint64_t n_cap;
	uint64_t *sums;
	uint64_t sum_size;
	uint64_t n_ones;
	uint64_t *hash;
	uint64_t hash_size;
	uint64_t hash_mod;
	int64_t iter_idx;
} BitVec;

static inline uint32_t count_ones_bit32(uint32_t v){
	v = v - ((v >> 1) & 0x55555555U);                        // reuse input as temporary
	v = (v & 0x33333333U) + ((v >> 2) & 0x33333333U);        // temp
	return (((v + (v >> 4)) & 0xF0F0F0FU) * 0x1010101U) >> 24; // count
}

#define ONES_STEP_4 0x1111111111111111ULL
#define ONES_STEP_8 0x0101010101010101ULL

static inline int count_ones_bit64(const uint64_t x){
	register uint64_t byte_sums = x - ((x & 0xa * ONES_STEP_4) >> 1);
	byte_sums = (byte_sums & 3 * ONES_STEP_4) + ((byte_sums >> 2) & 3 * ONES_STEP_4);
	byte_sums = (byte_sums + (byte_sums >> 4)) & 0x0f * ONES_STEP_8;
	return byte_sums * ONES_STEP_8 >> 56;
}

static inline BitVec* init_bitvec(uint64_t n_bit){
	BitVec *bitv;
	if(n_bit == 0) n_bit = 64 * 8;
	bitv = (BitVec*)malloc(sizeof(BitVec));
	bitv->n_bit = 0;
	bitv->n_cap = (((n_bit + 63) / 64) + 7) / 8 * 64 * 8;
	bitv->bits  = (uint64_t*)malloc(bitv->n_cap / 8);
	memset(bitv->bits, 0, bitv->n_cap / 8);
	bitv->sums = NULL;
	bitv->hash = NULL;
	bitv->hash_size = 0;
	return bitv;
}

static inline void clear_bitvec(BitVec *bitv){ bitv->n_bit = 0; }

static inline void zeros_bitvec(BitVec *bitv){ memset(bitv->bits, 0, bitv->n_cap / 8); }

static inline void ones_bitvec(BitVec *bitv){ memset(bitv->bits, 0xFFU, bitv->n_cap / 8); }

static inline void flip_bitvec(BitVec *bitv, uint64_t idx){ bitv->bits[idx>>6] ^= 1LLU << (idx&0x3FU); }

static inline void one_bitvec(BitVec *bitv, uint64_t idx){ bitv->bits[idx>>6] |= 1LLU << (idx&0x3FU); }

static inline void zero_bitvec(BitVec *bitv, uint64_t idx){ bitv->bits[idx>>6] &= ~(1LLU << (idx&0x3FU)); }

static inline uint64_t get_bitvec(BitVec *bitv, uint64_t idx){ return (bitv->bits[idx>>6] >> (idx&0x3FU)) & 0x01LLU; }

static inline void encap_bitvec(BitVec *bitv, uint64_t num){
	uint64_t cap;
	if(bitv->n_bit + num < bitv->n_cap) return;
	cap = bitv->n_cap;
	while(bitv->n_bit + num >= bitv->n_cap){
		if(bitv->n_cap < 1024 * 1024 * 8){
			bitv->n_cap <<= 1;
		} else bitv->n_cap += 1024 * 1024 * 8;
	}
	bitv->bits = (uint64_t*)realloc(bitv->bits, bitv->n_cap / 8 + 8);
	memset(((void*)bitv->bits) + cap / 8, 0, (bitv->n_cap - cap) / 8 + 8);
	bitv->bits[cap / 64] = 0x0000000000000001LLU;
}

static inline void one2bitvec(BitVec *bitv){ encap_bitvec(bitv, 1); one_bitvec(bitv, bitv->n_bit); bitv->n_bit ++; }

static inline void zero2bitvec(BitVec *bitv){ encap_bitvec(bitv, 1); zero_bitvec(bitv, bitv->n_bit); bitv->n_bit ++; }

static inline uint64_t next_one_bitvec(BitVec *bitv, uint64_t idx){
	while((bitv->bits[idx >> 6] >> (idx & 0x3F)) == 0){
		idx = ((idx >> 6) + 1) << 6;
		if(idx >= bitv->n_cap) return 0xFFFFFFFFFFFFFFFFLLU;
	}
	while(((bitv->bits[idx >> 6] >> (idx & 0x3F)) & 0xFF) == 0) idx += 8;
	while(((bitv->bits[idx >> 6] >> (idx & 0x3F)) & 0x01) == 0) idx ++;
	return idx;
}

static inline void index_bitvec(BitVec *bitv){
	uint64_t i, k, s, t, m;
	if(bitv->sums) free(bitv->sums);
	m = ((bitv->n_cap + 63) / 64 + 7) / 8;
	bitv->sums = (uint64_t*)malloc((m * 2 + 1) * 8);
	memset(bitv->sums, 0, (m * 2 + 1) * 8);
	t = 0;
	for(i=0;i<bitv->n_cap;i+=64*8){
		k = ((i>>6) >> 3) << 1;
		bitv->sums[k] = t;
		s = 0;
		s += count_ones_bit64(bitv->bits[(i>>6)+0]);
		bitv->sums[k+1] |= s << 0;
		s += count_ones_bit64(bitv->bits[(i>>6)+1]);
		bitv->sums[k+1] |= s << 9;
		s += count_ones_bit64(bitv->bits[(i>>6)+2]);
		bitv->sums[k+1] |= s << 18;
		s += count_ones_bit64(bitv->bits[(i>>6)+3]);
		bitv->sums[k+1] |= s << 27;
		s += count_ones_bit64(bitv->bits[(i>>6)+4]);
		bitv->sums[k+1] |= s << 36;
		s += count_ones_bit64(bitv->bits[(i>>6)+5]);
		bitv->sums[k+1] |= s << 45;
		s += count_ones_bit64(bitv->bits[(i>>6)+6]);
		bitv->sums[k+1] |= s << 54;
		s += count_ones_bit64(bitv->bits[(i>>6)+7]);
		t += s;
	}
	bitv->sums[((i>>6) >> 3) << 1] = t;
	bitv->n_ones = t;
	bitv->sum_size = m;
	bitv->hash_size = (bitv->n_cap / 64 / 8) / 2;
	if(bitv->hash_size == 0) bitv->hash_size = 1;
	bitv->hash_mod = (t + bitv->hash_size - 1) / bitv->hash_size;
	if(bitv->hash_mod == 0) bitv->hash_mod = 1;
	bitv->hash = (uint64_t*)malloc(sizeof(uint64_t) * bitv->hash_size);
	s = 0;
	t = 0;
	for(i=0;i<=m;i++){
		k = bitv->sums[i*2] / bitv->hash_mod;
		if(s < k){
			while(s < k){ bitv->hash[s] = t; s ++; }
			t = i? i - 1 : 0;
		}
	}
	bitv->hash[bitv->sums[m*2] / bitv->hash_mod] = t;
}

static inline uint64_t rank_bitvec(BitVec *bitv, uint64_t idx){
	uint64_t p, s, sum;
	p = (idx>>6)>>3;
	s = (idx >> 6) & 0x07U;
	sum = bitv->sums[p<<1];
	if(s) sum += (bitv->sums[(p<<1)+1] >> (9 * (s - 1))) & 0x1FFU;
	if(idx & 0x3FU) sum += count_ones_bit64(bitv->bits[idx>>6]<<(64-(idx&0x3FU)));
	return sum;
}

static inline uint8_t select_8bytes(uint64_t word, uint8_t n_one){
	uint8_t idx, n, m;
	n = count_ones_bit32((uint32_t)word);
	if(n >= n_one){
		n = 0;
		idx = 0;
		word = word & 0xFFFFFFFFU;
	} else {
		idx = 32;
		word = word >> 32;
	}
	while(1){
		m = byte_ones_table[(uint8_t)word];
		if(n + m >= n_one) break;
		n += m;
		idx += 8;
		word >>= 8;
	}
	m = byte_ones_table[(uint8_t)(word & 0xF)];
	if(n + m < n_one){
		idx += 4;
		word >>= 4;
		n += m;
	}
	while(word){
		idx ++;
		if(word & 0x01){
			n ++;
			if(n == n_one) break;
		}
		word >>= 1;
	}
	return idx;
}

static inline uint64_t select_bitvec(BitVec *bitv, uint64_t idx){
	uint64_t i, p, s, sum, t;
	p = bitv->hash[idx / bitv->hash_mod];
	while(p + 1 < bitv->sum_size && bitv->sums[(p + 1) << 1] < idx) p ++;
	sum = bitv->sums[p << 1];
	i = 0;
	t = sum;
	while(i < 7){
		s = (bitv->sums[(p << 1) + 1] >> (9 * i)) & 0x1FFU;
		if(sum + s >= idx) break;
		t = sum + s;
		i ++;
	}
	p = p * 8 + i;
	s = idx - t;
	return p * 64 + select_8bytes(bitv->bits[p], s);
}

static inline void begin_iter_bitvec(BitVec *bitv){ bitv->iter_idx = -1; }

static inline uint64_t iter_bitvec(BitVec *bitv){
	if((uint64_t)(bitv->iter_idx + 1) > bitv->n_cap) return 0xFFFFFFFFFFFFFFFFLLU;
	bitv->iter_idx = next_one_bitvec(bitv, bitv->iter_idx + 1);
	return (uint64_t)bitv->iter_idx;
}

static inline void free_bitvec(BitVec *bitv){
	free(bitv->bits);
	if(bitv->sums) free(bitv->sums);
	if(bitv->hash) free(bitv->hash);
	free(bitv);
}

#endif
