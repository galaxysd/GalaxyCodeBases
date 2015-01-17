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
 
#ifndef __ANYTAG_TAG_RJ_H
#define __ANYTAG_TAG_RJ_H

#include "hashset.h"
#include "file_reader.h"
#include "list.h"
#include "dna.h"

#define MAX_KMER_SIZE	59
#define MAX_RD_LEN	1024
#define MAX_TAG_CNT ((1<<10)-1)

typedef struct {
	uint64_t k1, k2:54, cnt:10;
} Kmer;

#define kmer_hashcode(K) u64hash_code(((K).k1 + (K).k2))
#define kmer_equals(K1, K2) ((K1).k1 == (K2).k1 && (K1).k2 == (K2).k2)
define_hashset(kmerhash, Kmer, kmer_hashcode, kmer_equals);

#define kmer2code(K) (((K).k1 ^ ((K).k1 >> 1)) ^ ((K).k1 >> 11))

#ifdef __CPLUSPLUS
extern "C" {
#endif

static inline int min_value_kmer_dir(const uint8_t *seq, uint32_t len){
	uint32_t i;
	uint8_t rev;
	for(i=0;i<len;i++){
		rev = (~seq[len-i-1]) & 0x03;
		if(seq[i] < rev) return 0;
		if(seq[i] > rev) return 1;
	}
	return 0;
}

#ifdef __CPLUSPLUS
}
#endif

#endif
