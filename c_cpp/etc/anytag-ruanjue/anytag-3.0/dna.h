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
 
#ifndef __DNA_RJ_H
#define __DNA_RJ_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bitvec.h"

static const uint8_t base_bit_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static const uint8_t base_bit4_table[256] = {
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,

	15,  1, 14,  2,  13, 15, 15,  4,  11, 15, 15, 12,  15,  3, 15, 15,
	15, 15,  5,  6,   8, 15,  7,  9,  15, 10, 15, 15,  15, 15, 15, 15,
	15,  1, 14,  2,  13, 15, 15,  4,  11, 15, 15, 12,  15,  3, 15, 15,
	15, 15,  5,  6,   8, 15,  7,  9,  15, 10, 15, 15,  15, 15, 15, 15,

	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,

	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
	15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15
};

static const uint8_t bit4_bit_table[16] = { 4, 0, 1, 4,  2, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4 };

static const char bit_base_table[6] = "ACGTN-";
static const char bit4_base_table[16] = "-ACMGRSVTWYHKDBN";

static inline uint64_t dna_xor2ones(uint64_t seq){
	return ((seq & 0xAAAAAAAAAAAAAAAALLU) >> 1) | (seq & 0x5555555555555555LLU);
}

static inline uint64_t dna_rev_seq(uint64_t seq, uint8_t seq_size){
	seq = ~seq;
	seq = ((seq & 0x3333333333333333LLU)<< 2) | ((seq & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
	seq = ((seq & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
	seq = ((seq & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq & 0xFF00FF00FF00FF00LLU)>> 8);
	seq = ((seq & 0x0000FFFF0000FFFFLLU)<<16) | ((seq & 0xFFFF0000FFFF0000LLU)>>16);
	seq = ((seq & 0x00000000FFFFFFFFLLU)<<32) | ((seq & 0xFFFFFFFF00000000LLU)>>32);
	return seq >> (64 - (seq_size<<1));
}

static inline uint64_t seq2kmer(char *seq, uint32_t ksize){
	uint64_t kmer;
	uint32_t i;
	kmer = 0;
	for(i=0;i<ksize;i++) kmer = (kmer << 2) | base_bit_table[(int)seq[i]];
	return kmer;
}

static inline uint64_t seq2revkmer(char *seq, uint32_t ksize){
	uint64_t kmer;
	uint32_t i;
	kmer = 0;
	for(i=0;i<ksize;i++) kmer = (kmer << 2) | ((~base_bit_table[(int)seq[ksize - 1 - i]]) & 0x03);
	return kmer;
}

static inline void kmer2seq(char *seq, uint64_t kmer, uint32_t ksize){
	uint32_t i;
	for(i=0;i<ksize;i++){
		seq[i] = bit_base_table[(kmer >> ((ksize - 1 - i) << 1)) & 0x03];
	}
	seq[i] = 0;
}

static inline void kmer2revseq(char *seq, uint64_t kmer, uint32_t ksize){
	uint32_t i;
	kmer = ~kmer;
	for(i=0;i<ksize;i++){
		seq[i] = bit_base_table[(kmer >> (i << 1)) & 0x03];
	}
	seq[i] = 0;
}

#define kmer_mask(ksize) (0xFFFFFFFFFFFFFFFFLLU >> ((32 - (ksize)) * 2))

#define beg_seq2kmers(seq, seqlen, ksize, kmask, kmer, idx) {	\
kmer = 0;	\
for(idx=0;idx<ksize-1;idx++) kmer = (((kmer) << 2) | base_bit_table[(int)(seq)[idx]]);	\
for(idx=0;idx<=seqlen-ksize;idx++){	\
	kmer = ((kmer << 2) | base_bit_table[(int)(seq)[idx + ksize - 1]]) & kmask;
#define end_seq2kmers } }

#define beg_seq2revkmers(seq, seqlen, ksize, kmask, kmer, idx) {	\
kmer = 0;	\
for(idx=0;idx<ksize-1;idx++) kmer = (((kmer) << 2) | base_bit_table[(int)(seq)[seqlen - 1 - idx]]);	\
for(idx=0;idx<=seqlen-ksize;idx++){	\
	kmer = ((kmer << 2) | base_bit_table[(int)(seq)[seqlen - idx - ksize]]) & kmask;
#define end_seq2kmers } }

static inline void reverse_dna(char *seq, int len){
	int i, j;
	char c;
	i = 0;
	j = len - 1;
	while(i < j){
		c = seq[i]; seq[i] = seq[j]; seq[j] = c;
		i ++; j --;
	}
	for(i=0;i<len;i++){
		switch(seq[i]){
			case 'a': seq[i] = 't'; break;
			case 'A': seq[i] = 'T'; break;
			case 'c': seq[i] = 'g'; break;
			case 'C': seq[i] = 'G'; break;
			case 'g': seq[i] = 'c'; break;
			case 'G': seq[i] = 'C'; break;
			case 't': seq[i] = 'a'; break;
			case 'T': seq[i] = 'A'; break;
		}
	}
}

#define bit2bits(bits, off, bit) { if(((off) & 0x1FU) == 0) (bits)[(off) >> 5] = 0; (bits)[(off) >> 5] |= ((uint64_t)(bit)) << (((~(off)) & 0x1FU) << 1); }

static inline void seq2bits(uint64_t *bits, uint64_t bitoff, char *seq, uint32_t seqlen){
	uint64_t i, c;
	for(i=0;i<seqlen;i++){
		c = base_bit_table[(int)seq[i]];
		if(c == 4) c = lrand48() & 0x03;
		bit2bits(bits, bitoff + i, c);
	}
}

static inline void revseq2bits(uint64_t *bits, uint64_t bitoff, char *seq, uint32_t seqlen){
	uint64_t i, c;
	for(i=0;i<seqlen;i++){
		c = base_bit_table[(int)seq[seqlen - i - 1]];
		if(c == 4) c = lrand48();
		c = (~c) & 0x03;
		bit2bits(bits, bitoff + i, c);
	}
}

#define bits2bit(bits, off) (((bits)[(off) >> 5] >> (((~(off)) & 0x1FU) << 1)) & 0x03U)
#define bits2revbit(bits, off) ((~((bits)[(off) >> 5] >> (((~(off)) & 0x1FU) << 1))) & 0x03U)

static inline void bits2seq(char *seq, uint64_t *bits, uint64_t off, uint32_t len){
	uint32_t i, c;
	for(i=0;i<len;i++){
		c = bits2bit(bits, off + i);
		seq[i] = bit_base_table[c];
	}
	seq[i] = 0;
}

static inline void bits2revseq(char *seq, uint64_t *bits, uint64_t off, uint32_t len){
	uint32_t i, c;
	for(i=0;i<len;i++){
		c = (bits[(off + i)>>5] >> (((~(off + i)) & 0x1FU) << 1)) & 0x03;
		seq[len - i - 1] = bit_base_table[(~c)&0x03];
	}
	seq[i] = 0;
}

static inline uint64_t sub32seqbits(uint64_t *src, uint64_t off){
	if((off & 0x1F) == 0){
		return src[off>>5];
	} else {
		return (src[off>>5] << ((off & 0x1F) << 1)) | (src[(off>>5)+1] >> ((32 - (off & 0x1F)) << 1));
	}
}

static inline uint64_t sub4seqbits(uint64_t *src, uint64_t off){
	if(((off) & 0x1FU) > 28){
		return (((src[off>>5] << 32) | (src[(off>>5) + 1] >> 32)) >> ((28 - ((off - 16) & 0x1FU)) << 1)) & 0xFFU;
	} else {
		return (src[off>>5] >> ((28 - (off & 0x1FU)) << 1)) & 0xFFU;
	}
}

typedef struct {
	uint64_t *bits;
	uint64_t size;
	uint64_t cap;
} BaseBank;

static inline BaseBank* init_basebank(){
	BaseBank *bnk;
	bnk = malloc(sizeof(BaseBank));
	bnk->size = 0;
	bnk->cap  = 256;
	bnk->bits = malloc(8 * (bnk->cap / 32));
	memset(bnk->bits, 0, 8 * (bnk->cap / 32));
	return bnk;
}

static inline void free_basebank(BaseBank *bnk){
	free(bnk->bits);
	free(bnk);
}

static inline void encap_basebank(BaseBank *bnk, uint64_t inc){
	uint64_t old;
	if(bnk->size + inc < bnk->cap) return;
	old = bnk->cap;
	while(bnk->size + inc > bnk->cap){
		if(bnk->cap < 0x3FFFFFFLLU){
			bnk->cap <<= 1;
		} else {
			bnk->cap += 0x3FFFFFFLLU;
		}
	}
	bnk->bits = realloc(bnk->bits, bnk->cap);
	memset(bnk->bits + (old / 32), 0, (bnk->cap - old) / 4);
}

static inline void clear_basebank(BaseBank *bnk){
	memset(bnk->bits, 0, ((bnk->size + 31) / 32) * 8);
	bnk->size = 0;
}

static inline size_t dump_basebank(BaseBank *bnk, FILE *out){
	size_t n;
	fwrite(&bnk->size, sizeof(uint64_t), 1, out);
	n = ((bnk->size + 31) / 32);
	fwrite(bnk->bits, sizeof(uint64_t), n, out);
	return n * 8 + 8;
}

static inline BaseBank* load_basebank(FILE *inp){
	BaseBank *bnk;
	size_t n;
	bnk = init_basebank();
	if(fread(&bnk->size, sizeof(uint64_t), 1, inp) != 1){ free_basebank(bnk); return NULL; }
	encap_basebank(bnk, 0);
	n = (bnk->size + 31) / 32;
	if(fread(bnk->bits, sizeof(uint64_t), n, inp) != n){ free_basebank(bnk); return NULL; }
	return bnk;
}

static inline void bits2basebank(BaseBank *bnk, uint64_t *bits, uint64_t off, uint64_t len){
	uint64_t offset;
	encap_basebank(bnk, len);
	for(offset=off;offset<off+len;offset++){
		bit2bits(bnk->bits, bnk->size, bits2bit(bits, offset));
		bnk->size ++;
	}
}

static inline void revbits2basebank(BaseBank *bnk, uint64_t *bits, uint64_t off, uint64_t len){
	uint64_t offset;
	encap_basebank(bnk, len);
	for(offset=off+len;offset>off;offset--){
		bit2bits(bnk->bits, bnk->size, (~bits2bit(bits, offset - 1)) & 0x03);
		bnk->size ++;
	}
}

static inline void seq2basebank(BaseBank *bnk, char *seq, uint64_t len){
	char *p;
	uint8_t c;
	p = seq;
	seq = seq + len;
	encap_basebank(bnk, len);
	while(p < seq){
		c = base_bit_table[(int)*p];
		if(c == 4) c = lrand48() & 0x03;
		bit2bits(bnk->bits, bnk->size, c);
		bnk->size ++;
		p ++;
	}
}

static inline void revseq2basebank(BaseBank *bnk, char *seq, uint64_t len){
	char *p;
	uint8_t c;
	p = seq + len;
	encap_basebank(bnk, len);
	while(p > seq){
		c = base_bit_table[(int)*p];
		if(c == 4) c = lrand48() & 0x03;
		c = (~c) & 0x03;
		bit2bits(bnk->bits, bnk->size, c);
		p --;
		bnk->size ++;
	}
}

static inline void seq_basebank(BaseBank *bnk, uint64_t off, uint64_t len, char *seq){
	uint64_t i;
	for(i=0;i<len;i++){
		seq[i] = bit_base_table[bits2bit(bnk->bits, off + i)];
	}
	seq[i] = 0;
}

static inline void revseq_basebank(BaseBank *bnk, uint64_t off, uint64_t len, char *seq){
	uint64_t i;
	for(i=0;i<len;i++){
		seq[i] = bit_base_table[(~bits2bit(bnk->bits, off + len - 1 - i)) & 0x03];
	}
	seq[i] = 0;
}

static inline uint64_t sub32_basebank(BaseBank *bnk, uint64_t off){ return sub32seqbits(bnk->bits, off); }

static inline uint64_t sub4_basebank(BaseBank *bnk, uint64_t off){ return sub4seqbits(bnk->bits, off); }

static inline uint32_t mismatch_basebank(BaseBank *bnk, uint64_t off1, uint64_t off2, uint32_t len){
	uint64_t seq1, seq2;
	uint32_t mm, i;
	mm = 0;
	for(i=0;i+32<=len;i+=32){
		seq1 = sub32seqbits(bnk->bits, off1 + i);
		seq2 = sub32seqbits(bnk->bits, off2 + i);
		mm += count_ones_bit64(dna_xor2ones(seq1 ^ seq2));
	}
	if(i < len){
		seq1 = sub32seqbits(bnk->bits, off1 + i);
		seq2 = sub32seqbits(bnk->bits, off2 + i);
		mm += count_ones_bit64((dna_xor2ones(seq1 ^ seq2)) >> ((32 - (len - i)) << 1));
	}
	return mm;
}

#endif
