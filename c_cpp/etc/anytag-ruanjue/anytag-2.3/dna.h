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

static inline void bit2bits(uint64_t *bits, uint64_t bitoff, uint8_t c){
	if(c == 4) c = lrand48() & 0x03;
	if((bitoff & 0x1FU) == 0) bits[bitoff >> 5] = 0;
	bits[bitoff >> 5] |= c << (((~bitoff) & 0x1FU) << 1);
}

static inline void seq2bits(uint64_t *bits, uint64_t bitoff, char *seq, uint32_t seqlen){
	uint64_t i, idx, off, c;
	for(i=0;i<seqlen;i++){
		idx = (i + bitoff) >> 5;
		off = (i + bitoff) & 0x1FU;
		c = base_bit_table[(int)seq[i]];
		if(c == 4) c = lrand48() & 0x03;
		if(off == 0) bits[idx] = 0;
		bits[idx] |= c << (((~off)&0x1FU) << 1);
	}
}

static inline void revseq2bits(uint64_t *bits, uint64_t bitoff, char *seq, uint32_t seqlen){
	uint64_t i, idx, off, c;
	for(i=0;i<seqlen;i++){
		idx = (i + bitoff) >> 5;
		off = (i + bitoff) & 0x1FU;
		c = base_bit_table[(int)seq[seqlen - i - 1]];
		if(c == 4) c = lrand48();
		c = (~c) & 0x03;
		if(off == 0) bits[idx] = 0;
		bits[idx] |= c << (((~off)&0x1FU) << 1);
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

#endif
