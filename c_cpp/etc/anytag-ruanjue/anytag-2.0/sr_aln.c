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
 
#define _GNU_SOURCE
#include "sr_aln.h"
#include "stdaln.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

uint16_t sr_rdseq_length(SR_SeqDB *sdb, uint32_t rid){ return sdb->rd_len? sdb->rd_len : get_u16list(sdb->rd_lens, rid); }

uint64_t sr_rdseq_offset(SR_SeqDB *sdb, uint32_t rid){
	uint64_t off;
	uint16_t *p;
	if(sdb->rd_len) return ((uint64_t)sdb->rd_len) * rid + SR_RDSEQ_PADDING;
	off = get_u64list(sdb->rd_offs, rid / 16);
	p = ref_u16list(sdb->rd_lens, (rid / 16) * 16);
	rid = rid % 16;
	while(rid){ off += *p; rid --; p ++; }
	return off;
}

SR_AlnAux* sr_init_aux(){
	SR_AlnAux *aux;
	aux = malloc(sizeof(SR_AlnAux));
	aux->hits   = init_sr_hitv(32);
	aux->ptrs   = init_sr_ptrv(32);
	aux->pids   = init_u64list(1024);
	aux->fkmers = init_u32list(32);
	aux->rkmers = init_u32list(32);
	return aux;
}

void sr_free_aux(SR_AlnAux *aux){
	free_sr_hitv(aux->hits);
	free_sr_ptrv(aux->ptrs);
	free_u64list(aux->pids);
	free_u32list(aux->fkmers);
	free_u32list(aux->rkmers);
	free(aux);
}

SR_SeqDB* sr_init_sdb(uint8_t kmer_size, uint8_t seed_kmer_num, uint16_t rd_len){
	SR_SeqDB *sdb;
	uint32_t i, j;
	if(kmer_size > MAX_KMER_SIZE) kmer_size = MAX_KMER_SIZE;
	if(kmer_size < MIN_KMER_SIZE) kmer_size = MIN_KMER_SIZE;
	sdb = malloc(sizeof(SR_SeqDB));
	sdb->rd_seqs = init_u64list(SR_RDSEQ_PADDING / 32 + 32);
	for(i=0;i<SR_RDSEQ_PADDING/32;i++) set_u64list(sdb->rd_seqs, i, 0);
	sdb->offset = SR_RDSEQ_PADDING;
	sdb->rd_offs = init_u64list(16);
	sdb->rd_lens = init_u16list(16);
	sdb->rd_filtered = init_bitvec(128);
	sdb->rd_len  = rd_len;
	sdb->max_rd_len = rd_len;
	sdb->n_rd    = 0;
	sdb->kmer_size = kmer_size;
	sdb->kmer_mask = (1LLU << (kmer_size << 1)) - 1U;
	sdb->strand_type    = 2;
	sdb->min_overlap    = kmer_size * 2;
	sdb->min_similarity = 0.9;
	sdb->max_mismatch   = 10;
	sdb->max_hits_per_query = 1024;
	sdb->min_complexity = 12;
	sdb->allow_gap      = 0;
	sdb->filter_dup     = 1;
	sdb->ap             = (AlnParam){10, 3, 2, aln_sm_nt, 16, 75};

	sdb->n_seed = (seed_kmer_num < 2) ? 2 : seed_kmer_num;
	if(sdb->n_seed > MAX_N_SEED) sdb->n_seed = MAX_N_IDX;
	sdb->n_idx = 0;
	for(i=0;i<sdb->n_seed-1;i++){
		for(j=i+1;j<sdb->n_seed;j++){
			sdb->idx_pats[sdb->n_idx][0] = i;
			sdb->idx_pats[sdb->n_idx][1] = j;
			sdb->idx_pats[sdb->n_idx][2] = 0;
			sdb->n_idx ++;
			sdb->idx_pats[sdb->n_idx][0] = i;
			sdb->idx_pats[sdb->n_idx][1] = j;
			sdb->idx_pats[sdb->n_idx][2] = 1;
			sdb->n_idx ++;
		}
	}
	sdb->indexs      = malloc(sizeof(sr_prmhash*) * sdb->n_idx);
	sdb->fast_indexs = malloc(sizeof(u32list*) * sdb->n_idx);
	sdb->links       = malloc(sizeof(u32list*) * sdb->n_idx);
	sdb->is_fast     = (sdb->kmer_size <= 7);
	for(i=0;i<sdb->n_idx;i++){
		sdb->indexs[i] = init_sr_prmhash(1023);
		sdb->fast_indexs[i] = sdb->is_fast? init_u64list(1LLU << ((sdb->kmer_size << 1) << 1)) : NULL;
		sdb->links[i] = init_u32list(1024);
	}

	return sdb;
}

void sr_set_align_parameters(SR_SeqDB *sdb, uint8_t strand, uint8_t min_overlap, float min_sm, uint8_t max_mm, int allow_gap){
	sdb->strand_type    = strand;
	sdb->min_overlap    = min_overlap;
	sdb->min_similarity = min_sm;
	sdb->max_mismatch   = max_mm;
	sdb->allow_gap      = allow_gap;
}

void sr_set_filter_parameters(SR_SeqDB *sdb, uint32_t min_complexity, uint32_t max_hits_per_query, int filter_dup){
	sdb->min_complexity = min_complexity;
	sdb->max_hits_per_query = max_hits_per_query;
	sdb->filter_dup = filter_dup;
}

uint32_t cal_seq_complexity(char *seq, uint16_t seqlen){
	uint32_t i, tri;
	uint64_t tris;
	tris = 0;
	tri = 0;
	tri = (tri << 2) | base_bit_table[(int)seq[0]];
	tri = (tri << 2) | base_bit_table[(int)seq[1]];
	for(i=2;i<seqlen;i++){
		tri = (tri << 2) | base_bit_table[(int)seq[i]];
		tris |= 1LLU << (tri & 0x3F);
	}
	return count_ones_bit64(tris);
}

int sr_push_sdb(SR_SeqDB *sdb, char *seq, uint16_t seqlen){
	if(sdb->rd_len && seqlen != sdb->rd_len) return 0;
	encap_u64list(sdb->rd_seqs, (sdb->offset + seqlen) / 32 + 2);
	seq2bits(as_array_u64list(sdb->rd_seqs), sdb->offset, seq, seqlen);
	if(sdb->rd_len == 0){
		if(seqlen > sdb->max_rd_len) sdb->max_rd_len = seqlen;
		push_u16list(sdb->rd_lens, seqlen);
		if(sdb->n_rd % 16 == 0) push_u64list(sdb->rd_offs, sdb->offset);
	}
	sdb->offset += seqlen;
	if(cal_seq_complexity(seq,seqlen) < sdb->min_complexity) one2bitvec(sdb->rd_filtered);
	else zero2bitvec(sdb->rd_filtered);
	sdb->n_rd ++;
	return 1;
}

void sr_ready_sdb(SR_SeqDB *sdb){
	encap_u64list(sdb->rd_seqs, (sdb->offset + 128) / 32 + 2);
	return;
	sdb = sdb;
}

void sr_reset_sdb(SR_SeqDB *sdb){
	uint32_t i;
	for(i=0;i<sdb->n_idx;i++){
		clear_sr_prmhash(sdb->indexs[i]);
		clear_u32list(sdb->links[i]);
	}
	sdb->n_rd = 0;
	sdb->offset = SR_RDSEQ_PADDING;
	clear_u64list(sdb->rd_offs);
	clear_u16list(sdb->rd_lens);
	clear_bitvec(sdb->rd_filtered);
}

void sr_free_sdb(SR_SeqDB *sdb){
	uint32_t i;
	free_u64list(sdb->rd_seqs);
	free_u64list(sdb->rd_offs);
	free_u16list(sdb->rd_lens);
	for(i=0;i<sdb->n_idx;i++){
		free_sr_prmhash(sdb->indexs[i]);
		if(sdb->is_fast) free_u64list(sdb->fast_indexs[i]);
		free_u32list(sdb->links[i]);
	}
	free(sdb->indexs);
	free(sdb->fast_indexs);
	free(sdb->links);
	free_bitvec(sdb->rd_filtered);
	free(sdb);
}

void sr_index_sdb(SR_SeqDB *sdb, uint32_t n_idx){
	SR_Primer P, *p;
	uint32_t i, j, rd_len, c;
	uint64_t offset, k;
	int exists;
	P.rid = 0;
	if(sdb->is_fast){
		memset(ref_u64list(sdb->fast_indexs[n_idx], 0), 0xFFU, (1LLU << ((sdb->kmer_size << 1) << 1)) * sizeof(uint64_t));
	}
	for(i=0;i<sdb->n_rd;i++){
		P.kmer = 0;
		if(get_bitvec(sdb->rd_filtered, i)){
			push_u32list(sdb->links[n_idx], i);
			continue;
		}
		offset = sr_rdseq_offset(sdb, i);
		rd_len = sr_rdseq_length(sdb, i);
		if(rd_len < (sdb->idx_pats[n_idx][1] + 1) * sdb->kmer_size){
			push_u32list(sdb->links[n_idx], i);
			continue;
		}
		if(sdb->idx_pats[n_idx][2]){
			for(j=0;j<sdb->kmer_size;j++){
				c = (~bits2bit(as_array_u64list(sdb->rd_seqs), offset + rd_len - 1 - (sdb->idx_pats[n_idx][0] * sdb->kmer_size + j))) & 0x03;
				P.kmer = (P.kmer << 2) | c;
			}
			for(j=0;j<sdb->kmer_size;j++){
				c = (~bits2bit(as_array_u64list(sdb->rd_seqs), offset + rd_len - 1 - (sdb->idx_pats[n_idx][1] * sdb->kmer_size + j))) & 0x03;
				P.kmer = (P.kmer << 2) | c;
			}
		} else {
			for(j=0;j<sdb->kmer_size;j++){
				c = bits2bit(as_array_u64list(sdb->rd_seqs), offset + sdb->idx_pats[n_idx][0] * sdb->kmer_size + j);
				P.kmer = (P.kmer << 2) | c;
			}
			for(j=0;j<sdb->kmer_size;j++){
				c = bits2bit(as_array_u64list(sdb->rd_seqs), offset + sdb->idx_pats[n_idx][1] * sdb->kmer_size + j);
				P.kmer = (P.kmer << 2) | c;
			}
		}
		if(sdb->is_fast){
			k = get_u64list(sdb->fast_indexs[n_idx], P.kmer);
			if(k >> 63){
				k = (1LLU << 32) | i;
				push_u32list(sdb->links[n_idx], i);
			} else {
				push_u32list(sdb->links[n_idx], k & 0xFFFFFFFFU);
				k = (((k >> 32) + 1) << 32) | i;
			}
			set_u64list(sdb->fast_indexs[n_idx], P.kmer, k);
		} else {
			p = prepare_sr_prmhash(sdb->indexs[n_idx], P, &exists);
			if(exists){
				p->cnt ++;
				push_u32list(sdb->links[n_idx], p->rid);
			} else {
				p->kmer = P.kmer;
				p->cnt  = 1;
				push_u32list(sdb->links[n_idx], i);
			}
			p->rid = i;
		}
	}
}

void sr_index2_sdb(SR_SeqDB *sdb, uint32_t n_idx){
	SR_Primer P, *p;
	uint32_t i, j, s, rd_len, fix, beg;
	uint64_t offset, k, k1, k2;
	int exists;
	if(sdb->idx_pats[n_idx][0] != 0) return;
	if(sdb->idx_pats[n_idx][2] != 0) return;
	P.rid = 0;
	if(sdb->is_fast){
		memset(ref_u64list(sdb->fast_indexs[n_idx], 0), 0xFFU, (1LLU << ((sdb->kmer_size << 1) << 1)) * sizeof(uint64_t));
	}
	if((sdb->idx_pats[n_idx][1] + 1) * sdb->kmer_size > sdb->max_rd_len) return;
	fix = sdb->max_rd_len - (sdb->idx_pats[n_idx][1] + 1) * sdb->kmer_size + 1;
	beg = sdb->idx_pats[n_idx][1] * sdb->kmer_size;
	for(i=0;i<sdb->n_rd;i++){
		offset = sr_rdseq_offset(sdb, i);
		rd_len = sr_rdseq_length(sdb, i);
		k1 = k2 = 0;
		for(s=0;s+1<sdb->kmer_size;s++){
			k1 = (k1 << 2) | bits2bit(ref_u64list(sdb->rd_seqs, 0), offset + s);
			k2 = (k2 << 2) | bits2bit(ref_u64list(sdb->rd_seqs, 0), offset + s + beg);
		}
		for(j=0;j<fix;j++){
			if(j+(sdb->idx_pats[n_idx][1]+1)*sdb->kmer_size>rd_len){
				push_u32list(sdb->links[n_idx], i * fix + j);
				continue;
			}
			k1 = ((k1 << 2) | bits2bit(ref_u64list(sdb->rd_seqs, 0), offset + s + j)) & sdb->kmer_mask;
			k2 = ((k2 << 2) | bits2bit(ref_u64list(sdb->rd_seqs, 0), offset + s + j + beg)) & sdb->kmer_mask;
			P.kmer = (k1 << (sdb->kmer_size << 1)) | k2;
			if(sdb->is_fast){
				k = get_u64list(sdb->fast_indexs[n_idx], P.kmer);
				if(k >> 63){
					k = (1LLU << 32) | i;
					push_u32list(sdb->links[n_idx], i * fix + j);
				} else {
					push_u32list(sdb->links[n_idx], k & 0xFFFFFFFFU);
					k = (((k >> 32) + 1) << 32) | (i * fix + j);
				}
				set_u64list(sdb->fast_indexs[n_idx], P.kmer, k);
			} else {
				p = prepare_sr_prmhash(sdb->indexs[n_idx], P, &exists);
				if(exists){
					p->cnt ++;
					push_u32list(sdb->links[n_idx], p->rid);
				} else {
					p->kmer = P.kmer;
					p->cnt  = 1;
					push_u32list(sdb->links[n_idx], i * fix + j);
				}
				p->rid = i * fix + j;
			}
		}
	}
}

void sr_aln_load_rd(SR_SeqDB *sdb, SR_AlnAux *aux, uint32_t rid, uint64_t *seqs, uint64_t seqoff, uint32_t seqlen){
	uint32_t i, kmer;
	aux->rid = rid;
	aux->rd_len = seqlen;
	for(i=0;i<aux->rd_len;i+=32){ aux->rd_seq[i>>5] = sub32seqbits(seqs, seqoff + i); }
	for(i=0;i<aux->rd_len;i++){
		if((i & 0x1F) == 0) aux->rv_seq[i>>5] = 0;
		aux->rv_seq[i>>5] |= ((~bits2bit(aux->rd_seq, (aux->rd_len - 1 - i))) & 0x03) << (((~i) & 0x1F) << 1);
	}
	clear_u32list(aux->fkmers);
	kmer = 0;
	for(i=0;i<aux->rd_len;i++){
		kmer = (kmer << 2) | bits2bit(aux->rd_seq, i);
		if(i >= sdb->kmer_size - 1){
			push_u32list(aux->fkmers, kmer & sdb->kmer_mask);
		}
	}
	clear_u32list(aux->rkmers);
	kmer = 0;
	for(i=0;i<aux->rd_len;i++){
		kmer = (kmer << 2) | ((~bits2bit(aux->rd_seq, aux->rd_len - 1 - i)) & 0x03);
		if(i >= sdb->kmer_size - 1){
			push_u32list(aux->rkmers, kmer & sdb->kmer_mask);
		}
	}
	clear_sr_ptrv(aux->ptrs);
}

uint32_t sr_aln_core_simp(SR_SeqDB *sdb, SR_AlnAux *aux, SR_AlnPtr *ptr){
	SR_AlnHit *hit;
	uint64_t seq1[32], seq2[32], rd_off2;
	uint32_t i, j, rd_len1, rd_len2, len;
	uint32_t max, mm;
	rd_len2 = sr_rdseq_length(sdb, ptr->ptr);
	rd_off2 = sr_rdseq_offset(sdb, ptr->ptr);
	rd_len1 = aux->rd_len;
	if(ptr->dir1){
		for(i=ptr->off1,j=0;i<rd_len1;i+=32,j++){
			seq1[j] = sub32seqbits(aux->rv_seq, i);
		}
	} else {
		for(i=ptr->off1,j=0;i<rd_len1;i+=32,j++){
			seq1[j] = sub32seqbits(aux->rd_seq, i);
		}
	}
	if(ptr->dir2){
		for(i=ptr->off2,j=0;i<rd_len2;i+=32,j++){
			seq2[j] = dna_rev_seq(sub32seqbits(ref_u64list(sdb->rd_seqs, 0), rd_len2 + rd_off2 - i - 32), 32);
		}
	} else {
		for(i=ptr->off2,j=0;i<rd_len2;i+=32,j++){
			seq2[j] = sub32seqbits(ref_u64list(sdb->rd_seqs, 0), i + rd_off2);
		}
	}
	rd_len1 -= ptr->off1;
	rd_len2 -= ptr->off2;
	len = (rd_len1 < rd_len2)? rd_len1 : rd_len2;
	mm = 0;
	max = len * (1 -  sdb->min_similarity) + 0.5;
	if(max > sdb->max_mismatch) max = sdb->max_mismatch;
	for(i=0;i+32<=len;i+=32){
		mm += count_ones_bit64(dna_xor2ones(seq1[i>>5] ^ seq2[i>>5]));
		if(mm > max) return 0;
	}
	if(i < len){
		mm += count_ones_bit64(dna_xor2ones((seq1[i>>5] ^ seq2[i>>5]) >> ((32 - (len - i)) << 1)));
	}
	if(mm > max) return 0;
	if(ptr->off1 >= ptr->off2){
		hit = next_ref_sr_hitv(aux->hits);
		hit->n_ol    = len;
		hit->n_mm    = mm;
		hit->is_gap  = 0;
		hit->rid1    = aux->rid;
		hit->rid2    = ptr->ptr;
		hit->dir1    = ptr->dir1;
		hit->dir2    = ptr->dir2;
		hit->off     = ptr->off1 - ptr->off2;
		hit->n_cigar = append_cigars(hit->cigars, 0, ALN_CIGAR_TYPE_INS, hit->off);
		hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_MAT, len);
		if(len == rd_len2) hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_INS, rd_len1 - len);
		else hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_DEL, rd_len2 - len);
	} else {
		hit = next_ref_sr_hitv(aux->hits);
		hit->n_ol    = len;
		hit->n_mm    = mm;
		hit->is_gap  = 0;
		hit->rid1    = ptr->ptr;
		hit->rid2    = aux->rid;
		hit->dir1    = ptr->dir2;
		hit->dir2    = ptr->dir1;
		hit->off     = ptr->off2 - ptr->off1;
		hit->n_cigar = append_cigars(hit->cigars, 0, ALN_CIGAR_TYPE_INS, hit->off);
		hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_MAT, len);
		if(len == rd_len1) hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_INS, rd_len2 - len);
		else hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_DEL, rd_len1 - len);
	}
	return 1;
}

static inline uint64_t sr_stdaln_lh3(AlnCigar* cigars, AlnParam *ap, uint64_t *seqs1, uint64_t *seqs2, uint64_t off1, uint64_t off2, uint32_t dir1, uint32_t dir2, uint32_t len1, uint32_t len2){
	AlnAln *aln;
	char seq1[MAX_RD_LEN+1], seq2[MAX_RD_LEN+1];
	path_t *p;
	int i, n_cigar, mm, mn;
	if(dir1) bits2revseq(seq1, seqs1, off1, len1);
	else bits2seq(seq1, seqs1, off1, len1);
	if(dir2) bits2revseq(seq2, seqs2, off2, len2);
	else bits2seq(seq2, seqs2, off2, len2);
	aln = aln_stdaln_aux(seq1, seq2, ap, 1, 0, len1, len2);
	p = aln->path + aln->path_len - 1;
	n_cigar = -1;
	mm = mn = 0;
	while(p >= aln->path){
		if(p->ctype == FROM_M){
			mn ++;
			if(seq1[p->i] != seq2[p->j]) mm ++;
		}
		if(n_cigar >= 0 && p->ctype == cigars[n_cigar].type) cigars[n_cigar].len ++;
		else {
			n_cigar ++;
			cigars[n_cigar].type = p->ctype;
			cigars[n_cigar].len  = 1;
		}
		p --;
	}
	n_cigar ++;
	for(i=0;i<n_cigar;i++){
		switch(cigars[i].type){
			case FROM_M: cigars[i].type = ALN_CIGAR_TYPE_MAT; break;
			case FROM_I: cigars[i].type = ALN_CIGAR_TYPE_DEL; break;
			case FROM_D: cigars[i].type = ALN_CIGAR_TYPE_INS; break;
		}
	}
	//fprintf(stdout, "%s\n%s\n%s\n", aln->out1, aln->outm, aln->out2);
	//fflush(stdout);
	aln_free_AlnAln(aln);
	return ((uint64_t)n_cigar) | (((uint64_t)mm) << 10) | (((uint64_t)mn) << 30);
}

int sr_aln_core_dp(SR_SeqDB *sdb, SR_AlnAux *aux, SR_AlnPtr *ptr){
	SR_AlnHit *hit;
	uint64_t ret;
	uint32_t n_cigar, mm, ma;
	ret = sr_stdaln_lh3(aux->cigars, &sdb->ap, aux->rd_seq, ref_u64list(sdb->rd_seqs, 0), 0, sr_rdseq_offset(sdb, ptr->ptr), ptr->dir1, ptr->dir2, aux->rd_len, sr_rdseq_length(sdb, ptr->ptr));
	n_cigar = ret & 0x3FF;
	mm = (ret >> 10) & 0xFFFFF;
	ma = (ret >> 30) & 0xFFFFF;
	if(n_cigar == 0 || n_cigar > MAX_N_CIGAR) return 0;
	if(ma < sdb->min_overlap || mm > sdb->max_mismatch || mm > (uint32_t)(ma * (1 - sdb->min_similarity) + 0.5)) return 0;
	if(aux->cigars[0].type == ALN_CIGAR_TYPE_DEL){
		if(ptr->fix) return 0;
		hit = next_ref_sr_hitv(aux->hits);
		hit->n_ol    = ma;
		hit->n_mm    = mm;
		hit->is_gap  = 1;
		hit->rid1    = ptr->ptr;
		hit->rid2    = aux->rid;
		hit->dir1    = ptr->dir2;
		hit->dir2    = ptr->dir1;
		hit->off     = aux->cigars[0].len;
		memcpy(hit->cigars, aux->cigars, n_cigar * sizeof(AlnCigar));
		hit->n_cigar = n_cigar;
		flip_cigars(hit->cigars, hit->n_cigar);
	} else {
		hit = next_ref_sr_hitv(aux->hits);
		hit->n_ol    = ma;
		hit->n_mm    = mm;
		hit->is_gap  = 1;
		hit->rid1    = aux->rid;
		hit->rid2    = ptr->ptr;
		hit->dir1    = ptr->dir1;
		hit->dir2    = ptr->dir2;
		hit->off     = (aux->cigars[0].type == ALN_CIGAR_TYPE_MAT)? 0 : aux->cigars[0].len;
		memcpy(hit->cigars, aux->cigars, n_cigar * sizeof(AlnCigar));
		hit->n_cigar = n_cigar;
	}
	return 1;
}

uint32_t sr_aln_rd_core(SR_SeqDB *sdb, SR_AlnAux *aux, SR_AlnPtr *ptr){
	if(aux->rid == ptr->ptr && ptr->dir1 == ptr->dir2) return 0;
	if(sr_aln_core_simp(sdb, aux, ptr)) return 1;
	if(sdb->allow_gap && sr_aln_core_dp(sdb, aux, ptr)) return 1;
	return 0;
}

int cmp_pids_func(const void *e1, const void *e2){
	uint32_t a, b;
	a = *((uint64_t*)e1);
	b = *((uint64_t*)e2);
	if(a == b) return 0;
	if(a < b) return 1;
	return -1;
}

uint32_t _query_index_and_add_ptrv(SR_SeqDB *sdb, SR_AlnAux *aux, uint32_t idx, SR_Primer PRM, uint32_t off1, uint32_t off2, uint32_t dir1, uint32_t dir2, uint32_t fix){
	SR_AlnPtr *ptr;
	SR_Primer *prm;
	uint64_t k;
	uint32_t e, r, c;
	if(sdb->is_fast){
		k = get_u64list(sdb->fast_indexs[idx], PRM.kmer);
		if(k >> 63){
			e = 0; r = c = 0;
		} else {
			e = 1; r = k & 0xFFFFFFFFU; c = k >> 32;
		}
	} else {
		prm = get_sr_prmhash(sdb->indexs[idx], PRM);
		if(prm){
			e = 1; r = prm->rid; c = prm->cnt;
		} else {
			e = 0; r = c = 0;
		}
	}
	if(e){
		ptr = next_ref_sr_ptrv(aux->ptrs);
		ptr->n_idx  = idx;
		ptr->fix    = fix;
		ptr->dir1   = dir1;
		ptr->dir2   = dir2;
		ptr->off1   = off1;
		ptr->off2   = off2;
		ptr->ptr    = r;
		ptr->cnt    = c;
	}
	return c;
}

uint32_t sr_align_sdb(SR_SeqDB *sdb, uint32_t rid, SR_AlnAux *aux){
	SR_Primer PRM;
	SR_AlnPtr *ptr;
	uint32_t i, j, len1, len2, n_pit, set;
	sr_aln_load_rd(sdb, aux, rid, ref_u64list(sdb->rd_seqs, 0), sr_rdseq_offset(sdb, rid), sr_rdseq_length(sdb, rid));
	PRM.rid = 0;
	PRM.cnt = 0;
	n_pit = 0;
	for(j=0;j+2*sdb->kmer_size<=aux->rd_len&&j+sdb->min_overlap<=aux->rd_len;j++){
		for(i=0;i<sdb->n_idx;i++){
			if(j + (sdb->idx_pats[i][1] + 1) * sdb->kmer_size > aux->rd_len) continue;
			len1 = sdb->idx_pats[i][0] * sdb->kmer_size;
			len2 = sdb->idx_pats[i][1] * sdb->kmer_size;
			if((sdb->strand_type & 0x01)){
				if(sdb->idx_pats[i][2] == 0){
					PRM.kmer = (((uint64_t)get_u32list(aux->fkmers, j + len1)) << (sdb->kmer_size << 1)) | get_u32list(aux->fkmers, j + len2);
					n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, j, 0, 0, 0, 0);
				} else {
					PRM.kmer = (((uint64_t)get_u32list(aux->rkmers, j + len1)) << (sdb->kmer_size << 1)) | get_u32list(aux->rkmers, j + len2);
					n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, j, 0, 1, 1, 0);
				}
			}
			if(sdb->strand_type & 0x02){
				if(sdb->idx_pats[i][2] == 0){
					PRM.kmer = (((uint64_t)get_u32list(aux->rkmers, j + len1)) << (sdb->kmer_size << 1)) | get_u32list(aux->rkmers, j + len2);
					n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, j, 0, 1, 0, 0);
				} else {
					PRM.kmer = (((uint64_t)get_u32list(aux->fkmers, j + len1)) << (sdb->kmer_size << 1)) | get_u32list(aux->fkmers, j + len2);
					n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, j, 0, 0, 1, 0);
				}
			}
		}
		if(n_pit > sdb->max_hits_per_query * 4) break;
	}
	n_pit = 0;
#define SORT_UNIQ
#ifdef SORT_UNIQ
	clear_u64list(aux->pids);
	for(i=0;i<count_sr_ptrv(aux->ptrs);i++){
		ptr = ref_sr_ptrv(aux->ptrs, i);
		if(sdb->filter_dup){
			while(ptr->ptr > aux->rid){
				push_u64list(aux->pids, (((uint64_t)i) << 32) | ptr->ptr);
				rid = ptr->ptr;
				ptr->ptr = get_u32list(sdb->links[ptr->n_idx], ptr->ptr);
				if(rid == ptr->ptr) break;
			}
		} else {
			while(1){
				push_u64list(aux->pids, (((uint64_t)i) << 32) | ptr->ptr);
				rid = ptr->ptr;
				ptr->ptr = get_u32list(sdb->links[ptr->n_idx], ptr->ptr);
				if(rid == ptr->ptr) break;
			}
		}
	}
	qsort(ref_u64list(aux->pids, 0), count_u64list(aux->pids), sizeof(uint64_t), cmp_pids_func);
	rid = aux->rid;
	for(i=0;i<count_u64list(aux->pids);i++){
		ptr = ref_sr_ptrv(aux->ptrs, get_u64list(aux->pids, i) >> 32);
		ptr->ptr = get_u64list(aux->pids, i) & 0xFFFFFFFFU;
		if(ptr->ptr == rid) continue;
		rid = ptr->ptr;
		if(sr_aln_rd_core(sdb, aux, ptr)){
			n_pit ++;
			if(n_pit == sdb->max_hits_per_query) return n_pit;
		}
	}
	set = 0;
#else
	for(i=0;i<count_sr_ptrv(aux->ptrs);i++) ref_sr_ptrv(aux->ptrs, i)->num = 0;
	while(1){
		rid = 0;
		set = 0;
		for(i=0;i<count_sr_ptrv(aux->ptrs);i++){
			ptr = ref_sr_ptrv(aux->ptrs, i);
			if(ptr->num == ptr->cnt) continue;
			if(set == 0 || ptr->ptr > rid) { rid = ptr->ptr; set = 1; }
		}
		if(set == 0) break;
		set = 0;
		for(i=0;i<count_sr_ptrv(aux->ptrs);i++){
			ptr = ref_sr_ptrv(aux->ptrs, i);
			if(ptr->num == ptr->cnt) continue;
			if(sdb->filter_dup && ptr->ptr <= aux->rid){ ptr->num = ptr->cnt; continue; }
			if(ptr->ptr == rid){
				if(set == 0){
					set = 1;
					if(sr_aln_rd_core(sdb, aux, ptr)){
						n_pit ++;
						if(n_pit == sdb->max_hits_per_query) return n_pit;
					}
				}
				ptr->ptr = get_u32list(sdb->links[ptr->n_idx], ptr->ptr);
				if(ptr->ptr == rid){
					ptr->num = ptr->cnt;
				} else {
					ptr->num ++;
				}
			}
		}
	}
#endif
	return n_pit;
}

uint32_t sr_align_long_sdb(SR_SeqDB *sdb, SR_AlnAux *aux, uint32_t lrid, uint64_t *seqs, uint64_t seqoff, uint32_t seqlen){
	SR_Primer PRM;
	SR_AlnPtr *ptr;
	uint32_t i, j, len1, len2, n_pit;
	uint32_t rid;
	sr_aln_load_rd(sdb, aux, lrid, seqs, seqoff, seqlen);
	PRM.rid = 0;
	PRM.cnt = 0;
	n_pit = 0;
	for(j=0;j+2*sdb->kmer_size<=aux->rd_len&&j+sdb->max_rd_len<=aux->rd_len;j++){
		for(i=0;i<sdb->n_idx;i++){
			if(j + (sdb->idx_pats[i][1] + 1) * sdb->kmer_size > aux->rd_len) continue;
			if(j < sdb->idx_pats[i][0] * sdb->kmer_size) continue;
			len1 = sdb->idx_pats[i][0] * sdb->kmer_size;
			len2 = sdb->idx_pats[i][1] * sdb->kmer_size;
			if((sdb->strand_type & 0x01)){
				if(sdb->idx_pats[i][2] == 0){
					PRM.kmer = (((uint64_t)get_u32list(aux->fkmers, j + len1)) << (sdb->kmer_size << 1)) | get_u32list(aux->fkmers, j + len2);
					n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, j, 0, 0, 0, 1);
				} else {
					PRM.kmer = (((uint64_t)get_u32list(aux->rkmers, j + len1)) << (sdb->kmer_size << 1)) | get_u32list(aux->rkmers, j + len2);
					n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, j, 0, 1, 1, 1);
				}
			}
			if(sdb->strand_type & 0x02){
				if(sdb->idx_pats[i][2] == 0){
					PRM.kmer = (((uint64_t)get_u32list(aux->rkmers, j + len1)) << (sdb->kmer_size << 1)) | get_u32list(aux->rkmers, j + len2);
					n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, j, 0, 1, 0, 1);
				} else {
					PRM.kmer = (((uint64_t)get_u32list(aux->fkmers, j + len1)) << (sdb->kmer_size << 1)) | get_u32list(aux->fkmers, j + len2);
					n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, j, 0, 0, 1, 1);
				}
			}
		}
		if(n_pit > sdb->max_hits_per_query * 4) break;
	}
	n_pit = 0;
	clear_u64list(aux->pids);
	for(i=0;i<count_sr_ptrv(aux->ptrs);i++){
		ptr = ref_sr_ptrv(aux->ptrs, i);
		while(1){
			push_u64list(aux->pids, (((uint64_t)i) << 32) | ptr->ptr);
			rid = ptr->ptr;
			ptr->ptr = get_u32list(sdb->links[ptr->n_idx], ptr->ptr);
			if(rid == ptr->ptr) break;
		}
	}
	qsort(ref_u64list(aux->pids, 0), count_u64list(aux->pids), sizeof(uint64_t), cmp_pids_func);
	rid = 0xFFFFFFFFU;
	for(i=0;i<count_u64list(aux->pids);i++){
		ptr = ref_sr_ptrv(aux->ptrs, get_u64list(aux->pids, i) >> 32);
		ptr->ptr = get_u64list(aux->pids, i) & 0xFFFFFFFFU;
		if(ptr->ptr == rid) continue;
		rid = ptr->ptr;
		if(sr_aln_rd_core(sdb, aux, ptr)){
			n_pit ++;
			if(n_pit == sdb->max_hits_per_query) return n_pit;
		}
	}
	return n_pit;
}

uint32_t sr_align2_sdb(SR_SeqDB *sdb, uint32_t rid, SR_AlnAux *aux){
	SR_Primer PRM;
	SR_AlnPtr *ptr;
	uint32_t i, j, n_pit, fix;
	sr_aln_load_rd(sdb, aux, rid, ref_u64list(sdb->rd_seqs, 0), sr_rdseq_offset(sdb, rid), sr_rdseq_length(sdb, rid));
	PRM.rid = 0;
	PRM.cnt = 0;
	n_pit = 0;
	for(i=0;i<sdb->n_idx;i++){
		if(sdb->idx_pats[i][0] || sdb->idx_pats[i][2]) continue;
		for(j=0;j+sdb->idx_pats[i][1]<sdb->n_seed;j++){
			if((j + sdb->idx_pats[i][1] + 1) * sdb->kmer_size > aux->rd_len) continue;
			if(sdb->strand_type & 0x01){
				PRM.kmer = (((uint64_t)get_u32list(aux->fkmers, j * sdb->kmer_size)) << (sdb->kmer_size << 1)) |
					get_u32list(aux->fkmers, (j + sdb->idx_pats[i][1]) * sdb->kmer_size);
				n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, j * sdb->kmer_size, 0, 0, 0, 0);
				PRM.kmer = (((uint64_t)get_u32list(aux->fkmers, aux->rd_len - sdb->kmer_size - (j + sdb->idx_pats[i][1]) * sdb->kmer_size)) << (sdb->kmer_size << 1)) |
					get_u32list(aux->fkmers, aux->rd_len - sdb->kmer_size - j * sdb->kmer_size);
				n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, aux->rd_len - sdb->kmer_size - (j + sdb->idx_pats[i][1]) * sdb->kmer_size, 0, 0, 0, 0);
			}
			if(sdb->strand_type & 0x02){
				PRM.kmer = (((uint64_t)get_u32list(aux->rkmers, j * sdb->kmer_size)) << (sdb->kmer_size << 1)) |
					get_u32list(aux->rkmers, (j + sdb->idx_pats[i][1]) * sdb->kmer_size);
				n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, j * sdb->kmer_size, 0, 1, 0, 0);
				PRM.kmer = (((uint64_t)get_u32list(aux->rkmers, aux->rd_len - sdb->kmer_size - (j + sdb->idx_pats[i][1]) * sdb->kmer_size)) << (sdb->kmer_size << 1)) |
					get_u32list(aux->rkmers, aux->rd_len - sdb->kmer_size - j * sdb->kmer_size);
				n_pit += _query_index_and_add_ptrv(sdb, aux, i, PRM, aux->rd_len - sdb->kmer_size - (j + sdb->idx_pats[i][1]) * sdb->kmer_size, 0, 1, 0, 0);
			}
		}
		if(n_pit > sdb->max_hits_per_query * 4) break;
	}
	n_pit = 0;
	clear_u64list(aux->pids);
	for(i=0;i<count_sr_ptrv(aux->ptrs);i++){
		ptr = ref_sr_ptrv(aux->ptrs, i);
		if(sdb->filter_dup){
			fix = sdb->max_rd_len - (sdb->idx_pats[ptr->n_idx][1] + 1) * sdb->kmer_size + 1;
			while((ptr->ptr / fix) > aux->rid){
				push_u64list(aux->pids, (((uint64_t)i) << 32) | ptr->ptr);
				rid = ptr->ptr;
				ptr->ptr = get_u32list(sdb->links[ptr->n_idx], ptr->ptr);
				if(rid == ptr->ptr) break;
			}
		} else {
			while(1){
				push_u64list(aux->pids, (((uint64_t)i) << 32) | ptr->ptr);
				rid = ptr->ptr;
				ptr->ptr = get_u32list(sdb->links[ptr->n_idx], ptr->ptr);
				if(rid == ptr->ptr) break;
			}
		}
	}
	qsort(ref_u64list(aux->pids, 0), count_u64list(aux->pids), sizeof(uint64_t), cmp_pids_func);
	rid = aux->rid;
	for(i=0;i<count_u64list(aux->pids);i++){
		ptr = ref_sr_ptrv(aux->ptrs, get_u64list(aux->pids, i) >> 32);
		fix = sdb->max_rd_len - (sdb->idx_pats[ptr->n_idx][1] + 1) * sdb->kmer_size + 1;
		ptr->ptr = get_u64list(aux->pids, i) & 0xFFFFFFFFU;
		if((ptr->ptr / fix) == rid) continue;
		rid = ptr->ptr/ fix;
		ptr->off2 = ptr->ptr % fix;
		ptr->ptr = rid;
		if(sr_aln_rd_core(sdb, aux, ptr)){
			n_pit ++;
			if(n_pit == sdb->max_hits_per_query) return n_pit;
		}
	}
	return n_pit;
}

void sr_output_hit(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out){
	char str[MAX_RD_LEN + 1];
	cigars2string(hit->cigars, hit->n_cigar, str);
	fprintf(out, "%u\t%c\t%u\t%c\t%s", hit->rid1, "+-"[hit->dir1], hit->rid2, "+-"[hit->dir2], str);
	if(hit->rid1 < sdb->n_rd){
		if(hit->dir1) bits2revseq(str, as_array_u64list(sdb->rd_seqs), sr_rdseq_offset(sdb, hit->rid1), sr_rdseq_length(sdb, hit->rid1));
		else bits2seq(str, as_array_u64list(sdb->rd_seqs), sr_rdseq_offset(sdb, hit->rid1), sr_rdseq_length(sdb, hit->rid1));
		fprintf(out, "\t%s", str);
		if(hit->dir2) bits2revseq(str, as_array_u64list(sdb->rd_seqs), sr_rdseq_offset(sdb, hit->rid2), sr_rdseq_length(sdb, hit->rid2));
		else bits2seq(str, as_array_u64list(sdb->rd_seqs), sr_rdseq_offset(sdb, hit->rid2), sr_rdseq_length(sdb, hit->rid2));
		fprintf(out, "\t%s", str);
	}
	fprintf(out, "\n");
}

