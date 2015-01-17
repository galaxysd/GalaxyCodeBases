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

SR_AlnAux* sr_init_aux(){
	SR_AlnAux *aux;
	aux = malloc(sizeof(SR_AlnAux));
	aux->hits   = init_sr_hitv(32);
	aux->ptrs   = init_sr_ptrv(32);
	aux->hits_per_off = init_u32list(128);
	aux->pids   = init_u64list(1024);
	aux->fkmers = init_u32list(32);
	aux->rkmers = init_u32list(32);
	aux->rd_flags  = init_bitvec(128);
	aux->rd_tested = init_u32list(32);
	aux->auto_clear = 1;
	aux->strand_type = 0x03;
	aux->max_ol      = MAX_RD_LEN;
	aux->min_ol      = 1;
	aux->min_sm      = 0.7;
	aux->max_hits    = 1024;
	aux->max_hits_per_off  = 1024;
	aux->hit_id_range[0] = 0;
	aux->hit_id_range[1] = 0xFFFFFFFFU;
	aux->best_ols[0] = 0;
	aux->best_ols[1] = 0;
	aux->allow_gap = 0;
	aux->best_mode = 0;
	aux->rec_cigar = 1;
	return aux;
}

void sr_fit_aux2sdb(SR_AlnAux *aux, SR_SeqDB *sdb){
	encap_bitvec(aux->rd_flags, sdb->n_rd);
}

void sr_clear_aux(SR_AlnAux *aux){
	clear_sr_hitv(aux->hits);
}

void sr_block_aux(SR_AlnAux *aux, uint32_t rid){
	one_bitvec(aux->rd_flags, rid);
	push_u32list(aux->rd_tested, rid);
}

void sr_block_forever_aux(SR_AlnAux *aux, uint32_t rid){
	one_bitvec(aux->rd_flags, rid);
}

void sr_clear_rd_flags_aux(SR_AlnAux *aux){
	uint32_t i;
	for(i=0;i<aux->rd_tested->size;i++) zero_bitvec(aux->rd_flags, get_u32list(aux->rd_tested, i));
	clear_u32list(aux->rd_tested);
	zeros_u32list(aux->hits_per_off);
}

void sr_free_aux(SR_AlnAux *aux){
	free_sr_hitv(aux->hits);
	free_sr_ptrv(aux->ptrs);
	free_u64list(aux->pids);
	free_u32list(aux->fkmers);
	free_u32list(aux->rkmers);
	free_u32list(aux->hits_per_off);
	free_bitvec(aux->rd_flags);
	free_u32list(aux->rd_tested);
	free(aux);
}

SR_SeqDB* sr_init2_sdb(uint32_t kmer_size, uint32_t seed_kmer_num, uint32_t rd_len, uint32_t min_cpx, float hash_factor){
	SR_SeqDB *sdb;
	uint32_t i, j;
	if(kmer_size > MAX_KMER_SIZE) kmer_size = MAX_KMER_SIZE;
	if(kmer_size < MIN_KMER_SIZE) kmer_size = MIN_KMER_SIZE;
	sdb = malloc(sizeof(SR_SeqDB));
	sdb->rdseqs = init_basebank();
	encap_basebank(sdb->rdseqs, SR_RDSEQ_PADDING);
	sdb->rdseqs->size = SR_RDSEQ_PADDING;
	sdb->rd_filtered = init_bitvec(128);
	sdb->rd_len  = rd_len;
	sdb->n_rd    = 0;
	sdb->kmer_size = kmer_size;
	sdb->kmer_mask = (1LLU << (kmer_size << 1)) - 1U;
	sdb->min_complexity = min_cpx;
	sdb->ap             = (AlnParam){10, 2, 2, aln_sm_nt, 16, 75};
	sdb->hash_factor = hash_factor;

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
	sdb->indexs = malloc(sizeof(sr_prm_hash*) * sdb->n_idx);
	sdb->links  = malloc(sizeof(u32list*) * sdb->n_idx);
	for(i=0;i<sdb->n_idx;i++){
		sdb->indexs[i]  = NULL;
		sdb->links[i]   = init_u32list(64);
	}

	return sdb;
}

void sr_mask_seed_sdb(SR_SeqDB *sdb, uint32_t seed_idx){
	uint32_t i, j, idx;
	idx = 0;
	for(i=0;i<sdb->n_seed-1;i++){
		for(j=i+1;j<sdb->n_seed;j++){
			if(i == seed_idx || j == seed_idx){
				sdb->idx_pats[idx + 0][2] = 2;
				sdb->idx_pats[idx + 1][2] = 2;
			}
			idx += 2;
		}
	}
}

size_t sr_dump_sdb(SR_SeqDB *sdb, FILE *out){
	size_t i, size;
	size = 0;
	size += fwrite(sdb, sizeof(SR_SeqDB), 1, out) * sizeof(SR_SeqDB);
	size += dump_basebank(sdb->rdseqs, out);
	size += dump_bitvec(sdb->rd_filtered, out);
	for(i=0;i<sdb->n_idx;i++){
		size += dump_sr_prm_hash(sdb->indexs[i], out);
		size += dump_u32list(sdb->links[i], out);
	}
	return size;
}

SR_SeqDB* sr_load_sdb(FILE *inp){
	SR_SeqDB *sdb;
	size_t size, i;
	sdb = (SR_SeqDB*)malloc(sizeof(SR_SeqDB));
	if((size = fread(sdb, sizeof(SR_SeqDB), 1, inp)) != 1) return NULL;
	if((sdb->rdseqs = load_basebank(inp)) == NULL) return NULL;
	if((sdb->rd_filtered = load_bitvec(inp)) == NULL) return NULL;
	sdb->indexs = malloc(sizeof(sr_prmhash*) * sdb->n_idx);
	sdb->links  = malloc(sizeof(u32list*) * sdb->n_idx);
	for(i=0;i<sdb->n_idx;i++){
		if((sdb->indexs[i] = load_sr_prm_hash(inp)) == NULL) return NULL;
		if((sdb->links[i]  = load_u32list(inp)) == NULL) return NULL;
	}
	return sdb;
}

SR_SeqDB* sr_init_sdb(uint32_t kmer_size, uint32_t seed_kmer_num, uint32_t rd_len, uint32_t min_cpx){
	return sr_init2_sdb(kmer_size, seed_kmer_num, rd_len, min_cpx, 0.77);
}

uint32_t cal_seq_complexity(uint64_t *seq, uint64_t off, uint16_t len){
	uint32_t i, tri;
	uint64_t tris;
	tris = 0;
	tri = 0;
	tri = (tri << 2) | bits2bit(seq, off);
	tri = (tri << 2) | bits2bit(seq, off + 1);
	for(i=2;i<len;i++){
		tri = (tri << 2) | bits2bit(seq, off + i);
		tris |= 1LLU << (tri & 0x3F);
	}
	return count_ones_bit64(tris);
}

int sr_push_sdb(SR_SeqDB *sdb, char *seq, uint8_t seqlen){
	if(sdb->rd_len){
		if(sdb->rd_len != seqlen){
			fprintf(stderr, " -- %d bp != %d bp, \"%s\" in %s -- %s:%d --\n", seqlen, sdb->rd_len, seq, __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
		}
	} else sdb->rd_len = seqlen;
	seq2basebank(sdb->rdseqs, seq, seqlen);
	zero2bitvec(sdb->rd_filtered);
	sdb->n_rd ++;
	return 1;
}

uint32_t sr_mask_low_complexity(SR_SeqDB *sdb, uint32_t n_part, uint32_t part_idx){
	uint64_t *seqs, off;
	uint32_t i, s, e, len, cpx, ret;
	if(sdb->min_complexity < 2) return 0;
	s = (sdb->n_rd / n_part) * part_idx;
	e = (part_idx + 1 == n_part)? sdb->n_rd : (sdb->n_rd / n_part) * (part_idx + 1);
	ret = 0;
	seqs = sdb->rdseqs->bits;
	off = sr_rdseq_offset(sdb, s);
	for(i=s;i<e;i++){
		len = sr_rdseq_length(sdb, i);
		cpx = cal_seq_complexity(seqs, off, len);
		off += len;
		if(cpx < sdb->min_complexity){
			one_bitvec(sdb->rd_filtered, i);
			ret ++;
		}
	}
	return ret;
}

void sr_ready_sdb(SR_SeqDB *sdb){
	encap_basebank(sdb->rdseqs, SR_RDSEQ_PADDING);
}

void sr_reset_sdb(SR_SeqDB *sdb){
	uint32_t i;
	for(i=0;i<sdb->n_idx;i++){
		if(sdb->indexs[i]){
			free_sr_prm_hash(sdb->indexs[i]);
			sdb->indexs[i] = NULL;
		}
		if(sdb->links[i]) clear_u32list(sdb->links[i]);
	}
	sdb->n_rd = 0;
	clear_basebank(sdb->rdseqs);
	sdb->rdseqs->size = SR_RDSEQ_PADDING;
	zeros_bitvec(sdb->rd_filtered);
}

void sr_free_sdb(SR_SeqDB *sdb){
	uint32_t i;
	free_basebank(sdb->rdseqs);
	free_bitvec(sdb->rd_filtered);
	for(i=0;i<sdb->n_idx;i++){
		if(sdb->indexs[i]) free_sr_prm_hash(sdb->indexs[i]);
		if(sdb->links[i]) free_u32list(sdb->links[i]);
	}
	free(sdb->indexs);
	free(sdb->links);
	free(sdb);
}

void sr_index_sdb(SR_SeqDB *sdb, uint32_t idx){
	SR_Primer P, *p;
	sr_prmhash *hash;
	uint32_t i, j, rd_len, c, n;
	uint64_t offset, kmer;
	int exists;
	if(idx >= sdb->n_idx) return;
	if(sdb->idx_pats[idx][2] > 1) return;
	if(sdb->indexs[idx]){
		free_sr_prm_hash(sdb->indexs[idx]);
		sdb->indexs[idx] = NULL;
	}
	free_u32list(sdb->links[idx]);
	if(sdb->rd_len < (sdb->idx_pats[idx][1] + 1) * sdb->kmer_size) return;
	hash = init2_sr_prmhash(13, sdb->hash_factor);
	P.off = P.cnt = 0;
	for(i=0;i<sdb->n_rd;i++){
		P.kmer1 = P.kmer2 = 0;
		if(get_bitvec(sdb->rd_filtered, i)){ continue; }
		offset = sr_rdseq_offset(sdb, i);
		rd_len = sr_rdseq_length(sdb, i);
		if(sdb->idx_pats[idx][2]){
			for(j=0;j<sdb->kmer_size;j++){
				c = (~bits2bit(sdb->rdseqs->bits, offset + rd_len - 1 - (sdb->idx_pats[idx][0] * sdb->kmer_size + j))) & 0x03;
				P.kmer1 = (P.kmer1 << 2) | c;
			}
			for(j=0;j<sdb->kmer_size;j++){
				c = (~bits2bit(sdb->rdseqs->bits, offset + rd_len - 1 - (sdb->idx_pats[idx][1] * sdb->kmer_size + j))) & 0x03;
				P.kmer2 = (P.kmer2 << 2) | c;
			}
		} else {
			for(j=0;j<sdb->kmer_size;j++){
				c = bits2bit(sdb->rdseqs->bits, offset + sdb->idx_pats[idx][0] * sdb->kmer_size + j);
				P.kmer1 = (P.kmer1 << 2) | c;
			}
			for(j=0;j<sdb->kmer_size;j++){
				c = bits2bit(sdb->rdseqs->bits, offset + sdb->idx_pats[idx][1] * sdb->kmer_size + j);
				P.kmer2 = (P.kmer2 << 2) | c;
			}
		}
		kmer = (((uint64_t)P.kmer1) << (sdb->kmer_size << 1)) | P.kmer2;
		p = prepare_sr_prmhash(hash, P, &exists);
		if(exists){
			p->cnt ++;
		} else {
			p->kmer1 = P.kmer1;
			p->kmer2 = P.kmer2;
			p->off   = 0;
			p->cnt   = 1;
		}
	}
	sdb->indexs[idx] = sr_prmhash2sr_prm_hash(hash);
	free_sr_prmhash(hash);
	reset_iter_sr_prm_hash(sdb->indexs[idx]);
	n = 0;
	while((p = ref_iter_sr_prm_hash(sdb->indexs[idx]))){
		p->off = n;
		n += p->cnt;
		p->cnt = 0;
	}
	sdb->links[idx] = init_u32list(n);
	sdb->links[idx]->size = n;
	for(i=0;i<sdb->n_rd;i++){
		if(get_bitvec(sdb->rd_filtered, i)){ continue; }
		P.kmer1 = P.kmer2 = 0;
		offset = sr_rdseq_offset(sdb, i);
		rd_len = sr_rdseq_length(sdb, i);
		if(sdb->idx_pats[idx][2]){
			for(j=0;j<sdb->kmer_size;j++){
				c = (~bits2bit(sdb->rdseqs->bits, offset + rd_len - 1 - (sdb->idx_pats[idx][0] * sdb->kmer_size + j))) & 0x03;
				P.kmer1 = (P.kmer1 << 2) | c;
			}
			for(j=0;j<sdb->kmer_size;j++){
				c = (~bits2bit(sdb->rdseqs->bits, offset + rd_len - 1 - (sdb->idx_pats[idx][1] * sdb->kmer_size + j))) & 0x03;
				P.kmer2 = (P.kmer2 << 2) | c;
			}
		} else {
			for(j=0;j<sdb->kmer_size;j++){
				c = bits2bit(sdb->rdseqs->bits, offset + sdb->idx_pats[idx][0] * sdb->kmer_size + j);
				P.kmer1 = (P.kmer1 << 2) | c;
			}
			for(j=0;j<sdb->kmer_size;j++){
				c = bits2bit(sdb->rdseqs->bits, offset + sdb->idx_pats[idx][1] * sdb->kmer_size + j);
				P.kmer2 = (P.kmer2 << 2) | c;
			}
		}
		kmer = (((uint64_t)P.kmer1) << (sdb->kmer_size << 1)) | P.kmer2;
		p = get_sr_prm_hash(sdb->indexs[idx], P);
		set_u32list(sdb->links[idx], p->off + p->cnt, i);
		p->cnt ++;
	}
}

void sr_clear_index_sdb(SR_SeqDB *sdb, uint32_t idx){
	if(sdb->indexs[idx]){
		free_sr_prm_hash(sdb->indexs[idx]);
		sdb->indexs[idx] = NULL;
	}
	free_u32list(sdb->links[idx]);
	sdb->links[idx] = init_u32list(64);
}

void _sr_aux_load_core(SR_AlnAux *aux, uint32_t kmer_size, uint32_t kmer_mask){
	uint32_t i, kmer;
	encap_u32list(aux->hits_per_off, aux->rd_len * 2);
	zeros_u32list(aux->hits_per_off);
	for(i=0;i<aux->rd_len;i++){
		if((i & 0x1F) == 0) aux->rv_seq[i>>5] = 0;
		aux->rv_seq[i>>5] |= ((~bits2bit(aux->rd_seq, (aux->rd_len - 1 - i))) & 0x03) << (((~i) & 0x1F) << 1);
	}
	clear_sr_ptrv(aux->ptrs);
	clear_u32list(aux->fkmers);
	kmer = 0;
	for(i=0;i<aux->rd_len;i++){
		kmer = (kmer << 2) | bits2bit(aux->rd_seq, i);
		if(i >= kmer_size - 1){
			push_u32list(aux->fkmers, kmer & kmer_mask);
		}
	}
	clear_u32list(aux->rkmers);
	kmer = 0;
	for(i=0;i<aux->rd_len;i++){
		kmer = (kmer << 2) | ((~bits2bit(aux->rd_seq, aux->rd_len - 1 - i)) & 0x03);
		if(i >= kmer_size - 1){
			push_u32list(aux->rkmers, kmer & kmer_mask);
		}
	}

}

void sr_aux_load(SR_AlnAux *aux, SR_SeqDB *sdb, uint32_t rid, uint64_t *seqs, uint64_t seqoff, uint32_t seqlen){
	uint32_t i;
	aux->rid = rid;
	aux->rd_len = seqlen;
	for(i=0;i<aux->rd_len;i+=32) aux->rd_seq[i>>5] = sub32seqbits(seqs, seqoff + i);
	_sr_aux_load_core(aux, sdb->kmer_size, sdb->kmer_mask);
}

void sr_aux_load2(SR_AlnAux *aux, SR_SeqDB *sdb, uint32_t rid, char *seq, uint32_t seqlen){
	uint32_t i;
	aux->rid = rid;
	aux->rd_len = seqlen;
	for(i=0;i<aux->rd_len;i++) bit2bits(aux->rd_seq, i, base_bit_table[(int)seq[i]]);
	_sr_aux_load_core(aux, sdb->kmer_size, sdb->kmer_mask);
}

uint32_t sr_aln_core_simp(SR_SeqDB *sdb, SR_AlnAux *aux, SR_AlnPtr *ptr){
	SR_AlnHit *hit;
	uint64_t seq1[32], seq2[32], rd_off2;
	uint32_t i, j, rd_len1, rd_len2, len;
	uint32_t max, mm;
	rd_len2 = sr_rdseq_length(sdb, ptr->ptr);
	rd_off2 = sr_rdseq_offset(sdb, ptr->ptr);
	rd_len1 = aux->rd_len;
	if(ptr->off1 <= ptr->off2){
		ptr->off2 -= ptr->off1; ptr->off1 = 0;
	} else {
		ptr->off1 -= ptr->off2; ptr->off2 = 0;
	}
	len = (rd_len1 - ptr->off1  < rd_len2 - ptr->off2)? rd_len1 - ptr->off1 : rd_len2 - ptr->off2;
	if(len < aux->min_ol) return 0;
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
			seq2[j] = dna_rev_seq(sub32seqbits(sdb->rdseqs->bits, rd_len2 + rd_off2 - i - 32), 32);
		}
	} else {
		for(i=ptr->off2,j=0;i<rd_len2;i+=32,j++){
			seq2[j] = sub32seqbits(sdb->rdseqs->bits, i + rd_off2);
		}
	}
	rd_len1 -= ptr->off1;
	rd_len2 -= ptr->off2;
	mm = 0;
	max = len * (1 - aux->min_sm) + 0.5;
	for(i=0;i+32<=len;i+=32){
		mm += count_ones_bit64(dna_xor2ones(seq1[i>>5] ^ seq2[i>>5]));
		if(mm > max) return 0;
	}
	if(i < len){
		mm += count_ones_bit64(dna_xor2ones((seq1[i>>5] ^ seq2[i>>5]) >> ((32 - (len - i)) << 1)));
	}
	if(mm > max) return 0;
	{
		hit = next_ref_sr_hitv(aux->hits);
		hit->n_ol    = len;
		hit->n_mm    = mm;
		hit->is_gap  = 0;
		hit->rid     = ptr->ptr;
		hit->dir1    = ptr->dir1;
		hit->dir2    = ptr->dir2;
		hit->off     = ptr->off1 - ptr->off2;
		hit->n_cigar = 0;
		if(aux->rec_cigar){
			if(ptr->off1 < ptr->off2){
				hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_DEL, ptr->off2 - ptr->off1);
			} else if(ptr->off1 > ptr->off2){
				hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_INS, ptr->off1 - ptr->off2);
			}
			hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_MAT, len);
			if(rd_len1 < rd_len2){
				hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_DEL, rd_len2 - rd_len1);
			} else if(rd_len1 > rd_len2){
				hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_INS, rd_len1 - rd_len2);
			}
		}
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
	aln = aln_stdaln_aux(seq1, seq2, ap, len1, len2, 0, 0);
	aln->start1 --;
	aln->start2 --;
	mm = mn = 0;
	n_cigar = 0;
	if(aln->start1 <= aln->start2){
		n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_INS, aln->start2 - aln->start1);
		n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_MAT, aln->start1);
		mn += aln->start1;
		for(i=0;i<aln->start1;i++) if(seq1[i] != seq2[i + aln->start2 - aln->start1]) mm ++;
	} else {
		n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_DEL, aln->start1 - aln->start2);
		n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_MAT, aln->start2);
		mn += aln->start2;
		for(i=0;i<aln->start2;i++) if(seq2[i] != seq1[i + aln->start1 - aln->start2]) mm ++;
	}
	p = aln->path + aln->path_len - 1;
	while(p >= aln->path){
		switch(p->ctype){
			case FROM_M: mn ++; if(seq1[p->i] != seq2[p->j]) mm ++; n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_MAT, 1); break;
			case FROM_I: n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_DEL, 1); break;
			case FROM_D: n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_INS, 1); break;
		}
		p --;
	}
	if(len1 - aln->end1 <= len2 - aln->end2){
		n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_MAT, len1 - aln->end1);
		mn += len1 - aln->end1;
		for(i=aln->end1;i<(int)len1;i++) if(seq1[i] != seq2[i + aln->end2 - aln->end1]) mm ++;
		n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_DEL, len2 - aln->end2 - (len1 - aln->end1));
	} else {
		n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_MAT, len2 - aln->end2);
		mn += len2 - aln->end2;
		for(i=aln->end2;i<(int)len2;i++) if(seq2[i] != seq1[i + aln->end1 - aln->end2]) mm ++;
		n_cigar = append_cigars(cigars, n_cigar, ALN_CIGAR_TYPE_INS, len1 - aln->end1 - (len2 - aln->end2));
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
	ret = sr_stdaln_lh3(aux->cigars, &sdb->ap, aux->rd_seq, sdb->rdseqs->bits, 0, sr_rdseq_offset(sdb, ptr->ptr), ptr->dir1, ptr->dir2, aux->rd_len, sr_rdseq_length(sdb, ptr->ptr));
	n_cigar = ret & 0x3FF;
	mm = (ret >> 10) & 0xFFFFF;
	ma = (ret >> 30) & 0xFFFFF;
	if(n_cigar == 0 || n_cigar > MAX_N_CIGAR) return 0;
	if(ptr->fix && (aux->cigars[0].type == ALN_CIGAR_TYPE_DEL || aux->cigars[n_cigar-1].type == ALN_CIGAR_TYPE_DEL)) return 0;
	if(ma < aux->min_ol || ma > aux->max_ol || mm > (uint32_t)(ma * (1 - aux->min_sm) + 0.5)) return 0;
	{
		hit = next_ref_sr_hitv(aux->hits);
		hit->n_ol    = ma;
		hit->n_mm    = mm;
		hit->is_gap  = 1;
		hit->rid     = ptr->ptr;
		hit->dir1    = ptr->dir1;
		hit->dir2    = ptr->dir2;
		hit->off     = (aux->cigars[0].type == ALN_CIGAR_TYPE_MAT)? 0 : ((aux->cigars[0].type == ALN_CIGAR_TYPE_DEL)? - aux->cigars[0].len : aux->cigars[0].len);
		if(aux->rec_cigar){
			memcpy(hit->cigars, aux->cigars, n_cigar * sizeof(AlnCigar));
			hit->n_cigar = n_cigar;
		} else hit->n_cigar = 0;
	}
	return 1;
}

uint32_t sr_aln_rd_core(SR_SeqDB *sdb, SR_AlnAux *aux, SR_AlnPtr *ptr){
	if(sr_aln_core_simp(sdb, aux, ptr)) return 1;
	if(aux->allow_gap && sr_aln_core_dp(sdb, aux, ptr)) return 1;
	return 0;
}

uint32_t sr_align_sdb(SR_SeqDB *sdb, SR_AlnAux *aux){
	SR_Primer PRM, *prm;
	SR_AlnPtr *ptr, PTR;
	SR_AlnHit *hit;
	u32list *kmers;
	uint32_t i, j, s, m, klen, len1, len2, ol, off, n_pit, dir;
	int beg, end;
	n_pit = 0;
	if(aux->rd_len < aux->min_ol || sdb->rd_len < aux->min_ol) goto RET;
	PRM.off = PRM.cnt = 0;
	ptr = &PTR;
	aux->best_cnt[0] = 0;
	aux->best_cnt[1] = 0;
	for(i=0;i<sdb->n_idx;i++){
		if(sdb->indexs[i] == NULL || sdb->indexs[i]->count == 0) continue;
		len1 = sdb->idx_pats[i][0] * sdb->kmer_size;
		len2 = sdb->idx_pats[i][1] * sdb->kmer_size;
		if((sdb->idx_pats[i][1] - sdb->idx_pats[i][0] + 1) * sdb->kmer_size > aux->rd_len) continue;
		klen = (sdb->idx_pats[i][1] - sdb->idx_pats[i][0]) * sdb->kmer_size;
		for(s=0;s<2;s++){
			if(((aux->strand_type >> s) & 0x01) == 0) continue;
			kmers = (s ^ sdb->idx_pats[i][2])? aux->rkmers : aux->fkmers;
			for(j=0;j<=aux->rd_len-(klen+sdb->kmer_size);j++){
				beg = ((int)len1) - ((int)j);
				if(beg < 0) beg = 0;
				end = ((int)len1) - ((int)j) + ((int)aux->rd_len);
				if(end > (int)sdb->rd_len) end = sdb->rd_len;
				ol = end - beg;
				if(ol > aux->max_ol || ol < aux->min_ol) continue;
				dir = (s ^ sdb->idx_pats[i][2]) ^ ((j < len1)? 1 : 0);
				if(aux->best_mode && ol < aux->best_ols[dir]) continue;
				if(aux->best_mode == 2 && ol == aux->best_ols[dir] && aux->best_cnt[dir]) continue;
				PRM.kmer1 = kmers->buffer[j];
				PRM.kmer2 = kmers->buffer[j + klen];
				prm = get_sr_prm_hash(sdb->indexs[i], PRM);
				if(prm == NULL) continue;
				ptr->idx   = i;
				ptr->fix   = 0;
				ptr->dir1  = s ^ sdb->idx_pats[i][2];
				ptr->dir2  = sdb->idx_pats[i][2];
				ptr->off1  = j;
				ptr->off2  = len1;
				if(ptr->off1 >= ptr->off2){
					off = aux->rd_len + ptr->off1 - ptr->off2;
				} else {
					off = ptr->off2 - ptr->off1;
				}
				for(m=0;m<prm->cnt;m++){
					ptr->ptr = sdb->links[ptr->idx]->buffer[m + prm->off];
					if(aux->hits_per_off->buffer[off] >= aux->max_hits_per_off) break;
					if(ptr->ptr < aux->hit_id_range[0]) continue;
					if(ptr->ptr > aux->hit_id_range[1]) break;
					if(get_bitvec(aux->rd_flags, ptr->ptr)) continue;
					one_bitvec(aux->rd_flags, ptr->ptr);
					push_u32list(aux->rd_tested, ptr->ptr);
					if(sr_aln_rd_core(sdb, aux, ptr)){
						if(aux->best_mode){
							hit = ref_sr_hitv(aux->hits, aux->hits->size - 1);
							if(hit->n_ol > aux->best_ols[dir]){
								aux->best_ols[dir] = hit->n_ol;
								aux->best_cnt[dir] = 1;
							} else {
								aux->best_cnt[dir] ++;
							}
						}
						aux->hits_per_off->buffer[off] ++;
						n_pit ++;
						if(n_pit == aux->max_hits) break;
					}
				}
				if(n_pit && n_pit == aux->max_hits) break;
			}
			if(n_pit && n_pit == aux->max_hits) break;
		}
		if(n_pit && n_pit == aux->max_hits) break;
	}
RET:
	if(aux->auto_clear) sr_clear_rd_flags_aux(aux);
	return n_pit;
}

