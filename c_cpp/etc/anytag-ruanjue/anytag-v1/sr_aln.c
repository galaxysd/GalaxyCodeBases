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

SR_SeqDB* sr_init_sdb(char *out_prefix, uint32_t n_cpu, uint8_t kmer_size, uint16_t rd_len){
	SR_SeqDB *sdb;
	uint32_t i;
	int fd;
	if(kmer_size > MAX_KMER_SIZE) kmer_size = MAX_KMER_SIZE;
	if(kmer_size < MIN_KMER_SIZE) kmer_size = MIN_KMER_SIZE;
	sdb = malloc(sizeof(SR_SeqDB));
	sdb->rd_seqs = init_u64list(SR_RDSEQ_PADDING / 32 + 32);
	for(i=0;i<SR_RDSEQ_PADDING/32;i++) set_u64list(sdb->rd_seqs, i, 0);
	sdb->offset = SR_RDSEQ_PADDING;
	sdb->rd_offs = init_u64list(16);
	sdb->rd_lens = init_u16list(16);
	sdb->rd_visuable = init_bitvec(128);
	sdb->rd_filtered = init_bitvec(128);
	sdb->rd_len  = rd_len;
	sdb->n_rd    = 0;
	sdb->prefix  = out_prefix;
	if(n_cpu < 1) n_cpu = 1;
	sdb->n_cpu   = n_cpu;

	sdb->prefix = out_prefix;
	if(sdb->prefix){
		if(strcmp(sdb->prefix, "-") == 0){
			sdb->prefix = malloc(30);
			sprintf(sdb->prefix, "tmp.ruanjue.XXXXXX");
			if((fd = mkstemp(sdb->prefix)) == -1){
				fprintf(stdout, " -- Cannot create %s in %s -- %s:%d --\n", sdb->prefix, __FUNCTION__, __FILE__, __LINE__);
				fflush(stdout);
				abort();
			}
			close(fd);
			unlink(sdb->prefix);
			sdb->out = stdout;
		} else if((sdb->out = fopen(sdb->prefix, "w+")) == NULL){
			fprintf(stdout, " -- Cannot write %s in %s -- %s:%d --\n", sdb->prefix, __FUNCTION__, __FILE__, __LINE__);
			fflush(stdout);
			abort();
		}
		sdb->hits = NULL;
	} else {
		sdb->out = NULL;
		sdb->hits = init_sr_hitv(16);
	}

	sdb->kmer_size = kmer_size;
	sdb->kmer_mask = (1LLU << (kmer_size << 1)) - 1U;
	sdb->strand_type    = 2;
	sdb->min_overlap    = kmer_size * 2;
	sdb->min_similarity = 0.9;
	sdb->max_mismatch   = 10;
	sdb->max_hits_per_query = 1024;
	sdb->min_complexity = 10;
	sdb->allow_gap      = 0;
	sdb->gap_init       = 24;
	sdb->ap             = (AlnParam){10, 2, 2, aln_sm_nt, 16, 75};

	sdb->indexs = malloc(sizeof(sr_prmhash*) * sdb->n_cpu);
	sdb->auxs   = malloc(sizeof(SR_AlnAux*) * sdb->n_cpu);
	sdb->links  = init_u32list(1024);
	for(i=0;i<sdb->n_cpu;i++){
		sdb->indexs[i]        = NULL;
		sdb->auxs[i]          = malloc(sizeof(SR_AlnAux));
		sdb->auxs[i]->n_hit   = 0;
		sdb->auxs[i]->n_idx   = i;
		sdb->auxs[i]->hits    = init_sr_hitv(4);
		sdb->auxs[i]->incs    = init_sr_hitv(4);
		sdb->auxs[i]->kmers   = init_u32list(32);
		sdb->auxs[i]->rids    = init_u32list(64);
		sdb->auxs[i]->in      = NULL;
		sdb->auxs[i]->out     = NULL;
		sdb->auxs[i]->buf_in  = NULL;
		sdb->auxs[i]->buf_out = NULL;
		sdb->auxs[i]->tmp_in  = NULL;
		sdb->auxs[i]->tmp_out = NULL;
		memset(&sdb->auxs[i]->HIT, 0, sizeof(SR_AlnHit));
		sdb->auxs[i]->HIT.rid1 = 0xFFFFFFFFU;
	}

	sdb->idx1 = 0;
	sdb->idx2 = 1;
	sdb->idx_strand = 0;

	sdb->output = sr_output_hit;
	sdb->load   = sr_load_hit;
	sdb->dump   = sr_dump_hit;

	return sdb;
}

void sr_set_align_parameters(SR_SeqDB *sdb, uint8_t strand, uint8_t min_overlap, float min_sm, uint8_t max_mm, int allow_gap){
	sdb->strand_type    = strand;
	sdb->min_overlap    = min_overlap;
	sdb->min_similarity = min_sm;
	sdb->max_mismatch   = max_mm;
	sdb->allow_gap      = allow_gap;
}

void sr_set_filter_parameters(SR_SeqDB *sdb, uint32_t min_complexity, uint32_t max_hits_per_query){
	sdb->min_complexity = min_complexity;
	sdb->max_hits_per_query = max_hits_per_query;
}

void sr_set_output_func(SR_SeqDB *sdb, sr_output_hit_func output){ 
	if(output) sdb->output = output;
}

void sr_set_load_dump_func(SR_SeqDB *sdb, sr_load_hit_func load, sr_dump_hit_func dump){
	sdb->load = load;
	sdb->dump = dump;
}

int sr_push_sdb(SR_SeqDB *sdb, char *seq, uint16_t seqlen, int visuable){
	if(sdb->rd_len && seqlen != sdb->rd_len) return 0;
	encap_u64list(sdb->rd_seqs, (sdb->offset + seqlen) / 32 + 2);
	seq2bits(as_array_u64list(sdb->rd_seqs), sdb->offset, seq, seqlen);
	if(sdb->rd_len == 0){
		push_u16list(sdb->rd_lens, seqlen);
		if(sdb->n_rd % 16 == 0) push_u64list(sdb->rd_offs, sdb->offset);
	}
	sdb->offset += seqlen;
	push_u32list(sdb->links, 0xFFFFFFFFU);
	if(visuable) one2bitvec(sdb->rd_visuable);
	else zero2bitvec(sdb->rd_visuable);
	zero2bitvec(sdb->rd_filtered);
	sdb->n_rd ++;
	return 1;
}

void sr_ready_sdb(SR_SeqDB *sdb){
	return;
	sdb = sdb;
}

void sr_reset_sdb(SR_SeqDB *sdb){
	if(sdb->hits) clear_sr_hitv(sdb->hits);
	if(sdb->out) fflush(sdb->out);
	sdb->n_rd = 0;
	sdb->offset = SR_RDSEQ_PADDING;
	clear_u64list(sdb->rd_offs);
	clear_u16list(sdb->rd_lens);
	clear_u32list(sdb->links);
	clear_bitvec(sdb->rd_visuable);
	clear_bitvec(sdb->rd_filtered);
}

void sr_free_sdb(SR_SeqDB *sdb){
	uint32_t i;
	free_u64list(sdb->rd_seqs);
	free_u64list(sdb->rd_offs);
	free_u16list(sdb->rd_lens);
	for(i=0;i<sdb->n_cpu;i++){
		if(sdb->indexs[i]) free_sr_prmhash(sdb->indexs[i]);
		free_sr_hitv(sdb->auxs[i]->hits);
		free_sr_hitv(sdb->auxs[i]->incs);
		free_u32list(sdb->auxs[i]->kmers);
		free_u32list(sdb->auxs[i]->rids);
		free(sdb->auxs[i]);
	}
	free(sdb->indexs);
	free(sdb->auxs);
	free_u32list(sdb->links);
	free_bitvec(sdb->rd_visuable);
	free_bitvec(sdb->rd_filtered);
	if(sdb->hits) free_sr_hitv(sdb->hits);
	if(sdb->out == stdout){
		free(sdb->prefix);
	} else if(sdb->out) fclose(sdb->out);
	free(sdb);
}

void sr_indexing_sdb(SR_SeqDB *sdb, uint32_t n_idx){
	SR_Primer P, *p;
	sr_prmhash *index;
	uint32_t i, j, rd_len, c;
	uint64_t offset;
	int exists;
	if(sdb->indexs[n_idx]) free_sr_prmhash(sdb->indexs[n_idx]);
	index = sdb->indexs[n_idx] = init_sr_prmhash(1023);
	P.rid = 0;
	for(i=0;i<sdb->n_rd;i++){
		P.kmer = 0;
		if(get_bitvec(sdb->rd_filtered, i)) continue;
		offset = sr_rdseq_offset(sdb, i);
		rd_len = sr_rdseq_length(sdb, i);
		if(rd_len < (sdb->idx2 + 1) * sdb->kmer_size) continue;
		rd_len --;
		if(sdb->idx_strand){
			for(j=0;j<sdb->kmer_size;j++){
				c = (~bits2bit(as_array_u64list(sdb->rd_seqs), offset + rd_len - (sdb->idx1 * sdb->kmer_size + j))) & 0x03;
				P.kmer = (P.kmer << 2) | c;
			}
			for(j=0;j<sdb->kmer_size;j++){
				c = (~bits2bit(as_array_u64list(sdb->rd_seqs), offset + rd_len - (sdb->idx2 * sdb->kmer_size + j))) & 0x03;
				P.kmer = (P.kmer << 2) | c;
			}
		} else {
			for(j=0;j<sdb->kmer_size;j++){
				c = bits2bit(as_array_u64list(sdb->rd_seqs), offset + sdb->idx1 * sdb->kmer_size + j);
				P.kmer = (P.kmer << 2) | c;
			}
			for(j=0;j<sdb->kmer_size;j++){
				c = bits2bit(as_array_u64list(sdb->rd_seqs), offset + sdb->idx2 * sdb->kmer_size + j);
				P.kmer = (P.kmer << 2) | c;
			}
		}
		if((sr_kmer_code(P) % sdb->n_cpu) != n_idx) continue;
		p = prepare_sr_prmhash(index, P, &exists);
		if(exists){
			p->cnt ++;
			set_u32list(sdb->links, i, p->rid);
		} else {
			p->kmer = P.kmer;
			p->cnt  = 1;
			set_u32list(sdb->links, i, i);
		}
		p->rid = i;
	}
}

void sr_aln_load_rd(SR_SeqDB *sdb, uint32_t rid){
	SR_AlnAux *aux;
	uint32_t i;
	uint64_t off;
	SR_AlnHit *h;
	aux = sdb->auxs[rid % sdb->n_cpu];
	off = sr_rdseq_offset(sdb, rid);
	aux->rd_len = sr_rdseq_length(sdb, rid);
	for(i=0;i<aux->rd_len;i+=32){
		aux->rd_seq[i>>5] = sub32seqbits(as_array_u64list(sdb->rd_seqs), off + i);
	}
	if(aux->HIT.rid1 != 0xFFFFFFFFU){
		if(aux->HIT.rid1 < rid) sdb->dump(sdb, &aux->HIT, aux->out);
		else if(aux->HIT.rid1 > rid) return;
		else {
			h = next_ref_sr_hitv(aux->hits);
			*h = aux->HIT;
		}
	}
	while(sdb->load(sdb, &aux->HIT, aux->in)){
		if(aux->HIT.rid1 > rid) return;
		if(aux->HIT.rid1 < rid) sdb->dump(sdb, &aux->HIT, aux->out);
		else {
			h = next_ref_sr_hitv(aux->hits);
			*h = aux->HIT;
		}
	}
	aux->HIT.rid1 = 0xFFFFFFFFU;
}

void sr_aln_dump_rd(SR_SeqDB *sdb, uint32_t rid){
	SR_AlnAux *aux;
	SR_AlnHit *h1, *h2;
	uint32_t i, j;
	int cmp;
	aux = sdb->auxs[rid % sdb->n_cpu];
	i = 0;
	j = 0;
	bsort_sr_hits(ref_sr_hitv(aux->incs, 0), count_sr_hitv(aux->incs), NULL);
	while(1){
		if(i < count_sr_hitv(aux->hits)){
			h1 = ref_sr_hitv(aux->hits, i);
			if(j < count_sr_hitv(aux->incs)){
				h2 = ref_sr_hitv(aux->incs, j);
				cmp = sr_cmp_hits(*h1, *h2, NULL);
				if(cmp == 0){
					fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
					fflush(stdout);
					abort();
				} else if(cmp < 0){
					sdb->dump(sdb, h1, aux->out);
					i ++;
				} else {
					sdb->dump(sdb, h2, aux->out);
					j ++;
				}
			} else {
				sdb->dump(sdb, h1, aux->out);
				i ++;
			}
		} else {
			if(j < count_sr_hitv(aux->incs)){
				h2 = ref_sr_hitv(aux->incs, j);
				sdb->dump(sdb, h2, aux->out);
				j ++;
			} else {
				break;
			}
		}
	}
	clear_sr_hitv(aux->hits);
	clear_sr_hitv(aux->incs);
}

static inline uint64_t sr_stdaln_lh3(AlnCigar* cigars, SR_SeqDB *sdb, uint64_t off1, uint64_t off2, uint32_t dir1, uint32_t dir2, uint32_t len1, uint32_t len2){
	AlnAln *aln;
	char seq1[MAX_RD_LEN+1], seq2[MAX_RD_LEN+1];
	path_t *p;
	int i, n_cigar, mm, mn;
	if(dir1) bits2revseq(seq1, as_array_u64list(sdb->rd_seqs), off1, len1);
	else bits2seq(seq1, as_array_u64list(sdb->rd_seqs), off1, len1);
	if(dir2) bits2revseq(seq2, as_array_u64list(sdb->rd_seqs), off2, len2);
	else bits2seq(seq2, as_array_u64list(sdb->rd_seqs), off2, len2);
	aln = aln_stdaln_aux(seq1, seq2, &sdb->ap, 1, 0, len1, len2);
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

int sr_aln_hit_core_dp(SR_SeqDB *sdb, uint64_t off1, uint64_t off2, uint32_t len1, uint32_t len2, SR_AlnHit HIT){
	SR_AlnHit *hit;
	SR_AlnAux *aux;
	uint64_t ret;
	uint32_t n_cigar, mm, ma;
	ret = sr_stdaln_lh3(sdb->cigars, sdb, off1, off2, HIT.dir1, HIT.dir2, len1, len2);
	n_cigar = ret & 0x3FF;
	mm = (ret >> 10) & 0xFFFFF;
	ma = (ret >> 30) & 0xFFFFF;
	if(n_cigar == 0 || n_cigar > MAX_N_CIGAR) return 2;
	if(sdb->cigars[0].type == ALN_CIGAR_TYPE_DEL) return 2;
	HIT.off = (sdb->cigars[0].type == ALN_CIGAR_TYPE_MAT)? 0 : sdb->cigars[0].len;
	aux = sdb->auxs[HIT.rid1 % sdb->n_cpu];
	if(ma >= sdb->min_overlap && mm <= sdb->max_mismatch && mm <= (uint32_t)(ma * (1 - sdb->min_similarity) + 0.5)){
		hit = next_ref_sr_hitv(aux->incs);
		hit->rid1    = HIT.rid1;
		hit->rid2    = HIT.rid2;
		hit->dir1    = HIT.dir1;
		hit->dir2    = HIT.dir2;
		hit->off     = HIT.off;
		hit->n_ol    = ma;
		hit->n_mm    = mm;
		memcpy(hit->cigars, sdb->cigars, n_cigar * sizeof(AlnCigar));
		hit->n_cigar = n_cigar;
		return 0;
	} else return 2;
}

int sr_aln_hit_core_simp(SR_SeqDB *sdb, uint64_t off1, uint64_t off2, uint32_t len1, uint32_t len2, SR_AlnHit HIT){
	SR_AlnHit *hit;
	SR_AlnAux *aux;
	uint32_t i, mo, mm, max, len;
	uint64_t p1, p2;
	uint8_t c1, c2;
	len = len1 - HIT.off;
	if(len > len2) len = len2;
	if(HIT.dir1){
		off1 = off1 + len1 - (len + HIT.off);
	} else {
		off1 = off1 + HIT.off;
	}
	if(HIT.dir2){
		off2 = off2 + len2 - len;
	} else {
		off2 = off2;
	}
	mm = 0;
	max = len * (1 -  sdb->min_similarity) + 0.5;
	if(max > sdb->max_mismatch) max = sdb->max_mismatch;
	if(sdb->allow_gap){
		mo = 0;
		if(HIT.dir1 ^ HIT.dir2){
			for(i=0;i<len;i++){
				c1 = bits2bit(as_array_u64list(sdb->rd_seqs), off1 + i);
				c2 = (~bits2bit(as_array_u64list(sdb->rd_seqs), off2 + len - 1 - i)) & 0x03;
				if(c1 != c2){
					if(mm == 0) mo = i + 1;
					mm ++;
					if(mm > max) return 2 + mo;
				}
			}
		} else {
			for(i=0;i<len;i++){
				c1 = bits2bit(as_array_u64list(sdb->rd_seqs), off1 + i);
				c2 = bits2bit(as_array_u64list(sdb->rd_seqs), off2 + i);
				if(c1 != c2){
					if(mm == 0) mo = i + 1;
					mm ++;
					if(mm > max) return 2 + mo;
				}
			}
		}
	} else {
		if(HIT.dir1 ^ HIT.dir2){
			for(i=0;i+32<=len;i+=32){
				p1 = sub32seqbits(ref_u64list(sdb->rd_seqs, 0), off1 + i);
				p2 = sub32seqbits(ref_u64list(sdb->rd_seqs, 0), off2 + len - i - 32);
				p2 = dna_rev_seq(p2, 32);
				mm += count_ones_bit64(dna_xor2ones(p1 ^ p2));
				if(mm > max) return 2;
			}
			if(i < len){
				p1 = sub32seqbits(ref_u64list(sdb->rd_seqs, 0), off1 + i);
				p2 = sub32seqbits(ref_u64list(sdb->rd_seqs, 0), off2 + len - i - 32);
				p2 = dna_rev_seq(p2, 32);
				mm += count_ones_bit64(dna_xor2ones((p1 ^ p2) >> ((32 - (len - i)) << 1)));
				if(mm > max) return 2;
			}
		} else {
			for(i=0;i+32<=len;i+=32){
				p1 = sub32seqbits(ref_u64list(sdb->rd_seqs, 0), off1 + i);
				p2 = sub32seqbits(ref_u64list(sdb->rd_seqs, 0), off2 + i);
				mm += count_ones_bit64(dna_xor2ones(p1 ^ p2));
				if(mm > max) return 2;
			}
			if(i < len){
				p1 = sub32seqbits(ref_u64list(sdb->rd_seqs, 0), off1 + i);
				p2 = sub32seqbits(ref_u64list(sdb->rd_seqs, 0), off2 + i);
				mm += count_ones_bit64(dna_xor2ones((p1 ^ p2) >> ((32 - (len - i)) << 1)));
				if(mm > max) return 2;
			}
		}
	}
	aux = sdb->auxs[HIT.rid1 % sdb->n_cpu];
	if(1){
		hit = next_ref_sr_hitv(aux->incs);
		hit->rid1    = HIT.rid1;
		hit->rid2    = HIT.rid2;
		hit->dir1    = HIT.dir1;
		hit->dir2    = HIT.dir2;
		hit->off     = HIT.off;
		hit->n_ol    = len;
		hit->n_mm    = mm;
		hit->n_cigar = append_cigars(hit->cigars, 0, ALN_CIGAR_TYPE_INS, HIT.off);
		hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_MAT, len);
		if(len == len2) hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_INS, len1 - len);
		else hit->n_cigar = append_cigars(hit->cigars, hit->n_cigar, ALN_CIGAR_TYPE_DEL, len2 - len);
		return 0;
	} else return 1;
}

uint32_t sr_aln_rd_core(SR_SeqDB *sdb, uint32_t rid, uint32_t dir){
	SR_Primer P, *p;
	SR_AlnAux *aux;
	SR_AlnHit HIT;
	uint32_t i, j, vis1, vis2;
	uint32_t c, idx, len1, len2, len, mm, klen, ret;
	uint64_t kmer, off1, off2, n_idx;
	n_idx = rid % sdb->n_cpu;
	aux = sdb->auxs[rid % sdb->n_cpu];
	P.rid = 0;
	P.cnt = 0;
	ret = 0;
	clear_u32list(aux->kmers);
	kmer = 0;
	for(i=0;i<aux->rd_len;i++){
		if(dir){
			j = aux->rd_len - i - 1;
			c = (~(aux->rd_seq[j>>5] >> (((~j) & 0x1F) << 1))) & 0x03;
		} else {
			c =   (aux->rd_seq[i>>5] >> (((~i) & 0x1F) << 1))  & 0x03;
		}
		kmer = ((kmer << 2) | c) & sdb->kmer_mask;
		if(i + 1 < sdb->kmer_size) continue;
		push_u32list(aux->kmers, kmer);
	}
	klen = (sdb->idx2 - sdb->idx1) * sdb->kmer_size;
	memset(&HIT, 0, sizeof(SR_AlnHit));
	HIT.dir1 = dir;
	HIT.dir2 = sdb->idx_strand;
	off1 = sr_rdseq_offset(sdb, rid);
	len1 = sr_rdseq_length(sdb, rid);
	vis1 = get_bitvec(sdb->rd_visuable, rid);
	for(i=sdb->idx1*sdb->kmer_size;i+klen<count_u32list(aux->kmers);i++){
		if(i - sdb->idx1 * sdb->kmer_size + sdb->min_overlap > len1) break;
		P.kmer = (((uint64_t)get_u32list(aux->kmers, i)) << (sdb->kmer_size << 1)) | get_u32list(aux->kmers, i + klen);
		p = get_sr_prmhash(sdb->indexs[sr_kmer_code(P)%sdb->n_cpu], P);
		if(p == NULL) continue;
		if(p->cnt > sdb->max_hits_per_query) continue;
		idx = p->rid;
		HIT.rid1 = rid;
		HIT.off  = i - sdb->idx1 * sdb->kmer_size;
		clear_u32list(aux->rids);
		while(1){
			if(idx <= rid) break;
			if(vis1 || (vis2 = get_bitvec(sdb->rd_visuable, idx))){
				HIT.rid2 = idx;
				if(count_sr_hitv(aux->hits) < 8){
					if(locate_sr_hitv(aux->hits, HIT, 0) == count_sr_hitv(aux->hits)){
						if(locate_sr_hitv(aux->incs, HIT, 0) == count_sr_hitv(aux->incs)) push_u32list(aux->rids, idx);
					}
				} else {
					if(query_sr_hits(ref_sr_hitv(aux->hits, 0), count_sr_hitv(aux->hits), HIT, NULL) < 0){
						if(locate_sr_hitv(aux->incs, HIT, 0) == count_sr_hitv(aux->incs)) push_u32list(aux->rids, idx);
					}
				}
			}
			if(get_u32list(sdb->links, idx) == idx) break;
			idx = get_u32list(sdb->links, idx);
		}
		for(j=count_u32list(aux->rids);j>0;j--){
			idx = get_u32list(aux->rids, j - 1);
			off2 = sr_rdseq_offset(sdb, idx);
			len2 = sr_rdseq_length(sdb, idx);
			len = len1 - HIT.off;
			if(len > len2) len = len2;
			if(len < sdb->min_overlap) continue;
			HIT.rid2 = idx;
			HIT.off  = i - sdb->idx1 * sdb->kmer_size;
			mm = sr_aln_hit_core_simp(sdb, off1, off2, len1, len2, HIT);
			if(mm == 0) ret ++;
			else if(mm >= sdb->gap_init + 2 && sdb->allow_gap && sr_aln_hit_core_dp(sdb, off1, off2, len1, len2, HIT) == 0) ret ++;
		}
	}
	return ret;
}

thread_begin_def(msr);
SR_SeqDB *sdb;
int task;
thread_end_def(msr);

void _msr_task_2(SR_SeqDB *sdb, int t_idx){
	SR_AlnAux *aux;
	uint64_t i;
	int fd;
	aux = sdb->auxs[t_idx];
	aux->n_hit = 0;
	aux->HIT.rid1 = 0xFFFFFFFFU;
	if(sdb->prefix){
		if(aux->in){
			fclose(aux->in);
			unlink(aux->tmp_in);
			free(aux->tmp_in);
			aux->tmp_in = NULL;
		}
		if(aux->out){
			aux->in = aux->out;
			aux->tmp_in = aux->tmp_out;
			fseek(aux->in, 0, SEEK_SET);
		} else {
			aux->tmp_in = malloc(strlen(sdb->prefix) + 20);
			sprintf(aux->tmp_in, "%s.%03d.%d%d%02d.XXXXXX", sdb->prefix, 0, 0, 0, t_idx);
			if((fd = mkstemp(aux->tmp_in)) == -1){
				fprintf(stderr, " -- Fail to create temp file as '%s' in %s -- %s:%d --\n", aux->tmp_in, __FUNCTION__, __FILE__, __LINE__);
				fflush(stdout);
				abort();
			}
			close(fd);
			aux->in = fopen(aux->tmp_in, "r");
		}
		{
			aux->tmp_out = malloc(strlen(sdb->prefix) + 20);
			sprintf(aux->tmp_out, "%s.%03d.%d%d%02d.XXXXXX", sdb->prefix, sdb->idx_strand, sdb->idx1, sdb->idx2, t_idx);
			if((fd = mkstemp(aux->tmp_out)) == -1){
				fprintf(stderr, " -- Fail to create temp file as '%s' in %s -- %s:%d --\n", aux->tmp_out, __FUNCTION__, __FILE__, __LINE__);
				fflush(stdout);
				abort();
			}
			close(fd);
			aux->out = fopen(aux->tmp_out, "w+");
		}
	} else {
		if(aux->in){
			fclose(aux->in);
			free(aux->buf_in);
			aux->buf_in = NULL;
		}
		if(aux->out){
			fflush(aux->out);
			aux->in = aux->out;
			aux->buf_in = aux->buf_out;
		} else {
			aux->in = open_memstream(&aux->buf_in, &aux->n_in);
			fflush(aux->in);
		}
		fseek(aux->in, 0, SEEK_SET);
		aux->out = open_memstream(&aux->buf_out, &aux->n_out);
	}
	for(i=t_idx;i<sdb->n_rd;i+=sdb->n_cpu){
		if(get_bitvec(sdb->rd_filtered, i)) continue;
		sr_aln_load_rd(sdb, i);
		if(sdb->strand_type == 1){
			aux->n_hit += sr_aln_rd_core(sdb, i, sdb->idx_strand);
		} else {
			aux->n_hit += sr_aln_rd_core(sdb, i, 0);
			aux->n_hit += sr_aln_rd_core(sdb, i, 1);
		}
		sr_aln_dump_rd(sdb, i);
	}
}

uint32_t cal_seq_complexity(SR_SeqDB *sdb, uint32_t rid){
	uint64_t off;
	uint32_t i, tri, len, cpx;
	off = sr_rdseq_offset(sdb, rid);
	len = sr_rdseq_length(sdb, rid);
	memset(sdb->tris, 0, 64 * sizeof(uint8_t));
	tri = 0;
	tri = (tri << 2) | bits2bit(ref_u64list(sdb->rd_seqs, 0), off + 0);
	tri = (tri << 2) | bits2bit(ref_u64list(sdb->rd_seqs, 0), off + 1);
	for(i=2;i<len;i++){
		tri = (tri << 2) | bits2bit(ref_u64list(sdb->rd_seqs, 0), off + i);
		sdb->tris[tri & 0x3F] ++;
	}
	cpx = 0;
	for(i=0;i<64;i++){
		if(sdb->tris[i]) cpx ++;
	}
	return cpx;
}

uint32_t filter_low_complexity(SR_SeqDB *sdb, uint32_t t_idx){
	uint32_t i, ret;
	ret = 0;
	if(sdb->min_complexity){
		for(i=t_idx;i<sdb->n_rd;i+=sdb->n_cpu){
			if(cal_seq_complexity(sdb, i) < sdb->min_complexity){
				one_bitvec(sdb->rd_filtered, i);
				ret ++;
			}
		}
	}
	return ret;
}

thread_begin_func(msr);
thread_begin_loop(msr);
if(msr->task == 1){
	sr_indexing_sdb(msr->sdb, msr->t_idx);
} else if(msr->task == 2){
	_msr_task_2(msr->sdb, msr->t_idx);
} else if(msr->task == 3){
	filter_low_complexity(msr->sdb, msr->t_idx);
}
thread_end_loop(msr);
thread_end_func(msr);

uint64_t sr_aln_sdb(SR_SeqDB *sdb){
	thread_preprocess(msr);
	uint64_t i, cnt;
	uint32_t rid;
	int n;
	SR_AlnAux *aux;
	SR_AlnHit *h;
	cnt = 0;
	if(sdb->n_cpu > 1){
		thread_begin_init(msr, ((int)sdb->n_cpu));
		msr->sdb = sdb;
		msr->task = 0;
		thread_end_init(msr);
		thread_begin_iter(msr);
		msr->task = 3;
		thread_wake(msr);
		thread_end_iter(msr);
		thread_waitfor_all_idle(msr);
	} else {
		filter_low_complexity(sdb, 0);
	}
	for(sdb->idx1=0;sdb->idx1<ALN_MAX_IDX;sdb->idx1++){
		for(sdb->idx2=sdb->idx1+1;sdb->idx2<ALN_MAX_IDX;sdb->idx2++){
			for(sdb->idx_strand=0;sdb->idx_strand<2;sdb->idx_strand++){
				if(sdb->n_cpu == 1){
					sr_indexing_sdb(sdb, 0);
					_msr_task_2(sdb, 0);
				} else {
					thread_begin_iter(msr);
					msr->task = 1;
					thread_wake(msr);
					thread_end_iter(msr);
					thread_waitfor_all_idle(msr);
					thread_begin_iter(msr);
					msr->task = 2;
					thread_wake(msr);
					thread_end_iter(msr);
					thread_waitfor_all_idle(msr);
				}
				for(i=0;i<sdb->n_cpu;i++){
					free_sr_prmhash(sdb->indexs[i]);
					sdb->indexs[i] = NULL;
					cnt += sdb->auxs[i]->n_hit;
				}
			}
		}
	}
	if(sdb->n_cpu > 1){
		thread_begin_close(msr);
		thread_end_close(msr);
	}
	for(i=0;i<sdb->n_cpu;i++){
		if(sdb->auxs[i]->in){
			if(sdb->prefix){
				fclose(sdb->auxs[i]->in);
				unlink(sdb->auxs[i]->tmp_in);
				free(sdb->auxs[i]->tmp_in);
				sdb->auxs[i]->tmp_in = NULL;
			} else {
				fclose(sdb->auxs[i]->in);
				free(sdb->auxs[i]->buf_in);
				sdb->auxs[i]->buf_in = NULL;
			}
			sdb->auxs[i]->in = NULL;
		}
		if(sdb->auxs[i]->out) fseek(sdb->auxs[i]->out, 0, SEEK_SET);
		sdb->auxs[i]->HIT.rid1 = 0xFFFFFFFFU;
	}
	for(rid=0;rid<sdb->n_rd;rid++){
		aux = sdb->auxs[rid % sdb->n_cpu];
		if(aux->out){
			if(aux->HIT.rid1 != 0xFFFFFFFFU){
				if(aux->HIT.rid1 > rid) continue;
				if(sdb->prefix) sdb->output(sdb, &aux->HIT, sdb->out);
				else {
					h = next_ref_sr_hitv(sdb->hits);
					*h = aux->HIT;
				}
			}
			while((n = sdb->load(sdb, &aux->HIT, aux->out))){
				if(aux->HIT.rid1 > rid) break;
				if(sdb->prefix) sdb->output(sdb, &aux->HIT, sdb->out);
				else {
					h = next_ref_sr_hitv(sdb->hits);
					*h = aux->HIT;
				}
			}
			if(n != 1){
				if(sdb->prefix){
					fclose(aux->out);
					unlink(aux->tmp_out);
					free(aux->tmp_out);
					aux->tmp_out = NULL;
				} else {
					fclose(aux->out);
					free(aux->buf_out);
					aux->buf_out = NULL;
				}
				aux->out = NULL;
			}
		}
	}
	for(i=0;i<sdb->n_cpu;i++){
		aux = sdb->auxs[i];
		if(aux->out){
			if(sdb->prefix){
				fclose(aux->out);
				unlink(aux->tmp_out);
				free(aux->tmp_out);
				aux->tmp_out = NULL;
			} else {
				fclose(aux->out);
				free(aux->buf_out);
				aux->buf_out = NULL;
			}
			aux->out = NULL;
		}
	}
	if(sdb->prefix) fflush(sdb->out);
	return cnt;
}

void sr_output_hit(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out){
	char str[MAX_RD_LEN + 1];
	cigars2string(hit->cigars, hit->n_cigar, str);
	fprintf(out, "%u\t%c\t%u\t%c\t%s", hit->rid1, "+-"[hit->dir1], hit->rid2, "+-"[hit->dir2], str);
	if(1){
		if(hit->dir1) bits2revseq(str, as_array_u64list(sdb->rd_seqs), sr_rdseq_offset(sdb, hit->rid1), sr_rdseq_length(sdb, hit->rid1));
		else bits2seq(str, as_array_u64list(sdb->rd_seqs), sr_rdseq_offset(sdb, hit->rid1), sr_rdseq_length(sdb, hit->rid1));
		fprintf(out, "\t%s", str);
		if(hit->dir2) bits2revseq(str, as_array_u64list(sdb->rd_seqs), sr_rdseq_offset(sdb, hit->rid2), sr_rdseq_length(sdb, hit->rid2));
		else bits2seq(str, as_array_u64list(sdb->rd_seqs), sr_rdseq_offset(sdb, hit->rid2), sr_rdseq_length(sdb, hit->rid2));
		fprintf(out, "\t%s", str);
	}
	fprintf(out, "\n");
}

int sr_load_hit(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *in){
	return (fread(hit, sizeof(SR_AlnHit), 1, in) == 1);
	sdb = sdb;
}

int sr_dump_hit(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out){
	return (fwrite(hit, sizeof(SR_AlnHit), 1, out) == 1);
	sdb = sdb;
}
