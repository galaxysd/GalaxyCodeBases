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
 
#ifndef __SIMPLE_ASM_RJ_H
#define __SIMPLE_ASM_RJ_H

#include "string.h"
#include "sr_aln.h"
#include "list.h"
#include "dna.h"
#include "stdaln.h"

typedef struct {
	uint32_t seqid, len;
	uint32_t ctg_id, ctg_dir:1, ctg_off:20, used:1, rank:10;
} SimpSeqInfo;

define_list(seqv, SimpSeqInfo);

typedef struct {
	uint32_t len:31, closed:1;
	String   *seq;
	u32list  *sids;
} SimpContigInfo;

define_list(ctgv, SimpContigInfo);

typedef struct {
	SR_SeqDB *sra;
	u32list  *rids;
	seqv     *rds;
	ctgv     *ctgs;
	u64hash  *r2r;
	sr_hitv  *ols;
	uint32_t iter_idx;
	AlnParam ap;
	int      aln_ext_size;
} SimpAssembler;

static inline SimpAssembler* init_simpasm(uint32_t n_cpu, uint32_t kmer_size, uint32_t rd_len, uint32_t strand, uint32_t min_ol, float min_sm, uint32_t max_mm, int allow_gap){
	SimpAssembler *sa;
	sa = malloc(sizeof(SimpAssembler));
	sa->sra  = sr_init_sdb(NULL, n_cpu, kmer_size, rd_len);
	sr_set_align_parameters(sa->sra, strand, min_ol, min_sm, max_mm, allow_gap);
	sa->rds  = init_seqv(16);
	sa->ctgs = init_ctgv(16);
	sa->r2r  = init_u64hash(13);
	sa->ols  = init_sr_hitv(64);
	sa->ap   = (AlnParam){10, 2, 2, aln_sm_nt, 16, 75};
	sa->aln_ext_size = 10;
	return sa;
}

static inline void free_simpasm(SimpAssembler *sa){
	SimpContigInfo *ctg;
	uint32_t i;
	sr_free_sdb(sa->sra);
	free_seqv(sa->rds);
	for(i=0;i<count_ctgv(sa->ctgs);i++){
		ctg = ref_ctgv(sa->ctgs, i);
		free_string(ctg->seq);
		free_u32list(ctg->sids);
	}
	free_ctgv(sa->ctgs);
	free_u64hash(sa->r2r);
	free_sr_hitv(sa->ols);
	free(sa);
}

static inline void reset_simpasm(SimpAssembler *sa){
	SimpContigInfo *ctg;
	uint32_t i;
	clear_seqv(sa->rds);
	for(i=0;i<count_ctgv(sa->ctgs);i++){
		ctg = ref_ctgv(sa->ctgs, i);
		free_string(ctg->seq);
		free_u32list(ctg->sids);
	}
	clear_ctgv(sa->ctgs);
	sr_reset_sdb(sa->sra);
}

static inline void push_simpasm(SimpAssembler *sa, uint32_t seqid, char *seq, uint32_t seqlen, uint8_t rank){
	SimpSeqInfo *rd;
	SimpContigInfo *ctg;
	rd = next_ref_seqv(sa->rds);
	rd->seqid   = seqid;
	rd->len     = seqlen;
	rd->ctg_id  = count_ctgv(sa->ctgs);
	rd->ctg_dir = 0;
	rd->ctg_off = 0;
	rd->used    = 0;
	rd->rank    = rank;
	ctg = next_ref_ctgv(sa->ctgs);
	ctg->seq    = init_string(seqlen);
	ctg->sids   = init_u32list(2);
	ctg->len    = seqlen;
	ctg->closed = 0;
	append_string(ctg->seq, seq, seqlen);
	uc_string(ctg->seq);
	push_u32list(ctg->sids, count_seqv(sa->rds) - 1);
	sr_push_sdb(sa->sra, seq, seqlen, 1); // visuable = 1, this seq will be visuable in alignment
}

static inline int cmp_sr_alnhit(const void *e1, const void *e2){
	SR_AlnHit *h1, *h2;
	h1 = (SR_AlnHit*)e1;
	h2 = (SR_AlnHit*)e2;
	if(h1->n_ol == h2->n_ol){
		if(h1->n_mm == h2->n_mm) return 0;
		else if(h1->n_mm < h2->n_mm) return 1;
		else return -1;
	} else if(h1->n_ol < h2->n_ol) return 1;
	else return -1;
}

static inline void simple_reverse_contig(SimpAssembler *sa, uint32_t ctg_id){
	SimpContigInfo *ctg;
	SimpSeqInfo *rd;
	uint32_t i;
	ctg = ref_ctgv(sa->ctgs, ctg_id);
	if(ctg->closed) return;
	reverse_dna(ctg->seq->string, ctg->len);
	for(i=0;i<count_u32list(ctg->sids);i++){
		rd = ref_seqv(sa->rds, get_u32list(ctg->sids, i));
		rd->ctg_dir = ! rd->ctg_dir;
		rd->ctg_off = ctg->len - (rd->ctg_off + rd->len);
	}
}

static inline void simple_move_rids(SimpAssembler *sa, uint32_t dst, uint32_t src, int off){
	SimpContigInfo *ctg1, *ctg2;
	SimpSeqInfo *rd;
	int i;
	ctg1 = ref_ctgv(sa->ctgs, dst);
	ctg2 = ref_ctgv(sa->ctgs, src);
	for(i=0;i<(int)count_u32list(ctg2->sids);i++){
		rd = ref_seqv(sa->rds, get_u32list(ctg2->sids, i));
		rd->ctg_id = dst;
		rd->ctg_off += off;
		push_u32list(ctg1->sids, get_u32list(ctg2->sids, i));
	}
	if(off < 0){
		for(i=0;i<(int)count_u32list(ctg1->sids);i++){
			rd = ref_seqv(sa->rds, get_u32list(ctg1->sids, i));
			rd->ctg_off -= off;
		}
	}
}

static inline int simple_join_contigs(SimpAssembler *sa, SR_AlnHit *hit){
	SimpContigInfo *ctg1, *ctg2;
	SimpSeqInfo *rd1, *rd2;
	AlnAln *aa;
	int i, dir1, dir2, off1, off2, off3, off4, mm, mn, ret;
	ret = 0;
	rd1 = ref_seqv(sa->rds, hit->rid1);
	rd2 = ref_seqv(sa->rds, hit->rid2);
	dir1 = rd1->ctg_dir ^ hit->dir1;
	dir2 = rd1->ctg_dir ^ hit->dir2;
	if(dir1 == dir2) dir1 = dir2 = 0;
	if(dir1) simple_reverse_contig(sa, rd1->ctg_id);
	if(dir2) simple_reverse_contig(sa, rd2->ctg_id);
	ctg1 = ref_ctgv(sa->ctgs, rd1->ctg_id);
	ctg2 = ref_ctgv(sa->ctgs, rd2->ctg_id);
	dir1 = rd1->ctg_dir;
	dir2 = rd2->ctg_dir;
	off1 = rd1->ctg_off - sa->aln_ext_size + (int)hit->off;
	off2 = rd2->ctg_off - sa->aln_ext_size;
	if(off1 < 0) off1 = 0;
	if(off2 < 0) off2 = 0;
	if(off1 < off2){ off2 -= off1; off1 = 0; }
	else {off1 -= off2; off2 = 0; }
	aa = aln_stdaln_aux(ctg1->seq->string + off1, ctg2->seq->string + off2, &sa->ap, 1, 1, ctg1->len - off1, ctg2->len - off2);
	mn = mm = 0;
	off3 = off4 = -1;
	for(i=aa->path_len-1;i>=0;i--){
		if(aa->path[i].ctype == FROM_M){
			if(off3 == -1) off3 = aa->path_len - 1 - i;
			if(off4 == -1) off4 = aa->path_len - 1 - i;
			mn ++;
			if(ctg1->seq->string[aa->path[i].i-1+off1] != ctg2->seq->string[aa->path[i].j-1+off2]) mm ++;
		} else if(aa->path[i].ctype == FROM_D){
			if(off3 == -1) off3 = aa->path_len - 1 - i;
		} else {
			if(off4 == -1) off4 = aa->path_len - 1 - i;
		}
	}
	//fprintf(stdout, "\n\nHIT %d %d %d, %d %d %d + %d\n", hit->rid1, rd1->ctg_dir, rd1->ctg_off, hit->rid2, rd2->ctg_dir, rd2->ctg_off, hit->off);
	//sr_output_hit(sa->sra, hit, stdout, 1);
	//fprintf(stdout, "CTG %s\nCTG %s\n", ctg1->seq->string, ctg2->seq->string);
	//fprintf(stdout, "%s\n%s\n", ctg1->seq->string + off1, ctg2->seq->string + off2);
	//fprintf(stdout, "%d %d, %d %d\n", off1, off3, off2, off4);
	//fprintf(stdout, "%d/%d\n", mm, mn);
	//fprintf(stdout, "%s\n%s\n%s\n\n", aa->out1, aa->outm, aa->out2);
	//fflush(stdout);
	if(mn >= (int)sa->sra->min_overlap && mm <= mn * (1 - sa->sra->min_similarity)){
		off3 = off3 - off1;
		off4 = off4 - off2;
		simple_move_rids(sa, rd1->ctg_id, rd2->ctg_id, off4 - off3);
		if(off1){
			trunc_string(ctg1->seq, off1);
		} else {
			clear_string(ctg1->seq);
			if(off2) append_string(ctg1->seq, ctg2->seq->string, off2);
		}
		for(i=0;i<aa->path_len;i++){
			add_char_string(ctg1->seq, (aa->out1[i] == '-')? aa->out2[i] : aa->out1[i]);
		}
		ctg1->len = ctg1->seq->size;
		ctg2->closed = 1;
		//fprintf(stdout, "## %s\n", ctg1->seq->string);
		//fflush(stdout);
		ret = 1;
	} else {
		//fprintf(stdout, "EE\n");
		//fflush(stdout);
	}
	aln_free_AlnAln(aa);
	return ret;
}

static inline void simple_assemble(SimpAssembler *sa){
	SimpSeqInfo *rd1, *rd2;
	SR_AlnHit *hit;
	uint64_t key, *k;
	uint32_t i, rank_type, rtype;
	int exists;
	sr_ready_sdb(sa->sra);
	sr_aln_sdb(sa->sra);
	//for(i=0;i<count_sr_hitv(sa->sra->hits);i++){
		//hit = ref_sr_hitv(sa->sra->hits, i);
		//sr_output_hit(sa->sra, hit, stdout, 1);
	//}
	qsort(as_array_sr_hitv(sa->sra->hits), count_sr_hitv(sa->sra->hits), sizeof(SR_AlnHit), cmp_sr_alnhit);
	clear_u64hash(sa->r2r);
	clear_sr_hitv(sa->ols);
	for(i=0;i<count_sr_hitv(sa->sra->hits);i++){
		hit = ref_sr_hitv(sa->sra->hits, i);
		key = (hit->rid1 < hit->rid2)? ((((uint64_t)hit->rid1) << 32) | hit->rid2) : ((((uint64_t)hit->rid2) << 32) | hit->rid1);
		k = prepare_u64hash(sa->r2r, key, &exists);
		if(exists) continue;
		*k = key;
		push_sr_hitv(sa->ols, *hit);
	}
	for(rank_type=0;rank_type<3;rank_type++){
		for(i=0;i<count_sr_hitv(sa->ols);i++){
			hit = ref_sr_hitv(sa->ols, i);
			rd1 = ref_seqv(sa->rds, hit->rid1);
			rd2 = ref_seqv(sa->rds, hit->rid2);
			if(rd1->ctg_id == rd2->ctg_id) continue;
			if(rd1->rank == rd2->rank) rtype = 0;
			else if(rd1->rank + 1 == rd2->rank || rd1->rank == rd2->rank + 1) rtype = 1;
			else rtype = 2;
			if(rtype != rank_type) continue;
			simple_join_contigs(sa, hit);
		}
	}
}

static inline void begin_iter_simpasm(SimpAssembler *sa){ sa->iter_idx = 0; }

static inline SimpContigInfo* iter_simpasm(SimpAssembler *sa){
	SimpContigInfo *ctg;
	while(sa->iter_idx < count_ctgv(sa->ctgs)){
		ctg = ref_ctgv(sa->ctgs, sa->iter_idx ++);
		if(ctg->closed == 0) return ctg;
	}
	return NULL;
}

#endif
