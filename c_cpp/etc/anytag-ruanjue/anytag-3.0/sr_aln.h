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
 
#ifndef __SHORT_READS_ALIGN_RJ_H
#define __SHORT_READS_ALIGN_RJ_H
#include "list.h"
#include "sort.h"
#include "bitvec.h"
#include "dna.h"
#include "aln_cigar.h"
#include "hashset.h"
#include "static_hash.h"
#include "stdaln.h"
#include <stdio.h>

#define MAX_KMER_SIZE	16
#define MIN_KMER_SIZE	3
#define MAX_N_SEED	8
#define MAX_N_IDX	112
#define MAX_N_CIGAR	6
#define MAX_RD_LEN	4000
#define MAX_KMER_CNT	0xFFFFF

typedef struct {
	uint32_t kmer1, kmer2, off, cnt;
} SR_Primer;

#define sr_kmer_code(p) ((p).kmer1 + (((uint64_t)(p).kmer2) << 32))
#define sr_primer_hashcode(p) u64hashcode((((uint64_t)(p).kmer1) << 32) | (p).kmer2)
#define sr_primer_equals(p1, p2) (((p1).kmer1 == (p2).kmer1) && ((p1).kmer2 == (p2).kmer2))
define_hashset(sr_prmhash, SR_Primer, sr_primer_hashcode, sr_primer_equals);
define_static_hashset(sr_prm_hash, sr_prmhash, SR_Primer, sr_primer_hashcode, sr_primer_equals);

typedef struct {
	uint32_t rid;
	uint32_t dir1:1, dir2:1, n_ol:8, n_mm:4, is_gap:1, n_cigar:3;
	int32_t  off:14;
	AlnCigar cigars[MAX_N_CIGAR];
} SR_AlnHit;

define_list_core(sr_hitv, SR_AlnHit, uint32_t, 64);

typedef struct {
	uint32_t idx;
	uint32_t fix;
	uint32_t dir1, dir2;
	int      off1, off2;
	uint32_t ptr, num;
} SR_AlnPtr;

define_list(sr_ptrv, SR_AlnPtr);

typedef struct {
	uint32_t rid;
	uint64_t rd_seq[MAX_RD_LEN / 32 + 1], rv_seq[MAX_RD_LEN / 32 + 1];
	uint32_t rd_len;
	u32list  *fkmers, *rkmers;
	u32list  *hits_per_off;
	sr_ptrv  *ptrs;
	sr_hitv  *hits;
	u64list  *pids;
	AlnCigar cigars[2 * MAX_RD_LEN];
	BitVec   *rd_flags;
	u32list  *rd_tested;
	int      auto_clear;

	uint32_t strand_type; // 1: forward, 2: reverse, 3: both
	uint32_t max_ol;
	uint32_t min_ol;
	float    min_sm;
	uint32_t max_hits, max_hits_per_off;
	int      allow_gap;
	uint32_t hit_id_range[2];
	uint32_t best_ols[2];
	uint32_t best_cnt[2];
	int      best_mode; // 0: all hits; 1: only best hits; 2: only one best

	int      rec_cigar;
} SR_AlnAux;

typedef struct SR_SeqDB {
	BaseBank *rdseqs;
	BitVec   *rd_filtered;
	uint32_t rd_len;
	uint32_t n_rd;
	uint32_t min_complexity;
	float    hash_factor;

	AlnParam ap;

	uint32_t n_seed, n_idx;
	uint32_t idx_pats[MAX_N_IDX][3];
	uint32_t kmer_size;
	uint32_t kmer_mask;
	sr_prm_hash  **indexs;
	u32list     **links;
} SR_SeqDB;

#define SR_RDSEQ_PADDING	64

#ifdef __CPLUSPLUS
extern "C" {
#endif

	SR_AlnAux* sr_init_aux();

#define sr_set_aux_strand(aux, strand) (aux)->strand_type = (strand)
#define sr_set_aux_overlap(aux, min, max) { (aux)->min_ol = (min); (aux)->max_ol = (max); }
#define sr_set_aux_min_similarity(aux, sm) (aux)->min_sm = (sm)
#define sr_set_aux_gap(aux, gap) (aux)->allow_gap = (gap)
#define sr_set_aux_hit_id_region(aux, min, max) { (aux)->hit_id_range[0] = (min); (aux)->hit_id_range[1] = (max); }
#define sr_set_aux_best_mode(aux, mode) (aux)->best_mode = (mode)
#define sr_set_aux_best_ols(aux, left, right) { (aux)->best_ols[0] = (left); (aux)->best_ols[1] = (right); }
#define sr_set_aux_max_hits(aux, max) (aux)->max_hits = (max)
#define sr_set_aux_max_hits_per_off(aux, max) (aux)->max_hits_per_off = (max)
#define sr_set_aux_auto_clear(aux, is_clear) (aux)->auto_clear = (is_clear)
#define sr_set_aux_cigar(aux, need_cigar) (aux)->rec_cigar = (need_cigar)

	void sr_fit_aux2sdb(SR_AlnAux *aux, SR_SeqDB *sdb);

	void sr_clear_aux(SR_AlnAux *aux);

	void sr_clear_rd_flags_aux(SR_AlnAux *aux);

	void sr_block_aux(SR_AlnAux *aux, uint32_t rid);
	// Will block this rid forever
	void sr_block_forever_aux(SR_AlnAux *aux, uint32_t rid);

	void sr_free_aux(SR_AlnAux *aux);

	SR_SeqDB* sr_init_sdb(uint32_t kmer_size, uint32_t seed_kmer_num, uint32_t rd_len, uint32_t min_complexity);
	SR_SeqDB* sr_init2_sdb(uint32_t kmer_size, uint32_t seed_kmer_num, uint32_t rd_len, uint32_t min_complexity, float hash_factor);

	void sr_mask_seed_sdb(SR_SeqDB *sdb, uint32_t seed_idx);

#define sr_rdseq_length(sdb, rid) ((sdb)->rd_len)
#define sr_rdseq_offset(sdb, rid) ((sdb)->rd_len * ((uint64_t)(rid)) + SR_RDSEQ_PADDING)

	size_t sr_dump_sdb(SR_SeqDB *sdb, FILE *out);
	SR_SeqDB* sr_load_sdb(FILE *inp);

	void sr_reset_sdb(SR_SeqDB *sdb);
// Push reads
	int sr_push_sdb(SR_SeqDB *sdb, char *seq, uint8_t seqlen);
// Ready for aligning
	void sr_ready_sdb(SR_SeqDB *sdb);

	uint32_t sr_mask_low_complexity(SR_SeqDB *sdb, uint32_t n_part, uint32_t part_i);
// Index
#define sr_total_idx(sdb) ((sdb)->n_idx)
	void sr_index_sdb(SR_SeqDB *sdb, uint32_t idx);
	void sr_clear_index_sdb(SR_SeqDB *sdb, uint32_t idx);
// Align
	// Load query
	void sr_aux_load(SR_AlnAux *aux, SR_SeqDB *sdb, uint32_t seqid, uint64_t *seqs, uint64_t seqoff, uint32_t seqlen);
	void sr_aux_load2(SR_AlnAux *aux, SR_SeqDB *sdb, uint32_t seqid, char *seq, uint32_t seqlen);
	uint32_t sr_align_sdb(SR_SeqDB *sdb, SR_AlnAux *aux);

// Free
	void sr_free_sdb(SR_SeqDB *sdb);

#ifdef __CPLUSPLUS
}
#endif

#endif
