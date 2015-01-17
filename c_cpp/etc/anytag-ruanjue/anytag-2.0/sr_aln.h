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
#include "stdaln.h"
#include <stdio.h>

#define MAX_KMER_SIZE	16
#define MIN_KMER_SIZE	3
#define MAX_N_SEED	8
#define MAX_N_IDX	112
#define MAX_N_CIGAR	6
#define MAX_RD_LEN	1000
#define MAX_KMER_CNT	0xFFFFF

typedef struct {
	uint64_t kmer;
	uint64_t rid:44, cnt:20;
} SR_Primer;

#define sr_kmer_code(p) ((p).kmer + ((p).kmer << 32))
#define sr_primer_hashcode(p) u64hashcode((p).kmer)
#define sr_primer_equals(p1, p2) ((p1).kmer == (p2).kmer)
define_hashset(sr_prmhash, SR_Primer, sr_primer_hashcode, sr_primer_equals);

typedef struct {
	uint32_t rid1, rid2;
	uint32_t dir1:1, dir2:1, off:14, n_ol:8, n_mm:4, is_gap:1, n_cigar:3;
	AlnCigar cigars[MAX_N_CIGAR];
} SR_AlnHit;

#define sr_bigger_hits(h1, h2, ref) ((h1).rid1 > (h2).rid1)
#define sr_hits_num_val(h, ref) ((((uint64_t)(h).rid2) << 1) | ((h).dir1 ^ (h).dir2))
#define sr_cmp_hits(h1, h2, ref) num_cmp_script(h1, h2, ref, sr_hits_num_val)
#define sr_hit_equals(h1, h2) (sr_hits_num_val(h1, NULL) == sr_hits_num_val(h2, NULL))
define_list_core(sr_hitv, SR_AlnHit, uint32_t, 64);
define_list_ext( sr_hitv, SR_AlnHit, uint32_t, sr_hit_equals);
define_search_array(query_sr_hits, SR_AlnHit, sr_cmp_hits);
define_bubble_sort(bsort_sr_hits, SR_AlnHit, sr_cmp_hits);

typedef struct {
	uint32_t n_idx;
	uint32_t fix;
	uint32_t dir1, dir2;
	uint32_t off1, off2;
	uint32_t ptr, cnt, num;
} SR_AlnPtr;

define_list(sr_ptrv, SR_AlnPtr);
#define sr_cmp_ptrs(h1, h2, ref) (((h1).off1 == (h2).off1)? 0 : (((h1).off1 < (h2).off1)? 1 : -1))
define_bubble_sort(bsort_sr_ptrs, SR_AlnPtr, sr_cmp_ptrs);

typedef struct {
	uint32_t rid;
	uint64_t rd_seq[MAX_RD_LEN / 32 + 1], rv_seq[MAX_RD_LEN / 32 + 1];
	uint32_t rd_len;
	u32list  *fkmers, *rkmers;
	sr_ptrv  *ptrs;
	sr_hitv  *hits;
	u64list  *pids;
	AlnCigar cigars[2 * MAX_RD_LEN];
} SR_AlnAux;

typedef struct SR_SeqDB {
	u64list  *rd_seqs;
	u16list  *rd_lens;
	u64list  *rd_offs;
	BitVec   *rd_filtered;
	uint32_t rd_len, max_rd_len;
	uint32_t n_rd;

	uint32_t strand_type; // 1: forward, 2: reverse, 3: both
	uint32_t min_overlap;
	float    min_similarity;
	uint32_t max_mismatch;
	uint32_t max_hits_per_query;
	int      allow_gap;
	int      filter_dup;
	uint32_t min_complexity;
	AlnParam ap;

	uint32_t n_seed, n_idx;
	uint32_t idx_pats[MAX_N_IDX][3];
	uint32_t kmer_size;
	uint32_t kmer_mask;
	sr_prmhash  **indexs;
	uint32_t is_fast;
	u64list     **fast_indexs;
	u32list     **links;

	uint64_t offset;

} SR_SeqDB;

#define SR_RDSEQ_PADDING	64

#ifdef __CPLUSPLUS
extern "C" {
#endif

	SR_AlnAux* sr_init_aux();

	void sr_free_aux(SR_AlnAux *aux);

	SR_SeqDB* sr_init_sdb(uint8_t kmer_size, uint8_t seed_kmer_num, uint16_t rd_len);

	void sr_set_align_parameters(SR_SeqDB* sdb, uint8_t strand, uint8_t min_overlap, float min_similarity, uint8_t max_mismatch, int allow_gap);

	void sr_set_filter_parameters(SR_SeqDB* sdb, uint32_t min_complexity, uint32_t max_hits_per_query, int local_filter_dup);

	uint64_t sr_rdseq_offset(SR_SeqDB *sdb, uint32_t rid);

	uint16_t sr_rdseq_length(SR_SeqDB *sdb, uint32_t rid);

// Push reads #MARK1
	int sr_push_sdb(SR_SeqDB *sdb, char *seq, uint16_t seqlen);
// Ready for aligning
	void sr_ready_sdb(SR_SeqDB *sdb);
// Index
	void sr_index_sdb(SR_SeqDB *sdb, uint32_t n_idx);
// Align
	uint32_t sr_align_sdb(SR_SeqDB *sdb, uint32_t rid, SR_AlnAux *aux);
	// Align long reads
	uint32_t sr_align_long_sdb(SR_SeqDB *sdb, SR_AlnAux *aux, uint32_t lrid, uint64_t *seqs, uint64_t seqoff, uint32_t seqlen);

	void sr_output_hit(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out);

	void sr_reset_sdb(SR_SeqDB *sdb);
// Free
	void sr_free_sdb(SR_SeqDB *sdb);

#ifdef __CPLUSPLUS
}
#endif

#endif
