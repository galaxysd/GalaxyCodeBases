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
#include "thread.h"
#include "stdaln.h"
#include <stdio.h>

#define MAX_KMER_SIZE	16
#define MIN_KMER_SIZE	3
#define ALN_MAX_IDX	4
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
	uint32_t dir1:1, dir2:1, off:10, n_ol:10, n_mm:6, n_cigar:4;
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
	uint32_t n_hit;
	uint64_t rd_seq[32];
	uint32_t rd_len;
	sr_hitv  *hits;
	sr_hitv  *incs;
	u32list  *kmers;
	u32list  *rids;
	FILE     *in, *out;
	char     *tmp_in, *tmp_out;
	char     *buf_in, *buf_out;
	size_t   n_in, n_out;
	SR_AlnHit HIT;
} SR_AlnAux;

struct SR_SeqDB;

typedef void (*sr_output_hit_func)(struct SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out);
typedef int  (*sr_dump_hit_func)(struct SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out);
typedef int  (*sr_load_hit_func)(struct SR_SeqDB *sdb, SR_AlnHit *hit, FILE *in);

typedef struct SR_SeqDB {
	u64list  *rd_seqs;
	u16list  *rd_lens;
	u64list  *rd_offs;
	BitVec   *rd_visuable;
	BitVec   *rd_filtered;
	uint32_t rd_len;
	uint32_t n_rd;

	char     *prefix;

	uint32_t strand_type;
	uint32_t min_overlap;
	float    min_similarity;
	uint32_t max_mismatch;
	uint32_t max_hits_per_query;
	int      allow_gap;
	uint32_t min_complexity;
	uint32_t gap_init;
	AlnParam ap;

	uint32_t n_cpu;
	uint32_t idx1, idx2, idx_strand;
	uint32_t kmer_size;
	uint32_t kmer_mask;
	sr_prmhash  **indexs;
	u32list     *links;
	SR_AlnAux   **auxs;

	sr_hitv  *hits;
	FILE     *out;

	uint64_t offset;

	AlnCigar cigars[2 * MAX_RD_LEN];
	uint8_t tris[64];

	sr_output_hit_func output;
	sr_load_hit_func load;
	sr_dump_hit_func dump;
} SR_SeqDB;

#define SR_RDSEQ_PADDING	64

#ifdef __CPLUSPLUS
extern "C" {
#endif

	SR_SeqDB* sr_init_sdb(char *out_prefix, uint32_t n_cpu, uint8_t kmer_size, uint16_t rd_len);

	void sr_set_align_parameters(SR_SeqDB* sdb, uint8_t strand, uint8_t min_overlap, float min_similarity, uint8_t max_mismatch, int allow_gap);

	void sr_set_filter_parameters(SR_SeqDB* sdb, uint32_t min_complexity, uint32_t max_hits_per_query);

	uint64_t sr_rdseq_offset(SR_SeqDB *sdb, uint32_t rid);

	uint16_t sr_rdseq_length(SR_SeqDB *sdb, uint32_t rid);

	void sr_set_output_func(SR_SeqDB *sdb, sr_output_hit_func fun);
	void sr_set_load_dump_func(SR_SeqDB *sdb, sr_load_hit_func load, sr_dump_hit_func dump);
// Push reads #MARK1
	int sr_push_sdb(SR_SeqDB *sdb, char *seq, uint16_t seqlen, int visuable);
// Ready for aligning
	void sr_ready_sdb(SR_SeqDB *sdb);
// Align
	uint64_t sr_aln_sdb(SR_SeqDB *sdb);
// Output hits
	void sr_output_hit(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out);
	int  sr_load_hit(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *in);
	int  sr_dump_hit(SR_SeqDB *sdb, SR_AlnHit *hit, FILE *out);
// goto MARK1
	void sr_reset_sdb(SR_SeqDB *sdb);
// Free
	void sr_free_sdb(SR_SeqDB *sdb);

#ifdef __CPLUSPLUS
}
#endif

#endif
