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
 
#ifndef __LOCAL_ASSEMBLY_RJ_H
#define __LOCAL_ASSEMBLY_RJ_H

#include "string.h"
#include "aln_cigar.h"
#include "list.h"
#include "sr_aln.h"
#include "dna.h"
#include "heap.h"
#include "sort.h"

typedef struct {
	uint32_t nid;
	uint32_t n_ol, off;
	uint32_t hit_idx;
} Edge;

define_list(edgev, Edge);

typedef struct {
	uint32_t rid;
	uint32_t seqdir:1, len:16, anchor_type:2;
	int min_ins:16, max_ins:16;
	uint32_t bt, bt_idx, bt_off:16, bt_dir:1, visit:1;
	uint32_t in_asm:1, asm_mm:8;
	int      asm_off:18;
	edgev   *edges[2];
} Node;

define_list(nodev, Node);

typedef struct {
	uint32_t nid, n_ol, off;
	uint32_t src_nid, hit_idx;
} Trace;

define_list(tracev, Trace*);

typedef struct {
	uint32_t trace_dir;
	tracev *traces;
	uint32_t m_trace;
	Heap   *heap;
} Tracker;

typedef struct {
	uint32_t nid;
	uint32_t off;
	uint32_t hit_idx;
} Layout;

define_list(layv, Layout);

define_list(cigarv, AlnCigar);
define_list(cigarm, cigarv*);
define_list(strv, String*);

typedef struct {
	uint32_t kmer:20, pos:12;
} Kmer;

#define kmer_hashcode(k) (k).kmer
#define kmer_equals(k1, k2) ((k1).kmer == (k2).kmer)
define_hashset(khash, Kmer, kmer_hashcode, kmer_equals);

typedef struct {
	nodev     *nodes;
	uint32_t  m_node;
	SR_SeqDB  *sdb;
	SR_AlnAux *aux;
	Tracker   *trackers[2];
	uint32_t  p1, p2;
	layv      *lays;
	cigarm    *cms;
	strv      *srs;
	String    *cns_seq;
	u32list   *ords, *aux_ords;
	khash     *cns_index;
	uint32_t  cns_id, cns_len;
	uint32_t  max_ins, min_ins;
	uint32_t  gap_cutoff;
} LGraph;

#ifdef __CPLUSPLUS
extern "C" {
#endif

	LGraph* init_lgraph(uint8_t kmer_size, uint8_t rd_len, uint8_t min_ol, float min_sm, uint8_t max_mm, uint32_t gap_cutoff, uint32_t min_ins, uint32_t max_ins);

	void free_lgraph(LGraph *g);

	void reset_lgraph(LGraph *g);

	// STEP 1: collecting sequences
	// anchor_type: 0, left; 1: right
	// min_loc, max_loc: where should this read locate
	void push_lgraph(LGraph *g, uint32_t rid, char *seq, uint8_t seqlen, uint32_t seqdir, uint32_t anchor_type, int min_ins, int max_ins);

	// STEP 2: aligning reads vs reads
	void align_lgraph(LGraph *g);

	// STEP 3: searcing connective path in overlap graph
	int layout_lgraph(LGraph *g);

	// STEP 4: consensus
	void consensus_lgraph(LGraph *g);

	// STEP 5: realign, do consensus again (optional)
	void realign_lgraph(LGraph *g);

	void output_lgraph(LGraph *g, FILE *cnsf, FILE *samf, FILE *msaf);

#ifdef __CPLUSPLUS
}
#endif

#endif
