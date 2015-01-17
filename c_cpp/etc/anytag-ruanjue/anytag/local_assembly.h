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
#include "anytag_aux.h"

typedef struct {
	uint32_t nid;
	uint32_t n_ol:10, off:10, mut_idx:11, closed:1;
	uint32_t hit_idx;
} Edge;

define_list(edgev, Edge);

typedef struct {
	uint32_t rid;
	uint32_t seqdir:1, len:15, closed:1, anchor_type:2;
	int min_ins:16, max_ins:16;
	uint32_t snid:16, snoff:15, visit:1;
	uint32_t bt, bt_idx, bt_off, bt_dir;
	uint32_t last_n_edges[2];
	edgev   *edges[2];
} Node;

define_list(nodev, Node);

typedef struct {
	uint32_t nid;
	uint32_t off:16, len:15, dir:1;
	uint32_t hit_idx;
} Layout;

define_list(layv, Layout);

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
	uint32_t bases[4];
} comb_base_t;

define_list(cbv, comb_base_t);

define_list(cigarv, AlnCigar);

define_list(cigarm, cigarv*);

define_list(strm, b8list*);

typedef struct {
	nodev     *nodes;
	uint32_t  m_node;
	SR_SeqDB  *sdb;
	uint32_t  min_ol, ol_step, cur_ol[2];
	SR_AlnAux *aux, *long_aux;
	Tracker   *trackers[2];
	uint32_t  p1, p2;
	layv      *lays;
	String    *skeleton;
	SR_SeqDB  *sp_sdb;
	sr_hitv   *sp_hits;
	u32list   *ords[2];
	cbv       *cbases;
	cigarm    *cigars[2];
	strm      *msa[2];
	u32list   *bcovs, *pcovs;
	String    *cns_seq;
	uint32_t  cns_id, cns_len;
	uint32_t  n_sr, n_l, n_r, n_layout;
	uint32_t  min_bcov, min_pcov;
	uint32_t  check_ins, max_ins, min_ins;
	u64list   *bits;
	BitVec    *flags;
	char     chs[3][MAX_RD_LEN];
	ATOptions *opt;
} LGraph;

#ifdef __CPLUSPLUS
extern "C" {
#endif

	LGraph* init_lgraph(ATOptions *opt);
	//LGraph* init_lgraph(uint8_t n_seed, uint8_t kmer_size, uint8_t rd_len, uint8_t min_ol, float min_sm, uint8_t max_mm, uint32_t min_ins, uint32_t max_ins);

	void trunon_inserts_checking_lgraph(LGraph *g);
	void trunoff_inserts_checking_lgraph(LGraph *g);

	void free_lgraph(LGraph *g);

	void reset_lgraph(LGraph *g);

	// STEP 1: collecting sequences
	// anchor_type: 0, left; 1: right
	// min_loc, max_loc: where should this read locate
	void push_lgraph(LGraph *g, uint32_t rid, char *seq, uint8_t seqlen, uint32_t seqdir, uint32_t anchor_type, int min_ins, int max_ins);

	// STEP 2
	void ready_lgraph(LGraph *g);

	// STEP 3: searcing connective path in overlap graph
	int layout_lgraph(LGraph *g);

	uint32_t simplify_lgraph(LGraph *g);

	void print_dot_lgraph(LGraph *g, FILE *out);

	// STEP 4: generating skeleton of pesudo-Sanger
	void skeleton_lgraph(LGraph *g);
	void get_skeleton_lgraph(LGraph *g, String *skeleton);
	void set_skeleton_lgraph(LGraph *g, String *skeleton);

	// STEP 5: adding paired alignments to skeleton
	void push_aln_lgraph(LGraph *g, SR_AlnHit *sp1, SR_AlnHit *sp2);

	uint32_t align_skeleton_lgraph(LGraph *g, SR_SeqDB *sdb);

	// STEP 6: consensus
	void consensus_lgraph(LGraph *g);

	// STEP 7: output
	void output_lgraph(LGraph *g, FILE *seqf, FILE *ctgf);

#ifdef __CPLUSPLUS
}
#endif

#endif
