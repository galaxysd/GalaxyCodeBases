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
 
#ifndef __PSEUDO_SANGER_RJ_H
#define __PSEUDO_SANGER_RJ_H

#include "dna.h"
#include "bitvec.h"
#include "list.h"
#include "hashset.h"
#include "sort.h"
#include "aln_cigar.h"
#include "heap.h"
#include "string.h"
#include "stdaln.h"
#include "sr_aln.h"

#define MAX_KSIZE	16
#define MAX_READS	65536
#define MAX_TRACE	65535
//#define MAX_RD_LEN	256

#ifdef ENABLE_LOOP
#define MAX_LOOP	32
#endif

typedef struct {
	uint32_t nid:20, off:12;
} rdoff_t;
define_list(rdoffv, rdoff_t);

typedef struct {
	uint32_t kmer, off, cnt;
} rdkmer_t;
#define rdkmer_hashcode(K) u32hashcode((K).kmer)
#define rdkmer_equals(K1, K2) ((K1).kmer == (K2).kmer)
define_hashset(rdkhash, rdkmer_t, rdkmer_hashcode, rdkmer_equals);

typedef struct {
	uint32_t nid:16, off:8, mm:5, closed:1, flag:1, loop:1;
} edge_t;
define_list_core(edgev, edge_t, uint16_t, 32);

typedef struct {
	uint32_t rid;
	uint32_t seqdir:1, closed:1, bt_nid:20;
	uint32_t bt_eid:12, bt_off:18, visit:1, sort:1, n_edge:10;
	edgev *edges;
} node_t;
define_list(nodev, node_t);

typedef struct {
	uint32_t nid:16, bt:16;
	uint32_t bt_idx:12, uniq:1;
	int      off:19;
#ifdef ENABLE_LOOP
	uint32_t loop;
#endif
} trace_t;
define_list(tracev, trace_t);

typedef struct {
	uint32_t bases[4];
} comb_base_t;
define_list(cbv, comb_base_t);

typedef struct {
	uint32_t rid, idx:17, dir:1, fix:1;
	int      off:13;
} qc_aln_t;
define_list(qc_alnv, qc_aln_t);

typedef struct {
	uint32_t len, qc;
	tracev   *path;
	cbv      *cbs;
	String   *seq;

	BaseBank *qc_seqs;
	qc_alnv  *alns;
} cns_t;
define_list(cnsv, cns_t);

typedef struct {
	uint32_t ps_id, ps_len;
	uint64_t rd_len;
	uint32_t min_ins, max_ins;
	uint32_t n_sr, n_l, n_r;
	BaseBank *rdseqs;
	nodev    *nodes;
	uint32_t m_node;

	uint32_t ksize;
	uint32_t kmask;
	uint32_t min_ol;
	float    min_sm;
	uint32_t out_mode;
	rdkhash  *hash;
	rdoffv   *links;

	BitVec   *conns;
	u8list   *tested;
	tracev   *tree;
	u32list  *stack;
	Heap     *heap;
	tracev   *fast, *slow;
	String   *fast_seq, *slow_seq;
	uint32_t fast_break, slow_break;
	u32hash  *aux_hash;
	uuhash   *loop_hash;
	cnsv     *cnss;
	uint32_t m_cns;

	SR_AlnAux *aux;

	char gdb_seqs[200];
} Graph;

#ifdef __CPLUSPLUS
extern "C" {
#endif

	Graph* init_graph(uint32_t ksize, uint32_t min_ol, float min_sm);

	void free_graph(Graph *g);

	void reset_graph(Graph *g);

	void set_insert_graph(Graph *g, uint32_t min, uint32_t max);

	// Anchor_type: bit1,bit2,bit3, bit1: left; bit2: right; bit3: lay out of ps
	void push_graph(Graph *g, uint32_t rid, char *rdseq, uint8_t rdlen, uint32_t rddir, uint32_t anchor_type);

	void ready_graph(Graph *g);

	void index_graph(Graph *g);

	void align_graph(Graph *g);

	void simplify_graph(Graph *g);
	uint32_t shave_graph(Graph *g);

	void print_dot_graph(Graph *g, FILE *dotf);

	uint32_t allpaths_graph(Graph *g);
	uint32_t validate_paths_graph(Graph *g);
	void consensus_graph(Graph *g);

	void realign_graph(Graph *g, SR_SeqDB *sdb, u32list *lib_cnts, u32list *lib_ins, int do_cns);

	uint32_t output_contigs_graph(Graph *g, FILE *ctgf, FILE *msaf);

#ifdef __CPLUSPLUS
}
#endif

#endif
