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
 
#ifndef __ALL_PATH_RJ_H
#define __ALL_PATH_RJ_H

#include "sr_aln.h"
#include "string.h"
#include "sort.h"
#include "dna.h"
#include "heap.h"

define_native_list(u32v, uint32_t);

#define hit_rid(h, ref) ((h).rid1)
#define cmp_hit_rid(h1, h2, ref) num_cmp_script(h1, h2, ref, hit_rid)
define_search_array(locate_hit, SR_AlnHit, cmp_hit_rid);

#define MAX_CTG_SIZE	5000

typedef struct {
	uint32_t rid;
	uint32_t prev;
	uint32_t used:1,inasm:1, dir:1, off:21, len:8;
	uint32_t idx, prev_idx;
	u32v *links;
} Node;

define_list(nodev, Node);

define_list(cigarv, AlnCigar);
define_list(cigarvv, cigarv);

define_list(strv, String*);

typedef struct Trace {
	uint32_t p_off, h_off, nid1, nid2, idx;
} Trace;

define_list(tracev, Trace);

typedef struct {
	uint32_t acgtn[5];
} Compose;

define_list(compv, Compose);

typedef struct {
	SR_SeqDB *sdb;
	nodev *nodes;
	u32v  *path;
	cigarvv *cigars;
	AlnCigar pat[64], mat[64];
	String *cns;
	compv  *comps;
	strv   *matrix;
	u32list *begs;
	u32list *ends;
	u32list *fwds, *revs;
	tracev  *traces;
	Heap    *heap;
	uint32_t head;
	uint32_t pair_score;
	uint32_t cns_len;
	uint32_t p1, p2;
	uint32_t max_ins, min_ins;
	uint32_t max_n_rd;
	uint32_t sn;
	
	uint32_t exec_simp;

	uint32_t output_aln, output_dot, output_cns;
	char     aux_str[1024];
} MyPath;

#ifdef __CPLUSPLUS
extern "C" {
#endif

MyPath* init_mypath(uint32_t n_cpu, uint8_t kmer_size, uint8_t rd_len, uint8_t min_ol, float min_sm, uint8_t max_mm,  int allow_gap, uint32_t min_ins, uint32_t max_ins, uint32_t max_n_rd);

void free_mypath(MyPath *mp);

void reset_mypath(MyPath *mp);

void push_mypath(MyPath *mp, uint32_t rid, char *seq, uint8_t seqlen, int dir);

void aln_mypath(MyPath *mp);

int solve_mypath(MyPath *mp);

void do_consensus_mypath(MyPath *mp);

void output_consensus_mypath(MyPath *mp, FILE *lr_seqs, FILE *lr_ctgs);

#ifdef __CPLUSPLUS
}
#endif

#endif
