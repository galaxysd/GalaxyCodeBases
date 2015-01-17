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
 
#ifndef __DE_BRUIJN_GRAPH_RJ_H
#define __DE_BRUIJN_GRAPH_RJ_H

#include "hashset.h"
#include "list.h"
#include "dna.h"

typedef struct {
	uint64_t kmer:42, fwd:4, rev:4, var:16;
} Kmer;

#define kmer_hashcode(k) u64hashcode((k).kmer)
#define kmer_equals(k1, k2) ((k1).kmer == (k2).kmer)
define_hashset(khash, Kmer, kmer_hashcode, kmer_equals);

typedef struct {
	uint32_t tid;
	int16_t score;
	uint16_t type;
	uint32_t xlen;
} Link;

define_list(linkv, Link);

typedef struct {
	Kmer     *k;
	uint64_t kmer;
	uint32_t dir:1, closed:1, x:30;
	linkv *links;
} Node;

define_list(nlist, Node);

typedef struct {
	uint64_t kmer;
	uint32_t dir:1, closed:1, x:30;
	uint32_t nid;
} Nod;

#define nod_hashcode(t) u64hashcode((t).kmer | (((uint64_t)(t).x) << 42))
#define nod_equals(t1, t2) ((t1).dir == (t2).dir && (t1).x == (t2).x && (t1).kmer == (t2).kmer)
define_hashset(nhash, Nod, nod_hashcode, nod_equals);

typedef struct {
	uint32_t nid;
	uint8_t type, y_base;
	uint32_t prev;
	int16_t score;
} Trace;

define_list(tlist, Trace);

typedef struct {
	khash *g;
	uint32_t kmer_size;
	uint64_t kmer_mask;

	int max_score;
	int s_mis, s_mat;
	int s_ins_beg, s_ins_ext, s_ins_end;
	int s_del_beg, s_del_ext, s_del_end;

	uint32_t rid, rdir, rlen;
	u8list   *seqs;
	uint32_t roff;
	uint64_t kmer;

	nlist *nodes;
	uint32_t n_node;
	nhash *hash;
	tlist *traces;
	u32list *ends[2];
	u32list *stack;
	u32list *path;
} ReadTracker;

#endif
