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
 
#ifndef __SEQ_GRAPH_RJ_H
#define __SEQ_GRAPH_RJ_H

#include "bitvec.h"
#include "dna.h"
#include "aln_cigar.h"
#include "list.h"
#include "hashset.h"
#include "file_reader.h"
#include "sort.h"

#define MAX_RD_LEN	2000

define_list_core(u32v, uint32_t, uint16_t, 16);
define_list_ext(u32v, uint32_t, uint16_t, native_number_equals);

typedef struct {
	uint32_t rid;
	uint32_t roff:11, rdir:1, rlen:11, n_cigar:9;
	AlnCigar *cigars;
} ReadPath;

define_list_core(plist, ReadPath, uint32_t, 2);

typedef struct {
	uint32_t rids[2];
	uint32_t rlens[2];
	uint32_t rdirs[2];
	uint32_t roffs[2];
	uint32_t rgaps[2];
	uint32_t olens[2];
	uint32_t n_cigar;
	AlnCigar *cigars;
} Overlap;

define_list(olist, Overlap);

typedef struct {
	plist *paths;
	uint32_t len:30, vdir1:1, vdir2:1;
	uint32_t vid1, vid2;
} Edge;

define_list(elist, Edge);

typedef u32v Vertex;

define_list(vlist, Vertex);

typedef struct {
	elist   *edges;
	vlist   *verts;
	u32list *ecache;
	u32list *vcache;
	BitVec  *v_flags, *r2v_dir, *e_flags, *e_chgs;
	u32list *r2v_map;
	uint32_t n_rd;
	uint32_t rd_len;
	uint32_t max_rd_len;
	u64list *rd_seqs;
	u16list *rd_lens;
	u64list *rd_seq_offs;
	olist   *ols;

	AlnCigar *aux_cigars[3];
	BitVec   *flags[3];
	u32list  *goffs[3];
	u32list  *gends[3];
} SeqGraph;

typedef struct {
	SeqGraph *g;
	uint32_t rid, rdir;
	uint32_t rlen, plen;
	uint32_t vid, vdir;
	uint32_t nvid, nvdir;
	uint32_t eid, edir;
	uint32_t pid;
	uint32_t off, len;
} ReadTracker;

#ifdef __CPLUSPLUS
extern "C" {
#endif

#ifdef __CPLUSPLUS
}
#endif

#endif
