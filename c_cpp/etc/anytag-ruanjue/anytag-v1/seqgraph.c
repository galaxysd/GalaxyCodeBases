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
 
#include "seqgraph.h"

#define DEBUG

// -- Macro --

#define v2d_edge_dir(e, vid, vdir) (((e)->vid1 == vid && (e)->vdir1 == vdir)? 0 : (((e)->vid2 == vid && (e)->vdir2 == vdir)? 1 : 2))

#define full_cmp_readpath(p1, p2, ref) (((p1).rid == (p2).rid)?	\
		(((p1).rdir == (p2).rdir)? (((p1).roff == (p1).roff)? 0 : (((p1).roff < (p2).roff)? -1 : 1)) : ((p1).rdir? 1 : -1))	\
		: (((p1).rid < (p2).rid)? -1 : 1))

#define readpath_rid(p, ref) (p).rid
#define cmp_readpath(p1, p2, ref) num_cmp_script(p1, p2, ref, readpath_rid)
#define is_smaller_readpath(p1, p2, ref) ((p1).rid < (p2).rid)
#define is_equal_readpath(p1, p2, ref) ((p1).rid == (p2).rid)

define_search_array(locate_rid_paths, ReadPath, cmp_readpath);

static inline uint32_t read_length(SeqGraph *g, uint32_t rid){ return g->rd_len? g->rd_len : get_u16list(g->rd_lens, rid); }

static inline uint64_t read_offset(SeqGraph *g, uint32_t rid){
	uint64_t off;
	uint16_t *p;
	if(g->rd_len) return ((uint64_t)g->rd_len) * rid;
	off = get_u64list(g->rd_seq_offs, rid / 16);
	p = as_array_u16list(g->rd_lens) + (rid / 16) * 16;
	rid = rid % 16;
	while(rid){
		off += *p;
		rid --;
		p ++;
	}
	return off;
}

// -- ReadPath --

uint32_t rid2vid(SeqGraph *g, uint32_t rid, uint32_t *dir){
	uint32_t vid;
	vid  = (*dir)? (rid * 2 + 1) : (rid * 2);
	*dir = get_bitvec(g->r2v_dir, vid);
	vid  = get_u32list(g->r2v_map, vid);
	return vid;
}

void init_readtracker(ReadTracker *rt, SeqGraph *g, uint32_t rid, uint32_t dir, uint32_t rlen){
	rt->g      = g;
	rt->rid    = rid;
	rt->rdir   = dir;
	rt->rlen   = rlen;
	rt->plen   = 0;
	rt->nvdir  = dir;
	rt->nvid   = rid2vid(g, rid, &rt->nvdir);
}

void set_readtracker(ReadTracker *rt, SeqGraph *g, uint32_t rid, uint32_t rdir, uint32_t rlen, uint32_t plen, uint32_t vid, uint32_t vdir){
	rt->g      = g;
	rt->rid    = rid;
	rt->rdir   = rdir;
	rt->rlen   = rlen;
	rt->plen   = plen;
	rt->nvdir  = vdir;
	rt->nvid   = vid;
}

static inline int cmp_readpath_func(const void *e1, const void *e2){
	return cmp_readpath((*(ReadPath*)e1), (*(ReadPath*)e2), NULL);
}

int next_readtracker(ReadTracker *rt){
	Vertex *v;
	Edge *e;
	ReadPath P, *p;
	uint32_t i;
	uint32_t edir, loop, f, ret;
	int32_t j, pidx;
	if(rt->plen >= rt->rlen) return 0;
	rt->vid  = rt->nvid;
	rt->vdir = rt->nvdir;
	v = ref_vlist(rt->g->verts, rt->vid);
	ret = 0;
	loop = 1;
	P.rid  = rt->rid;
	for(i=0;i<count_u32v(v);i++){
		e = ref_elist(rt->g->edges, get_u32v(v, i));
		f = 0;
		if(e->vid1 == e->vid2){
			if(e->vdir1 == rt->vdir){
				edir = 0;
				if(e->vdir2 == rt->vdir) f = 1;
			} else if(e->vdir2 == rt->vdir){
				edir = 1;
			} else continue;
		} else if(e->vid1 == rt->vid){
			if(e->vdir1 != rt->vdir) continue;
			edir = 0;
		} else {
			if(e->vdir2 != rt->vdir) continue;
			edir = 1;
		}
		if(count_plist(e->paths) < 12){
			pidx = 0;
			while(pidx < (int)count_plist(e->paths) && ref_plist(e->paths, pidx)->rid < rt->rid) pidx ++;
		} else {
			pidx = locate_rid_paths(as_array_plist(e->paths), count_plist(e->paths), P, NULL);
			if(pidx < 0) continue;
			while(pidx && ref_plist(e->paths, pidx - 1)->rid == rt->rid) pidx --;
		}
		for(j=pidx;(unsigned)j<count_plist(e->paths);j++){
			p = ref_plist(e->paths, j);
			if(p->rid != rt->rid) break;
			if(f){
				edir = (p->rdir != rt->rdir);
			} else if((edir ^ p->rdir) != rt->rdir) continue;
			if(rt->rdir){
				if((p->roff + p->rlen) != (int)(rt->rlen - rt->plen)) continue;
			} else {
				if(p->roff != rt->plen) continue;
			}
			if(loop || p->rlen){
				rt->eid  = get_u32v(v, i);
				rt->edir = edir;
				rt->pid  = j;
				rt->off  = rt->plen;
				rt->len  = p->rlen;
				rt->plen = rt->plen + p->rlen;
				if(edir){
					rt->nvid  = e->vid1;
					rt->nvdir = !e->vdir1;
				} else {
					rt->nvid  = e->vid2;
					rt->nvdir = !e->vdir2;
				}
				ret = 1;
				loop = (e->vid1 == e->vid2);
			}
			if(p->rlen) return 1;
		}
	}
	return ret;
}

int find_readpath_edge(SeqGraph *g, uint32_t eid, uint32_t edir, uint32_t rid, uint32_t rdir, uint32_t rlen, uint32_t plen, uint32_t *len, uint32_t *p_idx){
	Edge *e;
	ReadPath P, *p;
	int pidx;
	e = ref_elist(g->edges, eid);
	if(count_plist(e->paths) < 12){
		pidx = 0;
		while(pidx < (int)count_plist(e->paths) && ref_plist(e->paths, pidx)->rid < rid) pidx ++;
	} else {
		P.rid = rid;
		pidx = locate_rid_paths(as_array_plist(e->paths), count_plist(e->paths), P, NULL);
		if(pidx < 0) return 0;
		while(pidx && ref_plist(e->paths, pidx - 1)->rid == rid) pidx --;
	}
	for(;pidx<(int)count_plist(e->paths);pidx++){
		p = ref_plist(e->paths, pidx);
		if(p->rid != rid) break;
		if((edir ^ p->rdir) != rdir) continue;
		if(rdir){
			if((p->roff + p->rlen) != (int)(rlen - plen)) continue;
		} else {
			if(p->roff != plen) continue;
		}
		if(len) *len = p->rlen;
		if(p_idx) *p_idx = pidx;
		return 1;
	}
	return 0;
}

void align_read_path(SeqGraph *g, Overlap *ol, uint32_t idx, ReadTracker *rt){
	Edge *e;
	ReadPath *p;
	uint32_t i;
	int m1, m2, m3;
	uint32_t x, y;
	if(ol->rgaps[idx] == 0) return;
	x = select_cigars_seqlen(ol->cigars, ol->n_cigar, rt->off, idx);
	y = select_cigars_seqlen(ol->cigars, ol->n_cigar, rt->off + rt->len, idx);
	if(y - x == rt->len) return;
	m1 = sub_cigars(g->aux_cigars[0], ol->cigars, ol->n_cigar, x, y - x);
	e = ref_elist(g->edges, rt->eid);
	p = ref_plist(e->paths, rt->pid);
	if(rt->edir) reverse_cigars(g->aux_cigars[0], m1);
	m2 = compile_cigars(g->aux_cigars[1], p->cigars, p->n_cigar, g->aux_cigars[0], m1, 0);
	if(m2 == 1 && g->aux_cigars[1][0].type == ALN_CIGAR_TYPE_MAT) return;
	for(i=0;i<count_plist(e->paths);i++){
		p = ref_plist(e->paths, i);
		m3 = apply_cigars(g->aux_cigars[2], p->cigars, p->n_cigar, g->aux_cigars[1], m2);
		if(m3 > p->n_cigar){
			p->cigars = realloc(p->cigars, sizeof(AlnCigar) * m3);
		}
		memcpy(p->cigars, g->aux_cigars[2], sizeof(AlnCigar) * m3);
		p->n_cigar = m3;
	}
	e->len += y - x - rt->len;
}

static inline void move_readpath(ReadPath *p1, ReadPath *p2, int dir){
	*p1 = *p2;
	p2->cigars  = NULL;
	if(dir){
		p1->rdir = !p2->rdir;
		if(p2->n_cigar == 0) return;
		reverse_cigars(p1->cigars, p1->n_cigar);
	}
}

void sort_edge_paths(SeqGraph *g, uint32_t eid){
	Edge *e;
	e = ref_elist(g->edges, eid);
	qsort(as_array_plist(e->paths), count_plist(e->paths), sizeof(ReadPath), cmp_readpath_func);
}

void change_r2v_map(SeqGraph *g, uint32_t rid, uint32_t rdir, uint32_t vid, uint32_t vdir){
	rid = rid * 2 + rdir;
	if(vdir) one_bitvec(g->r2v_dir, rid);
	else zero_bitvec(g->r2v_dir, rid);
	set_u32list(g->r2v_map, rid, vid);
}

// -- Edge + Vertex --

#define edge_changed(g, eid)

uint32_t lend_edge(SeqGraph *g){
	Edge *e;
	uint32_t eid;
	if(pop_u32list(g->ecache, &eid)){
		e = ref_elist(g->edges, eid);
		one_bitvec(g->e_flags, eid);
		edge_changed(g, eid);
	} else {
		eid = count_elist(g->edges);
		e = next_ref_elist(g->edges);
		one2bitvec(g->e_flags);
		one2bitvec(g->e_chgs);
	}
	e->paths = init_plist(2);
	return eid;
}

void return_edge(SeqGraph *g,uint32_t eid){
	Edge *e;
	ReadPath *p;
	uint32_t i;
	zero_bitvec(g->e_flags, eid);
	e = ref_elist(g->edges, eid);
	if(e->paths){
		for(i=0;i<count_plist(e->paths);i++){
			p = ref_plist(e->paths, i);
			if(p->cigars) free(p->cigars);
		}
		free_plist(e->paths);
	}
	e->paths = NULL;
	push_u32list(g->ecache, eid);
}

uint32_t lend_vertex(SeqGraph *g){
	Vertex *v;
	uint32_t vid;
	if(pop_u32list(g->vcache, &vid)){
		clear_u32v(ref_vlist(g->verts, vid));
		one_bitvec(g->v_flags, vid);
	} else {
		vid = count_vlist(g->verts);
		v = next_ref_vlist(g->verts);
		u32v_init(v, 2);
		one2bitvec(g->v_flags);
	}
	return vid;
}

void return_vertex(SeqGraph *g, uint32_t vid){
	zero_bitvec(g->v_flags, vid);
	push_u32list(g->vcache, vid);
}

int edge_v2v(SeqGraph *g, uint32_t vid1, uint32_t vdir1, uint32_t vid2, uint32_t vdir2, uint32_t *eid, uint32_t *edir){
	Vertex *v1;
	Edge *e;
	uint32_t i;
	v1 = ref_vlist(g->verts, vid1);
	for(i=0;i<count_u32v(v1);i++){
		e = ref_elist(g->edges, get_u32v(v1, i));
		if(e->vid1 == vid2 && e->vid2 == vid1 && e->vdir1 == vdir2 && e->vdir2 == vdir1){
			*eid  = get_u32v(v1, i);
			*edir = 1;
			return 1;
		}
		if(e->vid1 == vid1 && e->vid2 == vid2 && e->vdir1 == vdir1 && e->vdir2 == vdir2){
			*eid  = get_u32v(v1, i);
			*edir = 0;
			return 1;
		}
	}
	return 0;
}

int has_this_link(SeqGraph *g, uint32_t vid1, uint32_t vid2){
	Vertex *v;
	Edge *e;
	uint32_t i;
	if(get_bitvec(g->v_flags, vid1) == 0) return 0;
	v = ref_vlist(g->verts, vid1);
	for(i=0;i<count_u32v(v);i++){
		e = ref_elist(g->edges, get_u32v(v, i));
		if(e->vid1 == vid2 || e->vid2 == vid2) return 1;
	}
	return 0;
}

void reverse_edge(SeqGraph *g, uint32_t eid){
	Edge *e;
	ReadPath *p;
	uint32_t i, t;
	edge_changed(g, eid);
	e = ref_elist(g->edges, eid);
	t = e->vid1; e->vid1 = e->vid2; e->vid2 = t;
	t = e->vdir1; e->vdir1 = e->vdir2; e->vdir2 = t;
	for(i=0;i<count_plist(e->paths);i++){
		p = ref_plist(e->paths, i);
		p->rdir = !p->rdir;
		reverse_cigars(p->cigars, p->n_cigar);
	}
}

void remove_edge(SeqGraph *g, uint32_t eid){
	Vertex *v;
	Edge *e;
	e = ref_elist(g->edges, eid);
	v = ref_vlist(g->verts, e->vid1);
	delete_u32v(v, eid);
	v = ref_vlist(g->verts, e->vid2);
	delete_u32v(v, eid);
	return_edge(g, eid);
}

uint32_t find_edges_v2d(SeqGraph *g, uint32_t vid, uint32_t vdir, u32v *eids){
	Vertex *v;
	Edge *e;
	uint32_t i, eid, n;
	n = 0;
	v = ref_vlist(g->verts, vid);
	for(i=0;i<count_u32v(v);i++){
		eid = get_u32v(v, i);
		e = ref_elist(g->edges, eid);
		if(v2d_edge_dir(e, vid, vdir) != 2){
			push_u32v(eids, eid);
			n ++;
		}
	}
	return n;
}

int merge_equal_edges(SeqGraph *g, uint32_t eid1, uint32_t eid2, uint32_t dir){
	Vertex *v1, *v2;
	Edge *e1, *e2;
	ReadPath *p1, *p2;
	plist *paths;
	uint32_t i1, i2;
	if(eid1 == eid2) return 0;
	edge_changed(g, eid1);
	edge_changed(g, eid2);
	e1 = ref_elist(g->edges, eid1);
	e2 = ref_elist(g->edges, eid2);
	if(e1->len != e2->len) return 0;
	v1 = ref_vlist(g->verts, e1->vid1);
	v2 = ref_vlist(g->verts, e1->vid2);
#ifdef DEBUG
	fprintf(stdout, " -- MERGE %d[%d%c,%d%c] %d[%d%c,%d%c] '%c' --\n", eid1, e1->vid1, "+-"[e1->vdir1], e1->vid2, "+-"[e1->vdir2],
		eid2, e2->vid1, "+-"[e2->vdir1], e2->vid2, "+-"[e2->vdir2], "+-"[dir]);
	fflush(stdout);
#endif
	if(dir){
		if(e1->vid1 != e2->vid2 || e1->vdir1 != e2->vdir2 || e1->vid2 != e2->vid1 || e1->vdir2 != e2->vdir1){
			fprintf(stderr, " -- MERGE %d[%d%c,%d%c] %d[%d%c,%d%c] '%c' --\n", eid1, e1->vid1, "+-"[e1->vdir1], e1->vid2, "+-"[e1->vdir2],
				eid2, e2->vid1, "+-"[e2->vdir1], e2->vid2, "+-"[e2->vdir2], "+-"[dir]);
			fflush(stderr);
			//abort();
			return 0;
		}
	} else {
		if(e1->vid1 != e2->vid1 || e1->vdir1 != e2->vdir1 || e1->vid2 != e2->vid2 || e1->vdir2 != e2->vdir2){
			fprintf(stderr, " -- MERGE %d[%d%c,%d%c] %d[%d%c,%d%c] '%c' --\n", eid1, e1->vid1, "+-"[e1->vdir1], e1->vid2, "+-"[e1->vdir2],
				eid2, e2->vid1, "+-"[e2->vdir1], e2->vid2, "+-"[e2->vdir2], "+-"[dir]);
			fflush(stderr);
			//abort();
			return 0;
		}
	}
	delete_u32v(v1, eid2);
	delete_u32v(v2, eid2);
	paths = init_plist(count_plist(e1->paths) + count_plist(e2->paths));
	i1 = i2 = 0;
	while(i1 < count_plist(e1->paths) && i2 < count_plist(e2->paths)){
		p1 = ref_next_plist(paths);
		if(cmp_readpath(*ref_plist(e1->paths, i1), *ref_plist(e2->paths, i2), NULL) <= 0){
		//if(ref_plist(e1->paths, i1)->rid <= ref_plist(e2->paths, i2)->rid){
			p2 = ref_plist(e1->paths, i1++);
			move_readpath(p1, p2, 0);
		} else {
			p2 = ref_plist(e2->paths, i2++);
			move_readpath(p1, p2, dir);
		}
	}
	while(i1 < count_plist(e1->paths)){
		p1 = ref_next_plist(paths);
		p2 = ref_plist(e1->paths, i1);
		move_readpath(p1, p2, 0);
		i1 ++;
	}
	while(i2 < count_plist(e2->paths)){
		p1 = ref_next_plist(paths);
		p2 = ref_plist(e2->paths, i2);
		move_readpath(p1, p2, dir);
		i2 ++;
	}
	free_plist(e1->paths);
	e1->paths = paths;
	return_edge(g, eid2);
	return 1;
}

void move_edge_v2v(SeqGraph *g, uint32_t eid, uint32_t vid1, uint32_t vid2, uint32_t dir){
	//Vertex *v1, *v2;
	Edge *e;
	//uint32_t vid3, vdir2, vdir3, eid2, edir1, edir2;
	edge_changed(g, eid);
	e  = ref_elist(g->edges, eid);
#ifdef DEBUG
	fprintf(stdout, " -- move E%u(%u:%c -- %u;%c) V%u -> V%u in %s -- %s:%d --\n", eid, e->vid1, "+-"[e->vdir1], e->vid2, "+-"[e->vdir2], vid1, vid2, __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
#endif
	//v1 = ref_vlist(g->verts, vid1);
	//v2 = ref_vlist(g->verts, vid2);
	if(e->vid1 != vid2 && e->vid2 != vid2) push_u32v(ref_vlist(g->verts, vid2), eid);
	if(e->vid1 == e->vid2){
		e->vid1  = vid2;
		e->vid2  = vid2;
		e->vdir1 = e->vdir1 ^ dir;
		e->vdir2 = e->vdir2 ^ dir;
		//vid3  = vid2;
		//vdir2 = e->vdir1;
		//vdir3 = e->vdir2;
		//edir1 = 0;
	} else if(e->vid1 == vid1){
		e->vid1  = vid2;
		e->vdir1 = e->vdir1 ^ dir;
		//vid3  = e->vid2;
		//vdir2 = e->vdir1;
		//vdir3 = e->vdir2;
		//edir1 = 0;
	} else {
		e->vid2  = vid2;
		e->vdir2 = e->vdir2 ^ dir;
		//vid3  = e->vid1;
		//vdir2 = e->vdir2;
		//vdir3 = e->vdir1;
		//edir1 = 1;
	}
	//if(edge_v2v(g, vid2, vdir2, vid3, vdir3, &eid2, &edir2)) merge_equal_edges(g, eid, eid2, edir1 ^ edir2);
}

void map_r2v_edge(SeqGraph *g, uint32_t eid){
	Edge *e;
	ReadPath *p;
	uint32_t i, rlen;
	if(get_bitvec(g->e_flags, eid) == 0) return;
	edge_changed(g, eid);
	e = ref_elist(g->edges, eid);
	for(i=0;i<count_plist(e->paths);i++){
		p = ref_plist(e->paths, i);
		rlen = read_length(g, p->rid);
		if(p->roff == 0){
			if(p->rdir) change_r2v_map(g, p->rid, 0, e->vid2, e->vdir2);
			else change_r2v_map(g, p->rid, 0, e->vid1, e->vdir1);
		}
		if((uint32_t)(p->roff + p->rlen) == rlen){
			if(p->rdir) change_r2v_map(g, p->rid, 1, e->vid1, e->vdir1);
			else change_r2v_map(g, p->rid, 1, e->vid2, e->vdir2);
		}
	}
}

uint32_t split_edge(SeqGraph *g, uint32_t eid1, uint32_t pivot, uint32_t *new_eid){
	Edge *e1, *e2;
	ReadPath *p1, *p2;
	uint32_t vid, i, eid2, nc;
	int clen;
	eid2 = lend_edge(g);
	vid  = lend_vertex(g);
	edge_changed(g, eid1);
	edge_changed(g, eid2);
	e1 = ref_elist(g->edges, eid1);
	e2 = ref_elist(g->edges, eid2);
	push_u32v(ref_vlist(g->verts, vid), eid1);
	push_u32v(ref_vlist(g->verts, vid), eid2);
	if(e1->vid1 != e1->vid2){
		replace_u32v(ref_vlist(g->verts, e1->vid2), eid1, eid2);
	} else push_u32v(ref_vlist(g->verts, e1->vid1), eid2);
	e2->vid2  = e1->vid2;
	e2->vdir2 = e1->vdir2;
	e1->vid2  = vid;
	e1->vdir2 = 1;
	e2->vid1  = vid;
	e2->vdir1 = 0;
	e2->len   = e1->len - pivot;
	e1->len   = pivot;
	for(i=0;i<count_plist(e1->paths);i++){
		p1 = ref_plist(e1->paths, i);
		nc = sub_cigars(g->aux_cigars[0], p1->cigars, p1->n_cigar, pivot, -1);
		cigars_lengths(g->aux_cigars[0], nc, NULL, &clen, NULL);
		if(clen){
			p2 = next_ref_plist(e2->paths);
			p2->rid  = p1->rid;
			p2->rdir = p1->rdir;
			p2->rlen = clen;
			if(p1->rdir){
				p2->roff = p1->roff;
				p1->roff = p1->roff + clen;
			} else {
				p2->roff = p1->roff + p1->rlen - clen;
			}
			p2->cigars = malloc(sizeof(AlnCigar) * nc);
			p2->n_cigar = cat_cigars(p2->cigars, 0, g->aux_cigars[0], nc);
		}
		p1->n_cigar = sub_cigars(p1->cigars, p1->cigars, p1->n_cigar, 0, pivot);
		p1->rlen = p1->rlen - clen;
	}
	//map_r2v_edge(g, eid1);
	//map_r2v_edge(g, eid2);
	*new_eid = eid2;
	return vid;
}

int transfer_v2v(SeqGraph *g, uint32_t vid1, uint32_t vid2, uint32_t dir){
	Vertex *v1;
	uint32_t i, eid;
	if(vid1 == vid2) return 0;
	v1 = ref_vlist(g->verts, vid1);
	for(i=0;i<count_u32v(v1);i++){
		eid = get_u32v(v1, i);
		move_edge_v2v(g, eid, vid1, vid2, dir);
		map_r2v_edge(g, eid);
	}
	return_vertex(g, vid1);
#ifdef DEBUG
	fprintf(stdout, " -- VID_ALIAS V%u -> V%u in %s -- %s:%d --\n", vid1, vid2, __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
#endif
	return 1;
}

uint32_t insert_v2e(SeqGraph *g, uint32_t eid, uint32_t edir, uint32_t vid, uint32_t vdir, uint32_t pivot){
	Edge *e;
	uint32_t tmp_vid, neid;
	if(edir){
		e = ref_elist(g->edges, eid);
		pivot = e->len - pivot;
	}
	tmp_vid = split_edge(g, eid, pivot, &neid);
	transfer_v2v(g, tmp_vid, vid, edir ^ vdir);
	return neid;
}

void add_edge_v2v(SeqGraph *g, uint32_t eid, uint32_t vid1, uint32_t vdir1, uint32_t vid2, uint32_t vdir2){
	Edge *e;
	Vertex *v;
	e = ref_elist(g->edges, eid);
	e->vid1  = vid1;
	e->vdir1 = vdir1;
	e->vid2  = vid2;
	e->vdir2 = vdir2;
	v = ref_vlist(g->verts, vid1);
	push_u32v(v, eid);
	v = ref_vlist(g->verts, vid2);
	push_u32v(v, eid);
}

int connect_edges(SeqGraph *g, uint32_t *eids, BitVec **flags, uint32_t n_eid, uint32_t *new_eid){
	Edge *e1, *e2, *e;
	ReadPath *p1, *p2, *p;
	uint32_t dst, i, j, k;
	uint32_t nc, rid, rdir, roff, rlen, plen, len, pidx;
	if(n_eid < 2) return 0;
	dst = lend_edge(g);
	e   = ref_elist(g->edges, dst);
	e->len = 0;
	for(i=0;i<n_eid;i++) e->len += ref_elist(g->edges, eids[i])->len;
	e1  = ref_elist(g->edges, eids[0]);
	for(i=0;i<count_plist(e1->paths);i++){
		one_bitvec(flags[0], i);
		p1 = ref_plist(e1->paths, i);
		p  = next_ref_plist(e->paths);
		nc = cat_cigars(g->aux_cigars[0], 0, p1->cigars, p1->n_cigar);
		rid  = p1->rid;
		rlen = read_length(g, rid);
		rdir = p1->rdir;
		if(rdir){
			roff = p1->roff + p1->rlen;
			plen = rlen - p1->roff;
		} else {
			roff = p1->roff;
			plen = p1->roff + p1->rlen;
		}
		for(j=1;j<n_eid;j++){
			e2 = ref_elist(g->edges, eids[j]);
			if(plen < rlen && find_readpath_edge(g, eids[j], 0, rid, rdir, rlen, plen, &len, &pidx)){
				one_bitvec(flags[j], pidx);
				plen += len;
				p2 = ref_plist(e2->paths, pidx);
				nc = cat_cigars(g->aux_cigars[0], nc, p2->cigars, p2->n_cigar);
			} else {
				nc = append_cigars(g->aux_cigars[0], nc, ALN_CIGAR_TYPE_DEL, e2->len);
			}
		}
		p->rid  = rid;
		p->rdir = rdir;
		if(rdir){
			p->roff = (rlen - plen);
			p->rlen = roff - p->roff;
		} else {
			p->roff = roff;
			p->rlen = plen - roff;
		}
		p->cigars = malloc(sizeof(AlnCigar) * nc);
		p->n_cigar = cat_cigars(p->cigars, 0, g->aux_cigars[0], nc);
	}
	e1  = ref_elist(g->edges, eids[n_eid-1]);
	for(i=0;i<count_plist(e1->paths);i++){
		if(get_bitvec(flags[n_eid - 1], i)) continue;
		p1 = ref_plist(e1->paths, i);
		p  = next_ref_plist(e->paths);
		nc = 0;
		for(k=0;k<p1->n_cigar;k++) nc = append_cigars(g->aux_cigars[0], nc, p1->cigars[p1->n_cigar - k - 1].type, p1->cigars[p1->n_cigar - k -1].len);
		rid  = p1->rid;
		rlen = read_length(g, rid);
		rdir = 1 ^ p1->rdir;
		if(rdir){
			roff = p1->roff + p1->rlen;
			plen = rlen - p1->roff;
		} else {
			roff = p1->roff;
			plen = p1->roff + p1->rlen;
		}
		for(j=n_eid-2;;j--){
			e2 = ref_elist(g->edges, eids[j]);
			if(plen < rlen && find_readpath_edge(g, eids[j], 1, rid, rdir, rlen, plen, &len, &pidx)){
				one_bitvec(flags[j], pidx);
				plen += len;
				p2 = ref_plist(e2->paths, pidx);
				for(k=0;k<p2->n_cigar;k++) nc = append_cigars(g->aux_cigars[0], nc, p2->cigars[p2->n_cigar - k - 1].type, p2->cigars[p2->n_cigar - k -1].len);
			} else {
				nc = append_cigars(g->aux_cigars[0], nc, ALN_CIGAR_TYPE_DEL, e2->len);
			}
			if(j == 0) break;
		}
		p->rid  = rid;
		p->rdir = 1 ^ rdir;
		if(rdir){
			p->roff = (rlen - plen);
			p->rlen = roff - p->roff;
		} else {
			p->roff = roff;
			p->rlen = plen - roff;
		}
		reverse_cigars(g->aux_cigars[0], nc);
		p->cigars = malloc(sizeof(AlnCigar) * nc);
		p->n_cigar = cat_cigars(p->cigars, 0, g->aux_cigars[0], nc);
	}
	sort_edge_paths(g, dst);
	add_edge_v2v(g, dst, ref_elist(g->edges, eids[0])->vid1, ref_elist(g->edges, eids[0])->vdir1,
			ref_elist(g->edges, eids[n_eid-1])->vid2, ref_elist(g->edges, eids[n_eid-1])->vdir2);
	map_r2v_edge(g, dst);
	if(new_eid) *new_eid = dst;
	return 1;
}

int connect_e2e2e(SeqGraph *g, uint32_t eid1, uint32_t eid2, uint32_t eid3){
	Edge *e1, *e2, *e3;
	uint32_t edir1, edir2, edir3, eids[3],eid4;
	e1 = ref_elist(g->edges, eid1);
	e2 = ref_elist(g->edges, eid2);
	e3 = ref_elist(g->edges, eid3);
	if(e1->vid1 == e1->vid2) return 0;
	if(e3->vid1 == e3->vid2) return 0;
	edir2 = 0;
	if(e1->vid1 == e2->vid1) edir1 = 1;
	else if(e1->vid2 == e1->vid1) edir1 = 0;
	else return 0;
	if(e2->vid2 == e3->vid1) edir3 = 0;
	else if(e2->vid2 == e3->vid2) edir3 = 1;
	else return 0;
	if((edir1? e1->vdir1 : e1->vdir2) == (edir2? e2->vdir2 : e2->vdir1)) return 0;
	if((edir2? e2->vdir1 : e2->vdir2) == (edir3? e3->vdir2 : e3->vdir1)) return 0;
	if(edir1) reverse_edge(g, eid1);
	if(edir3) reverse_edge(g, eid3);
	encap_bitvec(g->flags[0], count_plist(e1->paths));
	zeros_bitvec(g->flags[0]);
	encap_bitvec(g->flags[2], count_plist(e3->paths));
	zeros_bitvec(g->flags[2]);
	eids[0] = eid1;
	eids[1] = eid2;
	eids[2] = eid3;
	if(connect_edges(g, eids, g->flags, 3, &eid4) == 0) return 0;
	sort_edge_paths(g, eid4);
	e1 = ref_elist(g->edges, eid1);
	e3 = ref_elist(g->edges, eid3);
	add_edge_v2v(g, eid4, e1->vid1, e1->vdir1, e3->vid2, e3->vdir2);
	remove_edge(g, eid1);
	remove_edge(g, eid3);
	return 1;
}

//TODO -- Graph --

SeqGraph* init_seqgraph(FileReader *rds, int fix_rd_len){
	SeqGraph *g;
	ReadPath *p;
	Edge *edge;
	Vertex *v1, *v2;
	uint32_t rid, len, eid, vid1, vid2;
	uint64_t offset;
	int n;
	char *seq;
	Sequence *read;
	g = malloc(sizeof(SeqGraph));
	g->edges   = init_elist(1024);
	g->verts   = init_vlist(1024);
	g->ecache  = init_u32list(12);
	g->vcache  = init_u32list(12);
	g->v_flags = init_bitvec(1024);
	g->e_flags = init_bitvec(1024);
	g->e_chgs  = init_bitvec(1024);
	g->r2v_dir = init_bitvec(1024);
	g->r2v_map = init_u32list(1024);
	g->n_rd    = 0;
	g->rd_len  = 0;
	g->max_rd_len = 0;
	g->rd_seqs = init_u64list(1024);
	g->rd_lens = init_u16list(1024);
	g->rd_seq_offs = init_u64list(1024);
	g->ols     = init_olist(1024);
	g->aux_cigars[0] = malloc(sizeof(AlnCigar) * 1024 * 1024);
	g->aux_cigars[1] = malloc(sizeof(AlnCigar) * 1024 * 1024);
	g->aux_cigars[2] = malloc(sizeof(AlnCigar) * 1024 * 1024);
	g->flags[0] = init_bitvec(32);
	g->flags[1] = init_bitvec(32);
	g->flags[2] = init_bitvec(32);
	offset = 0;
	rid = 0;
	read = NULL;
	while(fread_fasta(&read, rds)){
		len = read->seq.size;
		seq = read->seq.string;
		if(fix_rd_len && g->rd_len == 0) g->rd_len = len;
		if(len > g->max_rd_len) g->max_rd_len = len;
		vid1 = lend_vertex(g);
		vid2 = lend_vertex(g);
		v1   = ref_vlist(g->verts, vid1);
		v2   = ref_vlist(g->verts, vid2);
		eid  = lend_edge(g);
		edge = ref_elist(g->edges, eid);
		edge->vid1  = vid1;
		edge->vid2  = vid2;
		edge->vdir1 = 0;
		edge->vdir2 = 1;
		edge->len   = len;
		n = strlen(seq);
		p = next_ref_plist(edge->paths);
		p->rid     = rid;
		p->roff    = 0;
		p->rdir    = 0;
		p->rlen    = n;
		p->n_cigar = 1;
		p->cigars  = malloc(sizeof(AlnCigar) * 1);
		p->cigars[0].len  = n;
		p->cigars[0].type = ALN_CIGAR_TYPE_MAT;
		push_u32v(v1, eid);
		push_u32v(v2, eid);
		push_u32list(g->r2v_map, vid1);
		zero2bitvec(g->r2v_dir);
		push_u32list(g->r2v_map, vid2);
		one2bitvec(g->r2v_dir);
		encap_u64list(g->rd_seqs, (offset + len + 31) / 32);
		seq2bits(as_array_u64list(g->rd_seqs), offset, seq, len);
		if(!fix_rd_len){
			if(count_u16list(g->rd_lens) % 16 == 0) push_u64list(g->rd_seq_offs, offset);
			push_u16list(g->rd_lens, len);
		}
		offset += len;
		rid ++;
	}
	g->n_rd = rid;
	return g;
}

void free_seqgraph(SeqGraph *g){
	uint32_t i, j;
	Vertex *v;
	Edge *e;
	ReadPath *p;
	for(i=0;i<count_elist(g->edges);i++){
		e = ref_elist(g->edges, i);
		if(e->paths){
			for(j=0;j<count_plist(e->paths);j++){
				p = ref_plist(e->paths, j);
				if(p->cigars) free(p->cigars);
			}
			free_plist(e->paths);
		}
	}
	free_elist(g->edges);
	for(i=0;i<count_vlist(g->verts);i++){
		v = ref_vlist(g->verts, i);
		u32v_free(v);
	}
	free_vlist(g->verts);
	free_u32list(g->ecache);
	free_u32list(g->vcache);
	free_bitvec(g->v_flags);
	free_bitvec(g->e_flags);
	free_bitvec(g->e_chgs);
	free_bitvec(g->r2v_dir);
	free_u32list(g->r2v_map);
	free_u64list(g->rd_seqs);
	free_u16list(g->rd_lens);
	free_u64list(g->rd_seq_offs);
	free_olist(g->ols);
	free(g->aux_cigars[0]);
	free(g->aux_cigars[1]);
	free(g->aux_cigars[2]);
	free_bitvec(g->flags[0]);
	free_bitvec(g->flags[1]);
	free_bitvec(g->flags[2]);
	free(g);
}

//TODO  -- Output --

void output_vertexs_seqgraph(SeqGraph *g, FILE *out){
	Vertex *v;
	uint32_t i, j;
	for(i=0;i<count_vlist(g->verts);i++){
		if(get_bitvec(g->v_flags, i) == 0) continue;
		v = ref_vlist(g->verts, i);
		fprintf(out, "V\t%u\tT\t%u\t", i, count_u32v(v));
		for(j=0;j<count_u32v(v);j++){
			fprintf(out, "%u,", get_u32v(v, j));
		}
		fprintf(out, "\n");
	}
}

void output_edges_seqgraph(SeqGraph *g, FILE *out){
	Edge *e;
	ReadPath *p;
	u8list *seq, *aln, *cig;
	uint32_t i, j, m;
	seq = init_u8list(1024);
	aln = init_u8list(1024);
	cig = init_u8list(1024);
	for(i=0;i<count_elist(g->edges);i++){
		if(get_bitvec(g->e_flags, i) == 0) continue;
		e = ref_elist(g->edges, i);
		fprintf(out, "E\t%u\t%u\t%u\t%c\t%u\t%c\t%u\n", i, e->len, e->vid1, "+-"[e->vdir1], e->vid2, "+-"[e->vdir2], count_plist(e->paths));
		encap_u8list(aln, e->len + 1);
		for(j=0;j<count_plist(e->paths);j++){
			p = ref_plist(e->paths, j);
			encap_u8list(seq, p->rlen + 1);
			if(p->rdir){
				bits2revseq((char*)as_array_u8list(seq), as_array_u64list(g->rd_seqs), read_offset(g, p->rid) + p->roff, p->rlen);
			} else {
				bits2seq((char*)as_array_u8list(seq), as_array_u64list(g->rd_seqs), read_offset(g, p->rid) + p->roff, p->rlen);
			}
			encap_u8list(cig, p->n_cigar * 10);
			cigars2string(p->cigars, p->n_cigar, (char*)as_array_u8list(cig));
			m = cigars_seq2aln((char*)as_array_u8list(aln), p->cigars, p->n_cigar, 0, (char*)as_array_u8list(seq));
			set_u8list(aln, m, 0);
			fprintf(out, "R\t%u\t%c\t%u\t%u\t%s\n", p->rid, "+-"[p->rdir], p->roff, p->rlen, (char*)as_array_u8list(cig));
			fprintf(out, "S\t%s\n", as_array_u8list(aln));
		}
	}
	free_u8list(seq);
	free_u8list(aln);
	free_u8list(cig);
}

void output_dot_seqgraph(SeqGraph *g, FILE *out){
	Vertex *v;
	Edge *e;
	uint32_t i, j, f;
	uint32_t eid;
	fprintf(out, "digraph {\n");
	for(i=0;i<count_vlist(g->verts);i++){
		v = ref_vlist(g->verts, i);
		if(get_bitvec(g->v_flags, i) == 0) continue;
		fprintf(out, "V%u [label=\"{V%u | {", i, i);
		f = 0;
		for(j=0;j<count_u32v(v);j++){
			eid = get_u32v(v, j);
			e = ref_elist(g->edges, eid);
			if(e->vid1 == i){
				if(f) fprintf(out, " | ");
				f ++;
				fprintf(out, "<E%u%c> E%u%c", eid, "FR"[e->vdir1], eid, "+-"[e->vdir1]);
			}
			if(e->vid2 == i){
				if(f) fprintf(out, " | ");
				f ++;
				fprintf(out, "<E%u%c> E%u%c", eid, "FR"[e->vdir2], eid, "+-"[e->vdir2]);
			}
		}
		fprintf(out, "}}\" shape=record color=blue]\n");
	}
	for(i=0;i<count_elist(g->edges);i++){
		if(get_bitvec(g->e_flags, i) == 0) continue;
		e = ref_elist(g->edges, i);
		fprintf(out, "V%d:E%u%c -> V%d:E%u%c [label=\"E%u:%dbp\"]\n", e->vid1, i, "FR"[e->vdir1], e->vid2, i, "FR"[e->vdir2], i, e->len);
		/*
		for(j=0;j<count_plist(e->paths);j++){
			p = ref_plist(e->paths, j);
			fprintf(out, "V%d:E%u%c -> V%d:E%u%c [label=\"R%d%c:%d-%dbp\" color=yellow]\n", e->vid1, i, "FR"[e->vdir1], e->vid2, i, "FR"[e->vdir2],
				p->rid, "+-"[p->rdir], p->roff, p->roff + p->rlen);
		}
		*/
	}
	fprintf(out, "}\n");
	fflush(out);
}

void disp_vertex_readpaths(SeqGraph *g, uint32_t vid){
	Vertex *v;
	Edge *e;
	ReadPath *p;
	uint32_t i, eid, pid;
	v = ref_vlist(g->verts, vid);
	for(i=0;i<count_u32v(v);i++){
		eid = get_u32v(v, i);
		e = ref_elist(g->edges, eid);
		fprintf(stdout, "V%u E%u %u V%u%c -> V%u%c\n", vid, eid, e->len, e->vid1, "+-"[e->vdir1], e->vid2, "+-"[e->vdir2]);
		for(pid=0;pid<count_plist(e->paths);pid++){
			p = ref_plist(e->paths, pid);
			if(p->rid != 1) continue;
			fprintf(stdout, "- V%u E%u R%u %c %u+%u\n", vid, eid, p->rid, "+-"[p->rdir], p->roff, p->rlen);
		}
	}
	fflush(stdout);
}

void disp_overlap(SeqGraph *g, Overlap *ol){
	ReadTracker RT, *rt;
	uint32_t gend, goff;
	rt = &RT;
	init_readtracker(rt, g, ol->rids[0], ol->rdirs[0], ol->rlens[0]);
	while(next_readtracker(rt)){
		goff = select_cigars_seqlen(ol->cigars, ol->n_cigar, rt->off, 0);
		gend = select_cigars_seqlen(ol->cigars, ol->n_cigar, rt->off + rt->len, 0);
		fprintf(stdout, "A path   %u%c:%u+%u(%u)[%u%c -- %u%c(%u) -> %u%c]\n", rt->rid, "+-"[rt->rdir], rt->off, rt->len, rt->rlen,
			rt->vid, "+-"[rt->vdir], rt->eid, "+-"[rt->edir], ref_elist(g->edges, rt->eid)->len, rt->nvid, "+-"[rt->nvdir]);
		fprintf(stdout, "A region %u%c[%u -> %u]\n", rt->rid, "+-"[rt->rdir], goff, gend);
	}
	init_readtracker(rt, g, ol->rids[1], ol->rdirs[1], ol->rlens[1]);
	while(next_readtracker(rt)){
		goff = select_cigars_seqlen(ol->cigars, ol->n_cigar, rt->off, 1);
		gend = select_cigars_seqlen(ol->cigars, ol->n_cigar, rt->off + rt->len, 1);
		fprintf(stdout, "B path   %u%c:%u+%u(%u)[%u%c -- %u%c(%u) -> %u%c]\n", rt->rid, "+-"[rt->rdir], rt->off, rt->len, rt->rlen,
			rt->vid, "+-"[rt->vdir], rt->eid, "+-"[rt->edir], ref_elist(g->edges, rt->eid)->len, rt->nvid, "+-"[rt->nvdir]);
		fprintf(stdout, "B region %u%c[%u -> %u]\n", rt->rid, "+-"[rt->rdir], goff, gend);
	}
}

void disp_vertex(SeqGraph *g, uint32_t vid){
	Vertex *v;
	Edge *e;
	ReadPath *p;
	uint32_t i, j;
	v = ref_vlist(g->verts, vid);
	fprintf(stdout, "DISP_VERTEX V%u\n", vid);
	for(i=0;i<count_u32v(v);i++){
		e = ref_elist(g->edges, get_u32v(v, i));
		fprintf(stdout, " E%u:%ubp V%u%c -> V%u%c\n", get_u32v(v, i), e->len, e->vid1, "+-"[e->vdir1], e->vid2, "+-"[e->vdir2]);
		for(j=0;j<count_plist(e->paths);j++){
			p = ref_plist(e->paths, j);
			fprintf(stdout, "  R%u%c %u+%u\n", p->rid, "+-"[p->rdir], p->roff, p->rlen);
		}
	}
	fprintf(stdout, "END DISP_VERTEX\n");
}

void disp_edge(SeqGraph *g, uint32_t eid){
	Edge *e;
	ReadPath *p;
	uint32_t j;
	e = ref_elist(g->edges, eid);
	fprintf(stdout, " E%u:%ubp V%u%c -> V%u%c\n", eid, e->len, e->vid1, "+-"[e->vdir1], e->vid2, "+-"[e->vdir2]);
	for(j=0;j<count_plist(e->paths);j++){
		p = ref_plist(e->paths, j);
		fprintf(stdout, "  R%u%c %u+%u\n", p->rid, "+-"[p->rdir], p->roff, p->rlen);
	}
}

void disp_readtracker(ReadTracker *rt){
	fprintf(stdout, "RT rid  = %u\n", rt->rid);
	fprintf(stdout, "RT rdir = %u\n", rt->rdir);
	fprintf(stdout, "RT rlen = %u\n", rt->rlen);
	fprintf(stdout, "RT plen = %u\n", rt->plen);
	fprintf(stdout, "RT vid  = %u\n", rt->vid);
	fprintf(stdout, "RT vdir = %c\n", "+-"[rt->vdir]);
	fprintf(stdout, "RT nvid = %u\n", rt->nvid);
	fprintf(stdout, "RT nvdir= %c\n", "+-"[rt->nvdir]);
	fprintf(stdout, "RT eid  = %u\n", rt->eid);
	fprintf(stdout, "RT edir = %c\n", "+-"[rt->edir]);
	fprintf(stdout, "RT pid  = %u\n", rt->pid);
	fprintf(stdout, "RT off  = %u\n", rt->off);
	fprintf(stdout, "RT len  = %u\n", rt->len);
	fflush(stdout);
}

//TODO -- Checking --

void check_seqgraph_link_refs(SeqGraph *g){
	Edge *e;
	uint32_t i;
	for(i=0;i<count_elist(g->edges);i++){
		if(get_bitvec(g->e_flags, i) == 0) continue;
		e = ref_elist(g->edges, i);
		if(get_bitvec(g->v_flags, e->vid1) == 0){
			fprintf(stdout, " -- E%u:V%u lost in %s -- %s:%d --\n", i, e->vid1, __FUNCTION__, __FILE__, __LINE__);
			fflush(stdout);
			abort();
		}
		if(get_bitvec(g->v_flags, e->vid2) == 0){
			fprintf(stdout, " -- E%u:V%u lost in %s -- %s:%d --\n", i, e->vid2, __FUNCTION__, __FILE__, __LINE__);
			fflush(stdout);
			abort();
		}
	}
}

void check_readpath(SeqGraph *g, uint32_t rid){
	ReadTracker *rt, RT;
	rt = &RT;
	init_readtracker(rt, g, rid, 0, read_length(g, rid));
	while(next_readtracker(rt));
	if(rt->plen != rt->rlen){
		fprintf(stdout, "BAD READ PATH %u/%u (%u != %u)\n", rid, g->n_rd, rt->plen, rt->rlen);
		disp_readtracker(rt);
		disp_vertex(g, rt->vid);
		abort();
	}
}

void check_seqgraph_readpaths(SeqGraph *g){
	uint32_t rid;
	for(rid=0;rid<g->n_rd;rid++){
		check_readpath(g, rid);
	}
}

void check_readpaths(SeqGraph *g, uint32_t rid){
	ReadTracker *rt, RT;
	rt = &RT;
	init_readtracker(rt, g, rid, 0, read_length(g, rid));
	while(next_readtracker(rt));
	if(rt->plen != rt->rlen){
		fprintf(stdout, "BAD READ PATH %u/%u (%u != %u)\n", rid, rt->rid, rt->plen, rt->rlen);
		disp_readtracker(rt);
		abort();
	}
}

void check_link(SeqGraph *g, uint32_t vid1, uint32_t vid2){
	int a, b;
	a = has_this_link(g, vid1, vid2);
	b = has_this_link(g, vid2, vid1);
	if(a != b){
		fprintf(stdout, " -- %u and %u [%d, %d] in %s -- %s:%d --\n", vid1, vid2, a, b, __FUNCTION__, __FILE__, __LINE__);
		fflush(stdout);
		abort();
	}
}

void check_edge_link(SeqGraph *g, uint32_t eid){
	Edge *e;
	Vertex *v1, *v2;
	e = ref_elist(g->edges, eid);
	v1 = ref_vlist(g->verts, e->vid1);
	v2 = ref_vlist(g->verts, e->vid2);
	if(locate_u32v(v1, eid, 0) == count_u32v(v1)){
#ifdef DEBUG
		fprintf(stdout, " -- %u miss eid %u in %s -- %s:%d --\n", e->vid1, eid, __FUNCTION__, __FILE__, __LINE__);
		fflush(stdout);
		abort();
#endif
	}
	if(locate_u32v(v2, eid, 0) == count_u32v(v2)){
#ifdef DEBUG
		fprintf(stdout, " -- %u miss eid %u in %s -- %s:%d --\n", e->vid2, eid, __FUNCTION__, __FILE__, __LINE__);
		fflush(stdout);
		abort();
#endif
	}
}

void check_seqgraph_edges(SeqGraph *g){
	Vertex *v;
	Edge *e;
	uint32_t i;
	for(i=0;i<count_elist(g->edges);i++){
		if(get_bitvec(g->e_flags, i) == 0) continue;
		e = ref_elist(g->edges, i);
		if(get_bitvec(g->v_flags, e->vid1) == 0){
			fprintf(stdout, " -- E%u point to obsoleted V%u in %s -- %s:%d --\n", i, e->vid1, __FUNCTION__, __FILE__, __LINE__);
			fflush(stdout);
			abort();
		}
		v = ref_vlist(g->verts, e->vid1);
		if(locate_u32v(v, i, 0) == count_u32v(v)){
			fprintf(stdout, " -- E%u , V%u miss eid, in %s -- %s:%d --\n", i, e->vid1, __FUNCTION__, __FILE__, __LINE__);
			fflush(stdout);
			abort();
		}
		if(get_bitvec(g->v_flags, e->vid2) == 0){
			fprintf(stdout, " -- E%u point to obsoleted V%u in %s -- %s:%d --\n", i, e->vid2, __FUNCTION__, __FILE__, __LINE__);
			fflush(stdout);
			abort();
		}
		v = ref_vlist(g->verts, e->vid2);
		if(locate_u32v(v, i, 0) == count_u32v(v)){
			fprintf(stdout, " -- E%u , V%u miss eid, in %s -- %s:%d --\n", i, e->vid2, __FUNCTION__, __FILE__, __LINE__);
			fflush(stdout);
			abort();
		}
	}
}

void check_edge_length(SeqGraph *g, uint32_t eid){
	Edge *e;
	ReadPath *p;
	uint32_t i;
	int tlen, rlen;
	e = ref_elist(g->edges, eid);
	for(i=0;i<count_plist(e->paths);i++){
		p = ref_plist(e->paths, i);
		cigars_lengths(p->cigars, p->n_cigar, &tlen, &rlen, NULL);
		if(tlen != (int)e->len || rlen != (int)p->rlen){
			fprintf(stdout, " -- E%u:%ubp R%u:%ubp but cigars %d %d  in %s -- %s:%d --\n", eid, e->len, p->rid, p->rlen, tlen, rlen, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	}
	zero_bitvec(g->e_chgs, eid);
}

void check_seqgraph_edge_lengths(SeqGraph *g){
	uint32_t eid;
	for(eid=0;eid<count_elist(g->edges);eid++){
		if(get_bitvec(g->e_flags, eid) && get_bitvec(g->e_chgs, eid)) check_edge_length(g, eid);
	}
}

//TODO -- Actions --

void overlapping_core(SeqGraph *g, Overlap *ol){
	ReadTracker RTS[2], *rts[2];
	uint32_t gends[2], goffs[2];
	uint32_t eids[2], edirs[2], flags[2], neid;
#ifdef DEBUG
	disp_overlap(g, ol);
	fflush(stdout);
#endif
	rts[0] = &RTS[0];
	rts[1] = &RTS[1];
	init_readtracker(rts[0], g, ol->rids[0], ol->rdirs[0], ol->rlens[0]);
	init_readtracker(rts[1], g, ol->rids[1], ol->rdirs[1], ol->rlens[1]);
	flags[0] = flags[1] = 1;
	while(1){
		if(flags[0]){
			if(next_readtracker(rts[0]) == 0){
#ifdef DEBUG
				fprintf(stdout, "BREAK A %u%c[%u <=> %u]\n", rts[0]->rid, "+-"[rts[0]->rdir], rts[0]->plen, rts[0]->rlen);
				fflush(stdout);
#endif
				break;
			}
			align_read_path(g, ol, 0, rts[0]);
		}
		if(flags[1]){
			if(next_readtracker(rts[1]) == 0){
#ifdef DEBUG
				fprintf(stdout, "BREAK B %u%c[%u <=> %u]\n", rts[1]->rid, "+-"[rts[1]->rdir], rts[1]->plen, rts[1]->rlen);
				fflush(stdout);
#endif
				break;
			}
			align_read_path(g, ol, 1, rts[1]);
		}
		flags[0] = flags[1] = 1;
		goffs[0] = select_cigars_seqlen(ol->cigars, ol->n_cigar, rts[0]->off, 0);
		goffs[1] = select_cigars_seqlen(ol->cigars, ol->n_cigar, rts[1]->off, 1);
		gends[0] = select_cigars_seqlen(ol->cigars, ol->n_cigar, rts[0]->off + rts[0]->len, 0);
		gends[1] = select_cigars_seqlen(ol->cigars, ol->n_cigar, rts[1]->off + rts[1]->len, 1);
#ifdef DEBUG
		Edge *e1, *e2;
		fprintf(stdout, "A%u%c:%u+%u(%u)[%u%c -- %u%c(%u) -> %u%c]\n", rts[0]->rid, "+-"[rts[0]->rdir], rts[0]->off, rts[0]->len, rts[0]->rlen,
			rts[0]->vid, "+-"[rts[0]->vdir], rts[0]->eid, "+-"[rts[0]->edir], ref_elist(g->edges, rts[0]->eid)->len, rts[0]->nvid, "+-"[rts[0]->nvdir]);
		fprintf(stdout, "regions A %u%c[%u -> %u]\n", rts[0]->rid, "+-"[rts[0]->rdir], goffs[0], gends[0]);
		e1 = ref_elist(g->edges, rts[0]->eid);
		fprintf(stdout, "edge %u%c [%u%c -> %u%c]\n", rts[0]->eid, "+-"[rts[0]->edir], e1->vid1, "+-"[e1->vdir1], e1->vid2, "+-"[e1->vdir2]);
		fprintf(stdout, "B%u%c:%u+%u(%u)[%u%c -- %u%c(%u) -> %u%c]\n", rts[1]->rid, "+-"[rts[1]->rdir], rts[1]->off, rts[1]->len, rts[1]->rlen,
			rts[1]->vid, "+-"[rts[1]->vdir], rts[1]->eid, "+-"[rts[1]->edir], ref_elist(g->edges, rts[1]->eid)->len, rts[1]->nvid, "+-"[rts[1]->nvdir]);
		fprintf(stdout, "regions B %u%c[%u -> %u]\n", rts[1]->rid, "+-"[rts[1]->rdir], goffs[1], gends[1]);
		e2 = ref_elist(g->edges, rts[1]->eid);
		fprintf(stdout, "edge %u%c [%u%c -> %u%c]\n", rts[1]->eid, "+-"[rts[1]->edir], e2->vid1, "+-"[e2->vdir1], e2->vid2, "+-"[e2->vdir2]);
		fflush(stdout);
#endif
		if(gends[0] <= goffs[1]){
			flags[1] = 0;
#ifdef DEBUG
			fprintf(stdout, "Next A\n");
			fflush(stdout);
#endif
		} else if(gends[1] <= goffs[0]){
			flags[0] = 0;
#ifdef DEBUG
			fprintf(stdout, "Next B\n");
			fflush(stdout);
#endif
		} else if(goffs[0] < goffs[1]){
			insert_v2e(g, rts[0]->eid, rts[0]->edir, rts[1]->vid, rts[1]->vdir, goffs[1] - goffs[0]);
			rts[0]->nvid  = rts[1]->vid; rts[0]->nvdir = rts[1]->vdir; rts[0]->plen  = rank_cigars_seqlen(ol->cigars, ol->n_cigar, goffs[1], 0);
			rts[1]->nvid  = rts[1]->vid; rts[1]->nvdir = rts[1]->vdir; rts[1]->plen  = rts[1]->off;
#ifdef DEBUG
			fprintf(stdout, "INSERT B -> A\n");
			fflush(stdout);
#endif
		} else if(goffs[1] < goffs[0]){
			insert_v2e(g, rts[1]->eid, rts[1]->edir, rts[0]->vid, rts[0]->vdir, goffs[0] - goffs[1]);
			rts[1]->nvid  = rts[0]->vid; rts[1]->nvdir = rts[0]->vdir; rts[1]->plen  = rank_cigars_seqlen(ol->cigars, ol->n_cigar, goffs[0], 1);
			rts[0]->nvid  = rts[0]->vid; rts[0]->nvdir = rts[0]->vdir; rts[0]->plen  = rts[0]->off;
#ifdef DEBUG
			fprintf(stdout, "INSERT A -> B\n");
			fflush(stdout);
#endif
		} else {
			eids[0]  = rts[0]->eid;
			eids[1]  = rts[1]->eid;
			edirs[0] = rts[0]->edir;
			edirs[1] = rts[1]->edir;
#ifdef DEBUG
			e1 = ref_elist(g->edges, eids[0]);
			e2 = ref_elist(g->edges, eids[1]);
			fprintf(stdout, " -- Will TRANSFER %u %u %c\n", rts[0]->vid, rts[1]->vid, "+-"[rts[0]->vdir ^ rts[1]->vdir]);
			fprintf(stdout, " -- Will MERGE %d[%d%c,%d%c] %d[%d%c,%d%c] '%c' --\n", eids[0], e1->vid1, "+-"[e1->vdir1], e1->vid2, "+-"[e1->vdir2],
				eids[1], e2->vid1, "+-"[e2->vdir1], e2->vid2, "+-"[e2->vdir2], "+-"[edirs[0] ^ edirs[1]]);
			fflush(stdout);
#endif
			if(transfer_v2v(g, rts[1]->vid, rts[0]->vid, rts[0]->vdir ^ rts[1]->vdir)){
#ifdef DEBUG
				fprintf(stdout, "TRANSFER B -> A\n");
				fflush(stdout);
#endif
				//if(rts[1]->nvid == rts[1]->vid){
					//rts[1]->nvid  = rts[0]->vid;
					//rts[1]->nvdir = rts[1]->nvdir ^ rts[0]->vdir ^ rts[1]->vdir;
				//}
				//if(rts[0]->nvid == rts[1]->vid){
					//rts[0]->nvid  = rts[0]->vid;
					//rts[0]->nvdir = rts[0]->nvdir ^ rts[0]->vdir ^ rts[1]->vdir;
				//}
				rts[1]->nvid  = rts[0]->nvid  = rts[0]->vid;
				rts[1]->nvdir = rts[0]->nvdir = rts[0]->vdir;
				rts[0]->plen  = rts[0]->off;
				rts[1]->plen  = rts[1]->off;
			}
			else {
				if(gends[0] < gends[1]){
					neid = insert_v2e(g, rts[1]->eid, rts[1]->edir, rts[0]->nvid, rts[0]->nvdir, gends[0] - goffs[1]);
					if(rts[1]->edir) eids[1] = neid;
					rts[1]->nvid  = rts[0]->nvid;
					rts[1]->nvdir = rts[0]->nvdir;
					rts[1]->plen  = rank_cigars_seqlen(ol->cigars, ol->n_cigar, gends[0], 1);
#ifdef DEBUG
					fprintf(stdout, "INSERT   A' -> B\n");
					fflush(stdout);
#endif
				} else if(gends[1] < gends[0]){
					neid = insert_v2e(g, rts[0]->eid, rts[0]->edir, rts[1]->nvid, rts[1]->nvdir, gends[1] - goffs[0]);
					if(rts[0]->edir) eids[0] = neid;
					rts[0]->nvid  = rts[1]->nvid;
					rts[0]->nvdir = rts[1]->nvdir;
					rts[0]->plen  = rank_cigars_seqlen(ol->cigars, ol->n_cigar, gends[1], 0);
#ifdef DEBUG
					fprintf(stdout, "INSERT   B' -> A\n");
					fflush(stdout);
#endif
				} else {
					if(transfer_v2v(g, rts[1]->nvid, rts[0]->nvid, rts[0]->nvdir ^ rts[1]->nvdir)){
						rts[1]->nvid  = rts[0]->nvid;
						rts[1]->nvdir = rts[0]->nvdir;
					}
#ifdef DEBUG
					fprintf(stdout, "TRANSFER B' -> A'\n");
					fflush(stdout);
#endif
				}
				merge_equal_edges(g, eids[0], eids[1], edirs[0] ^ edirs[1]);
			}
		}
#ifdef DEBUG
		fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		fflush(stdout);
#endif
	}
}

void overlapping_seqgraph(SeqGraph *g, FileReader *ols){
	Overlap OL, *ol;
	uint32_t col;
	int lens[3], offs[2], ends[2], x, y, n;
	OL.cigars = malloc(sizeof(AlnCigar) * 1024);
	ol = &OL;
	while((n = fread_table(ols)) != -1){
		if(n < 5) continue;
		if(ols->line->string[0] == '#') continue;
		fprintf(stdout, "OLS %s\n", ols->line->string);
		fflush(stdout);

		col = 0;

		ol->rids[0]  = atoll(get_col_str(ols, col++));
		ol->rdirs[0] = (get_col_str(ols, col++)[0] == '-');
		ol->rids[1]  = atoll(get_col_str(ols, col++));
		ol->rdirs[1] = (get_col_str(ols, col++)[0] == '-');

		ol->n_cigar = string2cigars(ol->cigars, get_col_str(ols, col), get_col_len(ols, col));
		cigars_lengths(ol->cigars, ol->n_cigar, lens, lens + 1, lens + 2);
		ol->rlens[0] = read_length(g, ol->rids[0]);
		ol->rlens[1] = read_length(g, ol->rids[1]);
		if(lens[1] != (int)ol->rlens[0]){
			fprintf(stderr, " -- Inconsistent read length \"%s\", read[%d] %d != %d in %s -- %s:%d --\n",
				ols->line->string, ol->rids[0], ol->rlens[0], lens[1],  __FUNCTION__, __FILE__, __LINE__);
			continue;
		}
		if(lens[2] != (int)ol->rlens[1]){
			fprintf(stderr, " -- Inconsistent read length \"%s\", read[%d] %d != %d in %s -- %s:%d --\n",
				ols->line->string, ol->rids[1], ol->rlens[1], lens[2],  __FUNCTION__, __FILE__, __LINE__);
			continue;
		}
		offs[0] = select_cigars_seqlen(ol->cigars, ol->n_cigar, 0, 0);
		offs[1] = select_cigars_seqlen(ol->cigars, ol->n_cigar, 0, 1);
		ends[0] = select_cigars_seqlen(ol->cigars, ol->n_cigar, ol->rlens[0] - 1, 0);
		ends[1] = select_cigars_seqlen(ol->cigars, ol->n_cigar, ol->rlens[1] - 1, 1);
		x = (offs[0] > offs[1])? offs[0] : offs[1];
		y = (ends[0] < ends[1])? ends[0] : ends[1];
		if(x > y) continue;
		ol->olens[0] = rank_cigars_seqlen(ol->cigars, ol->n_cigar, y, 0) - rank_cigars_seqlen(ol->cigars, ol->n_cigar, x, 0) + 1;
		ol->olens[1] = rank_cigars_seqlen(ol->cigars, ol->n_cigar, y, 1) - rank_cigars_seqlen(ol->cigars, ol->n_cigar, x, 1) + 1;
		ol->rgaps[0] = y - x + 1 - ol->olens[0];
		ol->rgaps[1] = y - x + 1 - ol->olens[1];

		overlapping_core(g, ol);
		//fprintf(stdout, " -- %uE %uV in %s -- %s:%d --\n", (uint32_t)count_elist(g->edges), (uint32_t)count_vlist(g->verts), __FUNCTION__, __FILE__, __LINE__);
		//fflush(stdout);
	}
	free(ol->cigars);
}

int simplify_vertex_core(SeqGraph *g, uint32_t vid, u32v *fwds, u32v *revs){
	Vertex *v;
	Edge *e;
	uint32_t eids[2];
	uint32_t i, eid;
	int ret;
	v = ref_vlist(g->verts, vid);
	clear_u32v(fwds);
	clear_u32v(revs);
	ret = 0;
	for(i=0;i<count_u32v(v);i++){
		eid = get_u32v(v, i);
		e = ref_elist(g->edges, eid);
		if(e->vid1 == e->vid2) return 0;
		if(e->vid1 == vid) push_u32v(e->vdir1? revs : fwds, eid);
		if(e->vid2 == vid) push_u32v(e->vdir2? revs : fwds, eid);
	}
	if(count_u32v(fwds) == 1 && count_u32v(revs) == 1){
		eids[0] = get_u32v(fwds, 0);
		e = ref_elist(g->edges, eids[0]);
		if(e->vid1 == vid) reverse_edge(g, eids[0]);
		encap_bitvec(g->flags[0], count_plist(e->paths));
		eids[1] = get_u32v(revs, 0);
		e = ref_elist(g->edges, eids[1]);
		if(e->vid2 == vid) reverse_edge(g, eids[1]);
		encap_bitvec(g->flags[1], count_plist(e->paths));
		zeros_bitvec(g->flags[0]);
		zeros_bitvec(g->flags[1]);
		if(connect_edges(g, eids, g->flags, 2, &eid)){
			remove_edge(g, eids[0]);
			remove_edge(g, eids[1]);
			return_vertex(g, vid);
			ret = 3;
		} else ret = 1;
	}
	return ret;
}

void simplify_seqgraph(SeqGraph *g){
	u32v *fwds, *revs;
	uint32_t vid1;
	uint64_t ret, cad;
	ret = 0;
	cad = 0;
	fwds = init_u32v(4);
	revs = init_u32v(4);
	for(vid1=0;vid1<count_vlist(g->verts);vid1++){
		if(get_bitvec(g->v_flags, vid1) == 0) continue;
		switch(simplify_vertex_core(g, vid1, fwds, revs)){
			case 3: ret ++;
			case 1: cad ++; break;
		}
	}
	free_u32v(fwds);
	free_u32v(revs);
	fprintf(stdout, " -- merge %u/%u cont edges in %s -- %s:%d --\n", (unsigned)ret, (unsigned)cad, __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
}

uint32_t remove_tips_seqgraph(SeqGraph *g){
	Vertex *v;
	Edge *e;
	ReadPath *p, *pp;
	ReadTracker RT, *rt;
	uint32_t vid, vdir, nvid, nvdir, eid, edir, rdir, i, rlen, nc, ret;
	ret = 0;
	rt = &RT;
	for(vid=0;vid<count_vlist(g->verts);vid++){
		if(get_bitvec(g->v_flags, vid) == 0) continue;
		v = ref_vlist(g->verts, vid);
		if(count_u32v(v) == 0){ return_vertex(g, vid); continue; }
		if(count_u32v(v) > 1) continue;
		eid = get_u32v(v, 0);
		e = ref_elist(g->edges, eid);
		if(e->len >= g->max_rd_len) continue;
		for(i=0;i<count_plist(e->paths);i++){
			p = ref_plist(e->paths, i);
			rlen = read_length(g, p->rid);
			if(p->roff == 0 && p->rlen == rlen) break;
		}
		if(i < count_plist(e->paths)) continue;
		if(e->vid1 == e->vid2) continue;
		ret ++;
		if(e->vid1 == vid){
			vdir  = e->vdir1;
			nvid  = e->vid2;
			nvdir = e->vdir2;
			edir  = 0;
		} else {
			vdir  = e->vdir2;
			nvid  = e->vid1;
			nvdir = e->vdir1;
			edir  = 1;
		}
		for(i=0;i<count_plist(e->paths);i++){
			p = ref_plist(e->paths, i);
			rdir = p->rdir ^ edir;
			rlen = read_length(g, p->rid);
			set_readtracker(rt, g, p->rid, rdir, rlen, p->rlen, nvid, !nvdir);
			if(next_readtracker(rt) == 0){
				fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
				fflush(stdout);
				//abort();
			}
			pp = ref_plist(ref_elist(g->edges, rt->eid)->paths, rt->pid);
#ifdef DEBUG
			fprintf(stdout, " -- edir %c rid = %u %c %u+%u vid %u->%u pp->len %u in %s -- %s:%d --\n", "+-"[edir], p->rid, "+-"[p->rdir], p->roff, p->rlen, vid, rt->vid, pp->rlen, __FUNCTION__, __FILE__, __LINE__);
			fflush(stdout);
#endif
			if(rdir){
				if(pp->rdir){
					nc = append_cigars(g->aux_cigars[0], 0, ALN_CIGAR_TYPE_CLIP1, p->rlen);
					nc = cat_cigars(g->aux_cigars[0], nc, pp->cigars, pp->n_cigar);
				} else {
					nc = cat_cigars(g->aux_cigars[0], 0, pp->cigars, pp->n_cigar);
					nc = append_cigars(g->aux_cigars[0], nc, ALN_CIGAR_TYPE_CLIP1, p->rlen);
				}
				pp->rlen += p->rlen; change_r2v_map(g, p->rid, rdir, rt->vid, rt->vdir);
				if(pp->n_cigar < nc) pp->cigars = realloc(pp->cigars, sizeof(AlnCigar) * nc);
				pp->n_cigar = cat_cigars(pp->cigars, 0, g->aux_cigars[0], nc);
			} else {
				if(pp->rdir){
					nc = cat_cigars(g->aux_cigars[0], 0, pp->cigars, pp->n_cigar);
					nc = append_cigars(g->aux_cigars[0], nc, ALN_CIGAR_TYPE_CLIP1, p->rlen);
				} else {
					nc = append_cigars(g->aux_cigars[0], 0, ALN_CIGAR_TYPE_CLIP1, p->rlen);
					nc = cat_cigars(g->aux_cigars[0], nc, pp->cigars, pp->n_cigar);
				}
				pp->roff = 0; pp->rlen += p->rlen; change_r2v_map(g, p->rid, rdir, rt->vid, rt->vdir);
				if(pp->n_cigar < nc) pp->cigars = realloc(pp->cigars, sizeof(AlnCigar) * nc);
				pp->n_cigar = cat_cigars(pp->cigars, 0, g->aux_cigars[0], nc);
			}
		}
		remove_edge(g, eid);
		return_vertex(g, vid);
	}
	fprintf(stdout, " -- remove %u tips in %s -- %s:%d --\n", ret, __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
	return ret;
}

uint32_t resolve_small_repeats_seqgraph(SeqGraph *g){
	Edge *e;
	ReadPath *p;
	ReadTracker RT;
	u32v *eids;
	uuhash *fwds, *revs;
	uuhash_t *uu, UU, *uu2;
	uint32_t eid, i, eid1, eid2, n, rdlen, ret;
	int exists;
	eids = init_u32v(4);
	fwds = init_uuhash(13);
	revs = init_uuhash(13);
	ret = 0;
	for(eid=0;eid<count_elist(g->edges);eid++){
		if(get_bitvec(g->e_flags, eid) == 0) continue;
		e = ref_elist(g->edges, eid);
		if(e->len >= g->max_rd_len) continue;
		if(e->vid1 == e->vid2) continue;
		clear_u32v(eids);
		if(find_edges_v2d(g, e->vid1, e->vdir1, eids) > 1) continue;
		if(find_edges_v2d(g, e->vid2, e->vdir2, eids) > 1) continue;
		if(find_edges_v2d(g, e->vid1, !e->vdir1, eids) < 2) continue;
		if(find_edges_v2d(g, e->vid2, !e->vdir2, eids) < 2) continue;
		clear_uuhash(fwds);
		clear_uuhash(revs);
		for(i=0;i<count_plist(e->paths);i++){
			p = ref_plist(e->paths, i);
			rdlen = read_length(g, p->rid);
			if(p->roff == 0 || (uint32_t)(p->roff + p->rlen) == rdlen) continue;
			set_readtracker(&RT, g, p->rid, 0, rdlen, p->roff + p->rlen, p->rdir? e->vid1 : e->vid2, !(p->rdir? e->vdir1 : e->vdir2));
			if(!next_readtracker(&RT)){
				fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
				fflush(stdout);
				abort();
			}
			eid1 = RT.eid;
			set_readtracker(&RT, g, p->rid, 1, rdlen, rdlen - p->roff, p->rdir? e->vid2 : e->vid1, !(p->rdir? e->vdir2 : e->vdir1));
			if(!next_readtracker(&RT)){
				fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
				fflush(stdout);
				abort();
			}
			eid2 = RT.eid;
			UU.key = eid1;
			UU.val = eid2;
			uu = prepare_uuhash(p->rdir? revs : fwds, UU, &exists);
			if(exists){
				if(uu->val != UU.val) uu->val = 0xFFFFFFFFU;
			} else {
				uu->key = eid1;
				uu->val = eid2;
			}
			UU.key = eid2;
			UU.val = eid1;
			uu = prepare_uuhash(p->rdir? fwds : revs, UU, &exists);
			if(exists){
				if(uu->val != UU.val) uu->val = 0xFFFFFFFFU;
			} else {
				uu->key = eid2;
				uu->val = eid1;
			}
		}
		reset_iter_uuhash(fwds);
		encap_bitvec(g->flags[1], count_plist(e->paths));
		zeros_bitvec(g->flags[1]);
		while((uu = ref_iter_uuhash(fwds))){
			if(uu->val >= uu->key) continue;
			UU.key = uu->val;
			UU.val = uu->key;
			uu2 = get_uuhash(revs, UU);
			if(uu2->val != UU.val) continue;
			if(connect_e2e2e(g, uu->val, eid, uu->key)) ret ++;
		}
		n = 0;
		e = ref_elist(g->edges, eid);
		for(i=0;i<count_plist(e->paths);i++){ if(get_bitvec(g->flags[1], i) == 0) n ++; }
		if(n == 0) remove_edge(g, eid);
	}
	free_u32v(eids);
	free_uuhash(fwds);
	free_uuhash(revs);
	fprintf(stdout, " -- resolve %u small repeats in %s -- %s:%d --\n", (unsigned)ret, __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
	return ret;
}

//TODO -- Program --

int main(int argc, char **argv){
	SeqGraph *g;
	FileReader *rds, *ols;
	FILE *out;
	char name[256];
	int fix_rd_len;
	if(argc < 3){
		fprintf(stdout, "Usage: seqgraph <rds_file> <ols_file> <is_fix_rd_len:0>\n");
		return 1;
	}
	if(argc > 3) fix_rd_len = atoi(argv[3]);
	else fix_rd_len = 0;
	if((rds = fopen_filereader(argv[1])) == NULL){
		fprintf(stdout, " -- Cannot open %s in %s -- %s:%d --\n", argv[1], __FUNCTION__, __FILE__, __LINE__);
		return 1;
	}
	if((ols = fopen_filereader(argv[2])) == NULL){
		fprintf(stdout, " -- Cannot open %s in %s -- %s:%d --\n", argv[2], __FUNCTION__, __FILE__, __LINE__);
		return 1;
	}
	g = init_seqgraph(rds, fix_rd_len);
	fclose_filereader(rds);
	overlapping_seqgraph(g, ols);
	fclose_filereader(ols);
	fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
	check_seqgraph_edges(g);
	fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
	check_seqgraph_link_refs(g);
	fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
	check_seqgraph_readpaths(g);
	fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
	simplify_seqgraph(g);
	resolve_small_repeats_seqgraph(g);
	remove_tips_seqgraph(g);
	fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
	check_seqgraph_edges(g);
	fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
	check_seqgraph_link_refs(g);
	fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
	check_seqgraph_readpaths(g);
	fprintf(stdout, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
	fflush(stdout);
	sprintf(name, "%s.pre.vertex", argv[1]);
	out = fopen(name, "w");
	output_vertexs_seqgraph(g, out);
	fclose(out);
	sprintf(name, "%s.pre.edge", argv[1]);
	out = fopen(name, "w");
	output_edges_seqgraph(g, out);
	fclose(out);
	sprintf(name, "%s.pre.dot", argv[1]);
	out = fopen(name, "w");
	output_dot_seqgraph(g, out);
	fclose(out);
	free_seqgraph(g);
	return 0;
}
