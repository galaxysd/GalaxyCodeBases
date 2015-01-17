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
 
#include "ps_assembly.h"

int _cmp_min_distance_nodes(uint32_t idx1, uint32_t idx2, void *ref){
	Graph *g;
	trace_t *t1, *t2;
	edge_t *e1, *e2;
	g = (Graph*)ref;
	t1 = ref_tracev(g->tree, idx1);
	t2 = ref_tracev(g->tree, idx2);
	e1 = ref_edgev(ref_nodev(g->nodes, t1->bt)->edges, t1->bt_idx);
	e2 = ref_edgev(ref_nodev(g->nodes, t2->bt)->edges, t2->bt_idx);
	cmp_2nums_proc(e1->off, e2->off);
	cmp_2nums_proc(e1->mm, e2->mm);
	cmp_2nums_proc(t1->off, t2->off);
	return 0;
}


Graph* init_graph(uint32_t ksize, uint32_t min_ol, float min_sm){
	Graph *g;
	g = malloc(sizeof(Graph));
	g->ps_id = 0;
	g->ps_len = 0;
	g->rd_len = 0;
	g->min_ins = 0;
	g->max_ins = 0;
	g->n_sr = 0;
	g->n_l = 0;
	g->n_r = 0;
	g->rdseqs = init_basebank();
	g->nodes = init_nodev(64);
	g->m_node = 0;
	if(ksize > 16) ksize = 16;
	if(ksize < 5) ksize = 5;
	g->ksize = ksize;
	g->kmask = 0xFFFFFFFFU >> ((16 - ksize) << 1);
	g->min_ol = min_ol;
	g->min_sm = min_sm;
	g->out_mode = 1;
	g->hash  = NULL;
	g->links = NULL;
	g->conns = init_bitvec(1024);
	g->tested = init_u8list(1024);
	g->tree  = init_tracev(1024);
	g->stack = init_u32list(1024);
	g->heap  = init_heap(_cmp_min_distance_nodes, g);
	g->fast  = init_tracev(64);
	g->slow  = init_tracev(64);
	g->fast_seq = init_string(1024);
	g->slow_seq = init_string(1024);
	g->aux_hash = init_u32hash(13);
	g->loop_hash = init_uuhash(13);
	g->cnss  = init_cnsv(6);
	g->m_cns = 0;
	g->aux   = sr_init_aux();
	return g;
}

void free_graph(Graph *g){
	uint32_t i;
	cns_t *cns;
	free_basebank(g->rdseqs);
	for(i=0;i<g->m_node;i++){ free_edgev(ref_nodev(g->nodes, i)->edges); }
	free_nodev(g->nodes);
	free_bitvec(g->conns);
	free_u8list(g->tested);
	if(g->hash){ free_rdkhash(g->hash); g->hash = NULL; }
	if(g->links){ free_rdoffv(g->links); g->links = NULL; }
	for(i=0;i<g->m_cns;i++){
		cns = ref_cnsv(g->cnss, i);
		free_tracev(cns->path);
		free_cbv(cns->cbs);
		free_string(cns->seq);
		free_basebank(cns->qc_seqs);
		free_qc_alnv(cns->alns);
	}
	free_cnsv(g->cnss);
	free_tracev(g->tree);
	free_u32list(g->stack);
	free_heap(g->heap);
	free_tracev(g->fast);
	free_tracev(g->slow);
	free_string(g->fast_seq);
	free_string(g->slow_seq);
	free_u32hash(g->aux_hash);
	free_uuhash(g->loop_hash);
	sr_free_aux(g->aux);
	free(g);
}

void reset_graph(Graph *g){
	clear_basebank(g->rdseqs);
	clear_nodev(g->nodes);
	clear_cnsv(g->cnss);
	zeros_bitvec(g->conns);
	//clear_u8list(g->tested);
	clear_uuhash(g->loop_hash);
	if(g->hash){ free_rdkhash(g->hash); g->hash = NULL; }
	if(g->links){ free_rdoffv(g->links); g->links = NULL; }
	g->rd_len = 0;
	g->n_l = 0;
	g->n_r = 0;
	g->n_sr = 0;
}

void set_insert_graph(Graph *g, uint32_t min, uint32_t max){
	g->min_ins = min;
	g->max_ins = max;
}

void push_graph(Graph *g, uint32_t rid, char *rdseq, uint8_t rdlen, uint32_t rddir, uint32_t anchor_type){
	node_t *node;
	node = next_ref_nodev(g->nodes);
	if(g->nodes->size > g->m_node){
		node->edges = init_edgev(12);
		g->m_node ++;
	} else {
		clear_edgev(node->edges);
	}
	node->rid = rid;
	node->seqdir = rddir;
	node->closed = 0;
	node->bt_nid = 0;
	node->bt_eid = 0;
	node->bt_off = 0;
	node->visit  = 0;
	node->sort   = 0;
	node->n_edge = 0;
	if(g->rd_len){
		if(rdlen != g->rd_len){
			fprintf(stderr, " -- Error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
		}
	} else {
		g->rd_len = rdlen;
	}
	if(anchor_type & 0x02) g->n_r ++;
	else if(anchor_type & 0x01) g->n_l ++;
	if(anchor_type & 0x04) node->closed = 1;
	seq2basebank(g->rdseqs, rdseq, rdlen);
}

void ready_graph(Graph *g){
	if(g->nodes->size < 2) return;
	g->n_sr = g->nodes->size - 2;
	g->ps_id = ref_nodev(g->nodes, 0)->rid >>1;
}

void index_graph(Graph *g){
	rdoff_t *rdf;
	rdkmer_t *k, K;
	uint64_t offset;
	uint32_t i, j, cnt;
	int exists;
	K.kmer = 0;
	K.off  = 0;
	K.cnt  = 0;
	g->hash = init_rdkhash(4481);
	for(i=0;i<g->nodes->size;i++){
		offset = ((uint64_t)i) * g->rd_len;
		K.kmer = 0;
		for(j=0;j<g->rd_len;j++){
			K.kmer = ((K.kmer << 2) | (bits2bit(g->rdseqs->bits, offset + j)) ) & g->kmask;
			if(j + 1 < g->ksize) continue;
			k = prepare_rdkhash(g->hash, K, &exists);
			if(exists) k->cnt ++;
			else { k->kmer = K.kmer; k->off = 0; k->cnt = 1; }
		}
	}
	cnt = 0;
	reset_iter_rdkhash(g->hash);
	while((k = ref_iter_rdkhash(g->hash))){
		k->off = cnt;
		cnt += k->cnt;
		k->cnt = 0;
	}
	g->links = init_rdoffv(cnt);
	for(i=0;i<g->nodes->size;i++){
		offset = ((uint64_t)i) * g->rd_len;
		K.kmer = 0;
		for(j=0;j<g->rd_len;j++){
			K.kmer = ((K.kmer << 2) | (bits2bit(g->rdseqs->bits, offset + j)) ) & g->kmask;
			if(j + 1 < g->ksize) continue;
			k = get_rdkhash(g->hash, K);
			rdf = ref_rdoffv(g->links, k->off + k->cnt);
			rdf->nid = i;
			rdf->off = j + 1 - g->ksize;
			k->cnt ++;
		}
	}
}

void _remove_closed_edges(node_t *n){
	edge_t *e1, *e2, e3;
	uint32_t j, k;
	k = 0;
	for(j=0;j<n->edges->size;j++){
		e1 = ref_edgev(n->edges, j);
		if(e1->closed) continue;
		if(k < j){
			e2 = ref_edgev(n->edges, k);
			e3 = *e2;
			*e2 = *e1;
			*e1 = e3;
		}
		k ++;
	}
	n->edges->size = k;
}

void align_graph(Graph *g){
	rdkmer_t *k;
	rdoff_t *f1, *f2;
	node_t *n1, *n2;
	edge_t *e;
	uint32_t i, j, mm, ol;
	encap_bitvec(g->conns, g->nodes->size * g->nodes->size);
	encap_u8list(g->tested, g->nodes->size * g->nodes->size);
	zeros_u8list(g->tested);
	reset_iter_rdkhash(g->hash);
	while((k = ref_iter_rdkhash(g->hash))){
		sort_array(g->links->buffer + k->off, k->cnt, rdoff_t, (a.off > b.off));
		for(i=0;i<k->cnt;i++){
			f1 = ref_rdoffv(g->links, k->off + i);
			n1 = ref_nodev(g->nodes, f1->nid);
			if(n1->closed) continue;
			for(j=i+1;j<k->cnt;j++){
				f2 = ref_rdoffv(g->links, k->off + j);
				if(f1->nid == f2->nid) continue;
				n2 = ref_nodev(g->nodes, f2->nid);
				if(n2->closed) continue;
				ol = g->rd_len - ((f1->off <= f2->off)? (f2->off - f1->off): (f1->off - f2->off));
				if(ol < g->min_ol) continue;
				if(get_u8list(g->tested, f1->nid * g->nodes->size + f2->nid) >= ol) continue;
				set_u8list(g->tested, f1->nid * g->nodes->size + f2->nid, ol);
				set_u8list(g->tested, f2->nid * g->nodes->size + f1->nid, ol);
				if(f1->off == f2->off){
					mm = mismatch_basebank(g->rdseqs, f1->nid * g->rd_len, f2->nid * g->rd_len, ol);
					if(mm > (uint32_t)(ol * (1 - g->min_sm) + 0.5)) continue;
					if(f1->nid < 2){
						if(f2->nid >= 2) n2->closed = 1;
					} else if(f2->nid < 2){
						if(f1->nid >= 2){ n1->closed = 1; break; }
					} else if(f1->nid < f2->nid){
						n2->closed = 1;
					} else {
						n1->closed = 1; break;
					}
				} else if(f1->off < f2->off){
					//if(f2->nid == 1 || f1->nid == 0) continue;
					mm = mismatch_basebank(g->rdseqs, f1->nid * g->rd_len, f2->nid * g->rd_len + f2->off - f1->off, ol);
					if(mm > (uint32_t)(ol * (1 - g->min_sm) + 0.5)) continue;
					e = next_ref_edgev(n2->edges);
					e->nid = f1->nid;
					e->off = g->rd_len - ol;
					e->mm  = mm;
					e->closed = 0;
					e->flag   = 0;
					e->loop   = 0;
				} else {
					//if(f1->nid == 1 || f2->nid == 0) continue;
					mm = mismatch_basebank(g->rdseqs, f1->nid * g->rd_len + f1->off - f2->off, f2->nid * g->rd_len, ol);
					if(mm > (uint32_t)(ol * (1 - g->min_sm) + 0.5)) continue;
					e = next_ref_edgev(n1->edges);
					e->nid = f2->nid;
					e->off = g->rd_len - ol;
					e->mm  = mm;
					e->closed = 0;
					e->flag   = 0;
					e->loop   = 0;
				}
			}
		}
	}
	for(i=0;i<g->nodes->size;i++){
		n1 = ref_nodev(g->nodes, i);
		for(j=0;j<n1->edges->size;j++){
			e = ref_edgev(n1->edges, j);
			n2 = ref_nodev(g->nodes, e->nid);
			if(n2->closed) e->closed = 1;
		}
		_remove_closed_edges(n1);
	}
	for(i=0;i<g->nodes->size;i++){
		n1 = ref_nodev(g->nodes, i);
		n1->n_edge = n1->edges->size;
		if(!n1->sort){
			sort_array(n1->edges->buffer, n1->edges->size, edge_t, (a.off > b.off));
			n1->sort = 1;
		}
		for(j=0;j<n1->edges->size;j++){
			e = ref_edgev(n1->edges, j);
			if(e->closed) continue;
			one_bitvec(g->conns, i * g->nodes->size + e->nid);
			one_bitvec(g->conns, e->nid * g->nodes->size + i);
		}
	}
}

void simplify_graph(Graph *g){
	node_t *n;
	edge_t *e1, *e2;
	uint32_t i, j, k;
#if 0
	node_t *nn;
	edge_t *e3, *e4;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->edges->size < 2) continue;
		sort_array(n->edges->buffer, n->edges->size, edge_t, (a.off > b.off));
		n->sort = 1;
		for(j=0;j+1<n->edges->size;j++){
			e1 = ref_edgev(n->edges, j);
			e2 = ref_edgev(n->edges, j + 1);
			if(get_bitvec(g->conns, e1->nid * g->nodes->size + e2->nid)) continue;
			for(k=j+2;k<n->edges->size;k++){
				e3 = ref_edgev(n->edges, k);
				if(get_bitvec(g->conns, e1->nid * g->nodes->size + e3->nid) == 0) continue;
				if(get_bitvec(g->conns, e2->nid * g->nodes->size + e3->nid) == 0) continue;
				one_bitvec(g->conns, e1->nid * g->nodes->size + e2->nid);
				one_bitvec(g->conns, e2->nid * g->nodes->size + e1->nid);
				nn = ref_nodev(g->nodes, e1->nid);
				nn->sort = 0;
				e4 = next_ref_edgev(nn->edges);
				e4->nid = e2->nid;
				e4->off = e2->off - e1->off;
				e4->closed = 0;
				e4->flag = 0;
				break;
			}
		}
	}
#endif
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		n->n_edge = n->edges->size;
		if(n->edges->size < 2) continue;
		if(n->sort == 0){
			sort_array(n->edges->buffer, n->edges->size, edge_t, (a.off > b.off));
			n->sort = 1;
		}
		for(j=n->edges->size-1;j>0;j--){
			e1 = ref_edgev(n->edges, j);
			if(e1->closed) continue;
			for(k=0;k<j;k++){
				e2 = ref_edgev(n->edges, k);
				if(e2->closed) continue;
				if(get_bitvec(g->conns, e1->nid * g->nodes->size + e2->nid)) break;
			}
			if(k < j){ e1->closed = 1; }
		}
		_remove_closed_edges(n);
	}
}

void _traceback_graph(Graph *g, uint32_t nid, uint32_t fast_bt, uint32_t slow_bt, uint32_t fast_eid, uint32_t slow_eid, uint32_t fast_off, uint32_t slow_off){
	node_t *n1, *n2;
	edge_t *e1, *e2;
	trace_t *p1, *p2;
	uint32_t nid1, nid2, *ref;
	int exists, off1, off2;
	g->fast_break = 0;
	g->slow_break = 0;
	clear_tracev(g->fast);
	clear_tracev(g->slow);
	clear_u32hash(g->aux_hash);
	p1 = next_ref_tracev(g->fast);
	p1->nid = nid;
	p1->bt  = fast_bt;
	p1->bt_idx = fast_eid;
	p1->off = fast_off;
	p2 = next_ref_tracev(g->slow);
	p2->nid = nid;
	p2->bt  = slow_bt;
	p2->bt_idx = slow_eid;
	p2->off = slow_off;
	put_u32hash(g->aux_hash, nid);
	while(1){
		p1 = peer_tracev(g->fast);
		off1 = p1->off;
		p2 = peer_tracev(g->slow);
		off2 = p2->off;
		if(off1 <= off2){
			nid2 = p2->bt;
			n2 = ref_nodev(g->nodes, nid2);
			e2 = ref_edgev(n2->edges, p2->bt_idx);
			if(e2->flag) g->slow_break = 1;
			p2 = next_ref_tracev(g->slow);
			p2->nid = nid2;
			p2->bt  = n2->bt_nid;
			p2->bt_idx = n2->bt_eid;
			p2->off = off2 - e2->off;
			ref = prepare_u32hash(g->aux_hash, nid2, &exists);
			if(exists){
				break;
			} else *ref = nid2;
		} else {
			nid1 = p1->bt;
			n1 = ref_nodev(g->nodes, nid1);
			e1 = ref_edgev(n1->edges, p1->bt_idx);
			if(e1->flag) g->fast_break = 1;
			p1 = next_ref_tracev(g->fast);
			p1->nid = nid1;
			p1->bt  = n1->bt_nid;
			p1->bt_idx = n1->bt_eid;
			p1->off = off1 - e1->off;
			ref = prepare_u32hash(g->aux_hash, nid1, &exists);
			if(exists){
				break;
			} else *ref = nid1;
		}
	}
	p1 = peer_tracev(g->fast);
	p2 = peer_tracev(g->fast);
	if(p1->nid != p2->nid){
		fprintf(stderr, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
	}
}

void _traces2seq(Graph *g, tracev *path, String *seq){
	trace_t *p;
	uint32_t i, length, off;
	int offset;
	clear_string(seq);
	offset = 0;
	length = 0;
	for(i=0;i<path->size;i++){
		p = ref_tracev(path, i);
		if(i == 0){ offset = p->off; }
		off = length + offset - p->off;
		encap_string(seq, length + g->rd_len - off);
		seq_basebank(g->rdseqs, p->nid * g->rd_len + off, g->rd_len - off, seq->string + length);
		length += g->rd_len - off;
	}
	seq->size = length;
}

int _is_similar_seqs(String *seq1, String* seq2, float min_sm){
	AlnAln *aln;
	path_t  *aln_p;
	uint32_t len, mm;
	aln = aln_stdaln(seq1->string, seq2->string, &aln_param_nt2nt, 1, 0);
	len = mm = 0;
	aln_p = aln->path + aln->path_len - 1;
	while(aln_p >= aln->path){
		if(aln_p->ctype == FROM_M){
			len ++;
			if(seq1->string[aln_p->i] != seq2->string[aln_p->j]) mm ++;
		}
		aln_p --;
	}
	aln_free_AlnAln(aln);
	return (mm <= (uint32_t)(min_sm * len + 0.5));
}

uint32_t _connective_score(tracev *path){
	uint32_t t, s;
	uint32_t i, len;
	t = s = 0;
	for(i=0;i+1<path->size;i++){
		len = ref_tracev(path, i)->off - ref_tracev(path, i + 1)->off;
		t += len * len;
		s ++;
	}
	return t / s;
}

uint32_t shave_graph(Graph *g){
	node_t *n;
	edge_t *e;
	trace_t *p, *p1, *p2;
	uint32_t i, j, pidx, nidx, len1, len2, min_len, score1, score2, ret;
	ret = 0;
	clear_tracev(g->tree);
	clear_heap(g->heap);
	push_heap(g->heap, 0);
	p = next_ref_tracev(g->tree);
	p->nid = 0;
	p->bt  = 0;
	p->bt_idx = 0;
	p->off = 0;
	while((pidx = pop_heap(g->heap)) != 0xFFFFFFFFU){
		p = ref_tracev(g->tree, pidx);
		nidx = p->nid;
		n = ref_nodev(g->nodes, nidx);
		if(n->visit){
			_traceback_graph(g, p->nid, n->bt_nid, p->bt, n->bt_eid, p->bt_idx, n->bt_off, p->off);
			// Loop
			if(g->fast->size == 1){
#ifdef ENABLE_LOOP
				if(peer_tracev(g->fast)->nid == 0 || g->loop_hash->count >= MAX_LOOP){
#endif
					ref_edgev(ref_nodev(g->nodes, p->bt)->edges, p->bt_idx)->flag = 1;
					ret ++;
#ifdef ENABLE_LOOP
				} else {
					ref_edgev(ref_nodev(g->nodes, p->bt)->edges, p->bt_idx)->loop = 1;
					kv_put_uuhash(g->loop_hash, (p->bt << 16) | p->bt_idx, g->loop_hash->count);
				}
#endif
				continue;
			}
			if(g->slow->size == 1){
#ifdef ENABLE_LOOP
				if(peer_tracev(g->slow)->nid == 0 || g->loop_hash->count >= MAX_LOOP){
#endif
					ref_edgev(ref_nodev(g->nodes, n->bt_nid)->edges, n->bt_eid)->flag = 1;
					ret ++;
#ifdef ENABLE_LOOP
				} else {
					ref_edgev(ref_nodev(g->nodes, n->bt_nid)->edges, n->bt_eid)->loop = 1;
					kv_put_uuhash(g->loop_hash, (n->bt_nid << 16) | n->bt_eid, g->loop_hash->count);
				}
#endif
				continue;
			}
			//if(g->fast_break || g->slow_break) continue;
			len1 = ref_tracev(g->fast, 0)->off - peer_tracev(g->fast)->off;
			score1 = _connective_score(g->fast);
			len2 = ref_tracev(g->slow, 0)->off - peer_tracev(g->slow)->off;
			score2 = _connective_score(g->slow);
			min_len = (len1 < len2)? len1 : len2;
			if(min_len <= g->rd_len){
				if(score1 > score2){
					p2 = ref_tracev(g->fast, 0);
					n->bt_nid = p->bt;
					n->bt_eid = p->bt_idx;
					n->bt_off = p->off;
				} else {
					p2 = ref_tracev(g->slow, 0);
				}
				ref_edgev(ref_nodev(g->nodes, p2->bt)->edges, p2->bt_idx)->flag = 1;
				ret ++;
			} else {
				reverse_tracev(g->fast);
				_traces2seq(g, g->fast, g->fast_seq);
				reverse_tracev(g->slow);
				_traces2seq(g, g->slow, g->slow_seq);
				if(_is_similar_seqs(g->fast_seq, g->slow_seq, g->min_sm)){
					if(score1 > score2){
						p2 = peer_tracev(g->fast);
						n->bt_nid = p->bt;
						n->bt_eid = p->bt_idx;
						n->bt_off = p->off;
					} else {
						p2 = peer_tracev(g->slow);
					}
					ref_edgev(ref_nodev(g->nodes, p2->bt)->edges, p2->bt_idx)->flag = 1;
					ret ++;
				}
			}
			continue;
		}
		n->visit = 1;
		n->bt_nid = p->bt;
		n->bt_eid = p->bt_idx;
		n->bt_off = p->off;
		if(p->nid == 1) continue;
		if(g->tree->size + n->edges->size >= MAX_TRACE) break;
		for(i=0;i<n->edges->size;i++){
			e = ref_edgev(n->edges, i);
			if(e->closed) continue;
			p1 = next_ref_tracev(g->tree);
			p1->nid = e->nid;
			p1->bt  = nidx;
			p1->bt_idx = i;
			p1->off = n->bt_off + e->off;
			push_heap(g->heap, g->tree->size - 1);
		}
	}
	if(ret == 0) return 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		for(j=0;j<n->edges->size;j++){
			e = ref_edgev(n->edges, j);
			if(e->flag) e->closed = 1;
		}
		_remove_closed_edges(n);
	}
	return ret;
}

void print_dot_graph(Graph *g, FILE *out){
	node_t *n;
	edge_t *e;
	uint32_t i, j;
	fprintf(out, "digraph G%u {\n", g->ps_id);
	fprintf(out, "rankdir=LR\n");
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		if(i != 1 && !n->visit) continue;
		if(i < 2){
			fprintf(out, "N%03u [label=\"N%03u(%u)\" shape=rect color=red]\n", i, i, n->bt_off);
		} else {
			fprintf(out, "N%03u [label=\"N%03u(%u)\"]\n", i, i, n->bt_off);
		}
	}
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		if(!n->visit) continue;
		for(j=0;j<n->n_edge;j++){
			e = ref_edgev(n->edges, j);
			//fprintf(out, "N%03u -> N%03u [label=\"%u\"color=blue]\n", i, e->nid, e->off);
			if(e->flag){
				fprintf(out, "N%03u -> N%03u [label=\"%u\"color=yellow]\n", i, e->nid, e->off);
			} else if(!e->closed){
				if(e->loop){
					fprintf(out, "N%03u -> N%03u [label=\"%u\"color=green]\n", i, e->nid, e->off);
				} else {
					fprintf(out, "N%03u -> N%03u [label=\"%u\"color=blue]\n", i, e->nid, e->off);
				}
			}
		}
	}
	fprintf(out, "}\n");
	fflush(out);
}

void print_dot_graph_raw(Graph *g, FILE *out){
	node_t *n;
	edge_t *e;
	uint32_t i, j;
	fprintf(out, "digraph G%u {\n", g->ps_id);
	fprintf(out, "rankdir=LR\n");
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		for(j=0;j<n->edges->size;j++){
			e = ref_edgev(n->edges, j);
			if(e->closed) continue;
			fprintf(out, "N%03u -> N%03u [label=\"%u\"color=blue]\n", i, e->nid, e->off);
		}
	}
	fprintf(out, "}\n");
	fflush(out);
}

uint32_t allpaths_graph(Graph *g){
	tracev *tree;
	u32list *stack;
	trace_t *p, *pp;
	node_t *n;
	edge_t *e;
	u32hash *hash;
	cns_t *cns;
	uint32_t i, pidx, nid, poff, ret;
#ifdef ENABLE_LOOP
	uint32_t loop;
#endif
	int uniq;
	ret = 0;
	tree = g->tree;
	clear_tracev(tree);
	stack = g->stack;
	clear_u32list(stack);
	hash = g->aux_hash;
	pp = next_ref_tracev(tree);
	pp->nid = 0;
	pp->bt  = MAX_TRACE;
	pp->bt_idx = 0;
	pp->off = 0;
#ifdef ENABLE_LOOP
	pp->loop = 0;
#endif
	// Enable uniq mode
	pp->uniq = (ref_nodev(g->nodes, 0)->edges->size == 1);
	// Disable uniq mode
	//pp->uniq = 0;
	push_u32list(stack, 0);
	// May cause loops, e->off > 0, so loops will run out of max_ins
	while(pop_u32list(stack, &pidx)){
		p = ref_tracev(tree, pidx);
		poff = p->off;
		nid  = p->nid;
		uniq = p->uniq;
#ifdef ENABLE_LOOP
		loop = p->loop;
#endif
		p = NULL;
		n = ref_nodev(g->nodes, nid);
		if(n->edges->size != 1) uniq = 0;
		for(i=0;i<n->edges->size;i++){
			e = ref_edgev(n->edges, i);
			//if(e->closed) continue;
			if(uniq){
				if(poff + e->off + g->rd_len > (uint32_t)(1.25 * g->max_ins)) continue;
			} else if(poff + e->off + g->rd_len > g->max_ins) continue;
			if(e->nid == 0) continue;
			if(e->nid == 1){
				if(uniq){
					if(poff + e->off + g->rd_len < (uint32_t)(0.75 * g->min_ins)) continue;
				} else if(poff + e->off + g->rd_len < g->min_ins) continue;
				cns = next_ref_cnsv(g->cnss);
				if(g->cnss->size > g->m_cns){
					cns->path = init_tracev(256);
					cns->cbs  = init_cbv(1024);
					cns->seq  = init_string(1024);
					cns->qc_seqs = init_basebank();
					cns->alns    = init_qc_alnv(256);
					g->m_cns ++;
				} else {
					clear_tracev(cns->path);
					clear_cbv(cns->cbs);
					clear_string(cns->seq);
					clear_basebank(cns->qc_seqs);
					clear_qc_alnv(cns->alns);
				}
				cns->len = poff + e->off + g->rd_len;
				cns->qc  = 0;
				clear_u32hash(hash);
				pp = next_ref_tracev(cns->path);
				pp->nid = 1;
				pp->bt  = pidx;
				pp->off = poff + e->off;
				pp->bt_idx = i;
#ifdef ENABLE_LOOP
				pp->loop = loop;
#endif
				pp = ref_tracev(tree, pidx);
				push_tracev(cns->path, *pp);
				put_u32hash(hash, pidx);
				while(pp->bt != MAX_TRACE){
					pp = ref_tracev(tree, pp->bt);
					push_tracev(cns->path, *pp);
				}
				if(pp->nid != 0){
					fprintf(stderr, " -- in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
				}
				reverse_tracev(cns->path);
				ret ++;
				continue;
			}
#ifdef ENABLE_LOOP
			if(e->loop && ((loop >> kv_get_uuhash(g->loop_hash, (nid << 16) | i)) & 0x01)) continue;
#endif
			if(tree->size >= MAX_TRACE) break;
			push_u32list(stack, tree->size);
			pp = next_ref_tracev(tree);
			pp->nid = e->nid;
			pp->bt  = pidx;
			pp->bt_idx = i;
			pp->off = poff + e->off;
			pp->uniq = uniq;
#ifdef ENABLE_LOOP
			pp->loop = e->loop? (loop | (1U << kv_get_uuhash(g->loop_hash, (nid << 16) | i))) : loop;
#endif
		}
	}
	return ret;
}

int find_inter_repeat(cns_t *cns1, cns_t *cns2, u32hash *hash){
	trace_t *p;
	uint32_t i, lst, cnt;
	clear_u32hash(hash);
	for(i=1;i+1<cns1->path->size;i++){ put_u32hash(hash, ref_tracev(cns1->path, i)->nid); }
	lst = 1;
	cnt = 0;
	for(i=1;i+1<cns2->path->size;i++){
		p = ref_tracev(cns2->path, i);
		if(exists_u32hash(hash, p->nid)){
			lst = i + 1;
		} else {
			if(lst == i) cnt ++;
			lst = 0;
		}
	}
	if(cnt > 1) return 1;
	clear_u32hash(hash);
	for(i=1;i+1<cns2->path->size;i++){ put_u32hash(hash, ref_tracev(cns2->path, i)->nid); }
	lst = 1;
	cnt = 0;
	for(i=1;i+1<cns1->path->size;i++){
		p = ref_tracev(cns1->path, i);
		if(exists_u32hash(hash, p->nid)){
			lst = i + 1;
		} else {
			if(lst == i) cnt ++;
			lst = 0;
		}
	}
	if(cnt > 1) return 1;
	return 0;
}

inline uint32_t _my_abs(int a){ return (a < 0)? - a : a; }

uint32_t validate_paths_graph(Graph *g){
	cns_t *cns, *cns1, *cns2;
	u32hash *hash;
	uint32_t i, j;
	for(i=0;i<g->cnss->size;i++){
		cns = ref_cnsv(g->cnss, i);
		cns->len = peer_tracev(cns->path)->off  + g->rd_len - ref_tracev(cns->path, 0)->off;
	}
	if(g->cnss->size > 1){
		for(i=0;i<g->cnss->size;i++){
			cns = ref_cnsv(g->cnss, i);
			_traces2seq(g, cns->path, cns->seq);
		}
		for(i=0;i+1<g->cnss->size;i++){
			cns1 = ref_cnsv(g->cnss, i);
			if(cns1->len == 0) continue;
			for(j=i+1;j<g->cnss->size;j++){
				cns2 = ref_cnsv(g->cnss, j);
				if(cns2->len == 0) continue;
				if(_is_similar_seqs(cns1->seq, cns2->seq, g->min_sm)){
					if(_my_abs(cns1->seq->size - (g->min_ins + g->max_ins) / 2) > _my_abs(cns2->seq->size - (g->min_ins + g->max_ins) / 2)){
						cns1->len = 0;
					} else {
						cns2->len = 0;
					}
				}
			}
		}
		sort_array(g->cnss->buffer, g->cnss->size, cns_t, (b.len > a.len));
		while(g->cnss->size && peer_cnsv(g->cnss)->len == 0) g->cnss->size --;
	}
	if(g->cnss->size > 1){ hash = init_u32hash(13); }
	else hash = NULL;
	for(i=0;i+1<g->cnss->size;i++){
		cns1 = ref_cnsv(g->cnss, i);
		for(j=i+1;j<g->cnss->size;j++){
			cns2 = ref_cnsv(g->cnss, j);
			if(cns1->len == 0 && cns2->len == 0) continue;
			if(find_inter_repeat(cns1, cns2, hash)){
				cns1->len = cns2->len = 0;
			}
		}
	}
	if(hash) free_u32hash(hash);
	sort_array(g->cnss->buffer, g->cnss->size, cns_t, (a.len < b.len));
	while(g->cnss->size && peer_cnsv(g->cnss)->len == 0) g->cnss->size --;
	return g->cnss->size;
}

void _make_consensus(cns_t *cns){
	comb_base_t *b;
	uint32_t j, k, c1, c2;
	encap_string(cns->seq, cns->len);
	for(j=0;j<cns->len;j++){
		b = ref_cbv(cns->cbs, j);
		c1 = c2 = 0;
		for(k=1;k<4;k++){
			if(b->bases[k] > b->bases[c1]){ c2 = c1; c1 = k; }
			else if(b->bases[k] > b->bases[c2]){ c2 = k; }
		}
		// Try call the same smaller base in multiple contigs
		if(c2 < c1 && b->bases[c2] > 1 && b->bases[c2] >= (uint32_t)(b->bases[c1] * 0.66666 + 0.5)) c1 = c2;
		cns->seq->string[j] = bit_base_table[c1];
	}
	cns->seq->string[j] = 0;
	cns->seq->size = cns->len;
}

void consensus_graph(Graph *g){
	cns_t *cns;
	trace_t *p;
	comb_base_t *b;
	uint64_t off;
	uint32_t i, j, k;
	for(i=0;i<g->cnss->size;i++){
		cns = ref_cnsv(g->cnss, i);
		clear_cbv(cns->cbs);
		encap_cbv(cns->cbs, cns->len);
		memset(cns->cbs->buffer, 0, cns->len * sizeof(comb_base_t));
		clear_string(cns->seq);
		for(j=0;j<cns->path->size;j++){
			p = ref_tracev(cns->path, j);
			off = p->nid * g->rd_len;
			for(k=0;k<g->rd_len;k++){
				b = ref_cbv(cns->cbs, p->off + k);
				b->bases[bits2bit(g->rdseqs->bits, off + k)] ++;
			}
		}
		_make_consensus(cns);
	}
}

static inline uint32_t _mean_distance(u32list *list){
	uint64_t sum, avg, std;
	uint32_t i, c;
	sum = 0;
	std = 0;
	for(i=0;i<list->size;i++){
		sum += list->buffer[i];
		std += list->buffer[i] * list->buffer[i];
	}
	avg = sum / list->size;
	if(1) return avg;
	std = sqrt((std - avg * avg) / list->size);
	sum = 0;
	c = 0;
	for(i=0;i<list->size;i++){
		if(list->buffer[i] + 2 * std < avg) continue;
		if(list->buffer[i] > avg + 2 * std) continue;
		c ++;
		sum += list->buffer[i];
	}
	return c? (sum / c) : 0;
}

void realign_graph(Graph *g, SR_SeqDB *sdb, u32list *lib_cnts, u32list *lib_ins, int do_cns){
	cns_t *cns;
	trace_t *p;
	node_t *n;
	SR_AlnHit *h;
	qc_aln_t *q, *q1, *q2;
	comb_base_t *b;
	uint64_t off;
	uint32_t i, j, k, s, e, m, c, ins;
	if(g->rd_len != sdb->rd_len){
		fprintf(stderr, " -- Error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stderr); abort();
	}
	sr_set_aux_strand(g->aux, 3);
	sr_set_aux_overlap(g->aux, g->rd_len * 0.75, g->rd_len);
	sr_set_aux_min_similarity(g->aux, g->min_sm);
	sr_set_aux_gap(g->aux, 0);
	if(lib_cnts->size > 2){
		for(m=c=0;m+1<lib_cnts->size;c+=lib_cnts->buffer[m++]){
			if(lib_ins->buffer[m] < 0.45 * lib_ins->buffer[0]) break;
		}
		sr_set_aux_hit_id_region(g->aux, 0, c);
	} else {
		sr_set_aux_hit_id_region(g->aux, 0, 0xFFFFFFFFU);
	}
	sr_set_aux_max_hits(g->aux, 10000);
	sr_set_aux_max_hits_per_off(g->aux, 2);
	sr_set_aux_cigar(g->aux, 0);
	sr_fit_aux2sdb(g->aux, sdb);
	for(i=0;i<g->cnss->size;i++){
		cns = ref_cnsv(g->cnss, i);
		clear_basebank(cns->qc_seqs);
		clear_qc_alnv(cns->alns);
		for(j=0;j<cns->path->size;j++){
			p = ref_tracev(cns->path, j);
			n = ref_nodev(g->nodes, p->nid);
			bits2basebank(cns->qc_seqs, g->rdseqs->bits, p->nid * ((uint64_t)g->rd_len), g->rd_len);
			q = next_ref_qc_alnv(cns->alns);
			q->rid = n->rid;
			q->idx = cns->alns->size - 1;
			q->dir = n->seqdir;
			q->off = p->off;
			q->fix = 1;
			sr_block_aux(g->aux, n->rid);
		}
		sr_clear_aux(g->aux);
		sr_aux_load2(g->aux, sdb, 0, cns->seq->string, cns->len);
		sr_align_sdb(sdb, g->aux);
		for(j=0;j<g->aux->hits->size;j++){
			h = ref_sr_hitv(g->aux->hits, j);
			q = next_ref_qc_alnv(cns->alns);
			q->rid = h->rid;
			q->idx = cns->alns->size - 1;
			if(h->dir1 ^ h->dir2){
				revbits2basebank(cns->qc_seqs, sdb->rdseqs->bits, sr_rdseq_offset(sdb, h->rid), g->rd_len);
				q->dir = 1;
			} else {
				bits2basebank(cns->qc_seqs, sdb->rdseqs->bits, sr_rdseq_offset(sdb, h->rid), g->rd_len);
				q->dir = 0;
			}
			if(h->dir1){
				q->off = cns->len - sdb->rd_len - h->off;
			} else {
				q->off = h->off;
			}
			q->fix = 0;
		}
		if(lib_cnts && lib_ins){
			sort_array(cns->alns->buffer, cns->alns->size, qc_aln_t, (a.rid > b.rid));
			clear_u32list(g->stack);
			for(j=1;j<cns->alns->size;j++){
				q1 = ref_qc_alnv(cns->alns, j - 1);
				q2 = ref_qc_alnv(cns->alns, j);
				if((q1->rid >> 1) != (q2->rid >>1)) continue;
				if(q1->dir == q2->dir) continue;
				if(q1->off < q2->off){
					if(q1->dir) continue;
					m = c = 0; while(c + lib_cnts->buffer[m] < q1->rid){ c += lib_cnts->buffer[m++]; }
					ins = get_u32list(lib_ins, m);
					ins += (q1->off) + (cns->len - (q2->off + g->rd_len));
				} else {
					if(q2->dir) continue;
					m = c = 0; while(c + lib_cnts->buffer[m] < q1->rid){ c += lib_cnts->buffer[m++]; }
					ins = get_u32list(lib_ins, m);
					ins += (q2->off) + (cns->len - (q1->off + g->rd_len));
				}
				push_u32list(g->stack, ins);
				q1->fix = 1;
				q2->fix = 1;
			}
			cns->qc = g->stack->size * 10000;
			if(g->stack->size >= 1){
				cns->qc += _mean_distance(g->stack);
			}
		}
		sort_array(cns->alns->buffer, cns->alns->size, qc_aln_t, (a.off > b.off));
		ins = 0;
		for(j=1;j<cns->alns->size;j++){
			q1 = ref_qc_alnv(cns->alns, j - 1);
			q2 = ref_qc_alnv(cns->alns, j);
			if(q2->off - q1->off > (int)ins) ins = q2->off - q1->off;
		}
		cns->qc += ins * 1000 * 10000;
		if(do_cns){
			cns->alns->size = apply_xchg_array(cns->alns->buffer, cns->alns->size, qc_aln_t, a.fix);
			zeros_cbv(cns->cbs);
			for(j=0;j<cns->alns->size;j++){
				q = ref_qc_alnv(cns->alns, j);
				off = (uint64_t)(q->idx) * g->rd_len;
				s = (q->off < 0)? - q->off : 0;
				e = (q->off + (int)g->rd_len > (int)cns->len)? cns->len - q->off : g->rd_len;
				for(k=s;k<e;k++){
					b = ref_cbv(cns->cbs, q->off + k);
					b->bases[bits2bit(cns->qc_seqs->bits, off + k)] ++;
				}
			}
			_make_consensus(cns);
		}
	}
}

uint32_t output_contigs_graph(Graph *g, FILE *cnsf, FILE *msaf){
	cns_t *cns;
	trace_t *p;
	node_t *n;
	qc_aln_t *q;
	uint32_t i, j, k, s, e;
	char *seqs;
	if(cnsf){
		for(i=0;i<g->cnss->size;i++){
			cns = ref_cnsv(g->cnss, i);
			fprintf(cnsf, ">ps%u_%u_%u_%010u n_sr=%u n_l=%u n_r=%u rdlen=%u n_rd=%u", g->ps_id, i, cns->len, cns->qc, g->n_sr, g->n_l, g->n_r, (uint32_t)g->rd_len, (uint32_t)cns->path->size);
			if(g->out_mode){
				fprintf(cnsf, " reads=");
				for(j=0;j<cns->path->size;j++){
					p = ref_tracev(cns->path, j);
					if(j) fprintf(cnsf, ",");
					n = ref_nodev(g->nodes, p->nid);
					fprintf(cnsf, "%u:%c:%u", n->rid, "+-"[n->seqdir], p->off);
				}
			}
			fprintf(cnsf, "\n%s\n", cns->seq->string);
		}
	}
	if(msaf){
		seqs = malloc(g->rd_len + 1);
		for(i=0;i<g->cnss->size;i++){
			cns = ref_cnsv(g->cnss, i);
			fprintf(msaf, ">ps%u_%u_%u_%010u n_sr=%u n_l=%u n_r=%u rdlen=%u n_rd=%u", g->ps_id, i, cns->len, cns->qc, g->n_sr, g->n_l, g->n_r, (uint32_t)g->rd_len, (uint32_t)cns->path->size);
			if(1){
				fprintf(msaf, " reads=");
				for(j=0;j<cns->path->size;j++){
					p = ref_tracev(cns->path, j);
					n = ref_nodev(g->nodes, p->nid);
					if(j) fprintf(msaf, ",");
					fprintf(msaf, "%u:%c:%u", n->rid, "+-"[n->seqdir], p->off);
				}
			}
			fprintf(msaf, "\n%s\n", cns->seq->string);
			for(j=0;j<cns->path->size;j++){
				p = ref_tracev(cns->path, j);
				n = ref_nodev(g->nodes, p->nid);
				seq_basebank(g->rdseqs, p->nid * g->rd_len, g->rd_len, seqs);
				for(k=0;k<g->rd_len;k++){
					if(seqs[k] != cns->seq->string[k + p->off]){
						seqs[k] = seqs[k] + 'a' - 'A';
					}
				}
				for(k=0;(int)k<p->off;k++) fprintf(msaf, "-");
				fprintf(msaf, "%s", seqs);
				for(k=p->off+g->rd_len;k<cns->len;k++) fprintf(msaf, "-");
				fprintf(msaf, "\t%010u\t%c\n", n->rid, "+-"[n->seqdir]);
			}
			if(cns->alns->size == 0) continue;
			fprintf(msaf, "\n%s\n", cns->seq->string);
			for(j=0;j<cns->alns->size;j++){
				q = ref_qc_alnv(cns->alns, j);
				seq_basebank(cns->qc_seqs, ((uint64_t)q->idx) * g->rd_len, g->rd_len, seqs);
				s = (q->off < 0)? - q->off : 0;
				e = (q->off + (int)g->rd_len > (int)cns->len)? cns->len - q->off : g->rd_len;
				for(k=s;k<e;k++){
					if(seqs[k] != cns->seq->string[k + q->off]){ seqs[k] = seqs[k] + 'a' - 'A'; }
				}
				seqs[e] = 0;
				for(k=0;(int)k<q->off;k++) fprintf(msaf, "-");
				fprintf(msaf, "%s", seqs + s);
				for(k=q->off+g->rd_len;k<cns->len;k++) fprintf(msaf, "-");
				fprintf(msaf, "\t%010u\t%c\n", q->rid, "+-"[q->dir]);
			}
		}
		free(seqs);
	}
	return g->cnss->size;
}

