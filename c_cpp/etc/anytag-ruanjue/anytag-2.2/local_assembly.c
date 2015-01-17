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
 
#include "local_assembly.h"

int cmp_trace(const void *e1, const void *e2, void *ref){
	Trace *t1, *t2;
	t1 = (Trace*)e1;
	t2 = (Trace*)e2;
	if(t1->n_ol == t2->n_ol) return 0;
	if(t1->n_ol <  t2->n_ol) return 1;
	else return -1;
	ref = ref; // useless
}

void init2_lgraph(LGraph *g, uint8_t kmer_size, uint8_t rd_len, uint8_t min_ol, float min_sm, uint8_t max_mm, uint32_t min_ins, uint32_t max_ins){
	g->sdb         = sr_init_sdb(kmer_size, 4, rd_len);
	sr_set_align_parameters(g->sdb, 1, min_ol, min_sm, max_mm, 0);
	sr_set_filter_parameters(g->sdb, 0, 1024, 0);
	g->aux         = sr_init_aux();
	g->nodes       = init_nodev(1024);
	g->m_node      = 0;
	g->seqnodes    = init_snv(1024);
	g->m_seqnode   = 0;
	g->paths       = init_vplist(12);
	g->trackers[0] = malloc(sizeof(Tracker));
	g->trackers[0]->trace_dir = 0;
	g->trackers[0]->traces    = init_tracev(1024);
	g->trackers[0]->m_trace   = 0;
	g->trackers[0]->heap      = init_heap(cmp_trace, g);
	g->trackers[1] = malloc(sizeof(Tracker));
	g->trackers[1]->trace_dir = 1;
	g->trackers[1]->traces    = init_tracev(1024);
	g->trackers[1]->m_trace   = 0;
	g->trackers[1]->heap      = init_heap(cmp_trace, g);
	g->p1          = 0;
	g->p2          = 1;
	g->lays        = init_layv(1024);
	g->skeleton    = init_string(1024);
	g->supps       = init_suppv(1024);
	g->sp_sdb      = NULL;
	g->cbases      = init_cbv(1024);
	g->bcovs       = init_u32list(1024);
	g->pcovs       = init_u32list(1024);
	g->cns_seq     = init_string(1024);
	g->ords        = init_u32list(1024);
	g->bits        = init_u64list(1024);
	g->cns_id      = 0;
	g->cns_len     = 0;
	g->check_ins   = 0;
	g->min_ins     = min_ins;
	g->max_ins     = max_ins;
	g->flags       = init_bitvec(1024);
}

LGraph* init_lgraph(uint8_t kmer_size, uint8_t rd_len, uint8_t min_ol, float min_sm, uint8_t max_mm, uint32_t min_ins, uint32_t max_ins){
	LGraph *g;
	g = malloc(sizeof(LGraph));
	init2_lgraph(g, kmer_size, rd_len, min_ol, min_sm, max_mm, min_ins, max_ins);
	return g;
}

void trunoff_inserts_checking_lgraph(LGraph *g){ g->check_ins = 0; }
void trunon_inserts_checking_lgraph(LGraph *g){ g->check_ins = 1; }

void free2_lgraph(LGraph *g){
	uint32_t i;
	SeqNode *sn;
	sr_free_sdb(g->sdb);
	sr_free_aux(g->aux);
	for(i=0;i<g->m_node;i++){
		free_edgev(ref_nodev(g->nodes, i)->edges[0]);
		free_edgev(ref_nodev(g->nodes, i)->edges[1]);
	}
	free_nodev(g->nodes);
	for(i=0;i<g->m_seqnode;i++){
		sn = ref_snv(g->seqnodes, i);
		free_layv(sn->lays);
		free_sev(sn->edges[0]);
		free_sev(sn->edges[1]);
		free_seqhitv(sn->hits);
	}
	free_snv(g->seqnodes);
	for(i=0;i<g->paths->size;i++) free_u32list((u32list*)ref_vplist(g->paths, i));
	free_vplist(g->paths);
	free_layv(g->lays);
	free_string(g->skeleton);
	free_u32list(g->ords);
	free_cbv(g->cbases);
	free_u32list(g->bcovs);
	free_u32list(g->pcovs);
	free_string(g->cns_seq);
	free_u64list(g->bits);
	for(i=0;i<g->trackers[0]->m_trace;i++) free(get_tracev(g->trackers[0]->traces, i));
	free_tracev(g->trackers[0]->traces);
	for(i=0;i<g->trackers[1]->m_trace;i++) free(get_tracev(g->trackers[1]->traces, i));
	free_tracev(g->trackers[1]->traces);
	free_heap(g->trackers[0]->heap);
	free_heap(g->trackers[1]->heap);
	free(g->trackers[0]);
	free(g->trackers[1]);
	free_bitvec(g->flags);
}

void free_lgraph(LGraph *g){
	free2_lgraph(g);
	free(g);
}

void reset_lgraph(LGraph *g){
	sr_reset_sdb(g->sdb);
	clear_nodev(g->nodes);
	clear_snv(g->seqnodes);
	clear_sr_hitv(g->aux->hits);
	clear_layv(g->lays);
	clear_string(g->cns_seq);
	clear_suppv(g->supps);
	clear_string(g->skeleton);
	clear_string(g->cns_seq);
	g->cns_len = 0;
	g->cns_id = 0;
}

void renew_lgraph(LGraph *g){
	uint32_t ksize, rd_len, min_ol, max_mm, check_ins, min_ins, max_ins;
	float min_sm;
	ksize = g->sdb->kmer_size;
	rd_len = g->sdb->rd_len;
	min_ol = g->sdb->min_overlap;
	max_mm = g->sdb->max_mismatch;
	min_ins = g->min_ins;
	max_ins = g->max_ins;
	min_sm = g->sdb->min_similarity;
	check_ins = g->check_ins;
	free2_lgraph(g);
	init2_lgraph(g, ksize, rd_len, min_ol, min_sm, max_mm, min_ins, max_ins);
	g->check_ins = check_ins;
}

void push_lgraph(LGraph *g, uint32_t rid, char *seq, uint8_t seqlen, uint32_t seqdir, uint32_t anchor_type, int min_ins, int max_ins){
	Node *node;
	node = next_ref_nodev(g->nodes);
	if(count_nodev(g->nodes) > g->m_node){
		node->edges[0] = init_edgev(64);
		node->edges[1] = init_edgev(64);
		g->m_node = count_nodev(g->nodes);
	} else {
		clear_edgev(node->edges[0]);
		clear_edgev(node->edges[1]);
	}
	node->closed  = 0;
	node->rid     = rid;
	node->seqdir  = seqdir;
	node->len     = seqlen;
	node->anchor_type = anchor_type;
	node->min_ins = min_ins;
	node->max_ins = max_ins;
	node->snid    = 0;
	node->snoff   = 0;
	node->visit   = 0;
	sr_push_sdb(g->sdb, seq, seqlen);
}

void align_lgraph(LGraph *g){
	SR_AlnHit *hit;
	Node *n1, *n2;
	Edge *e;
	uint32_t i, offb, offe, rid1, rid2, n_ol;
	g->n_sr = g->nodes->size;
	if(count_nodev(g->nodes) < 2) return;
	g->cns_id = ref_nodev(g->nodes, 0)->rid >> 1;
	sr_ready_sdb(g->sdb);
	for(i=0;i<g->sdb->n_idx;i+=2) sr_index_sdb(g->sdb, i);
	clear_sr_hitv(g->aux->hits);
	encap_bitvec(g->flags, g->nodes->size * g->nodes->size);
	zeros_bitvec(g->flags);
	for(i=0;i<g->sdb->n_rd;i++) sr_align_sdb(g->sdb, i, g->aux);
	sort_array(g->aux->hits->buffer, g->aux->hits->size, SR_AlnHit, (((int)a.off) - ((int)b.off)));
	for(i=0;i<count_sr_hitv(g->aux->hits);i++){
		hit = ref_sr_hitv(g->aux->hits, i);
		//if(hit->dir1 ^ hit->dir2) continue; // nerver happen
		//sr_output_hit(g->sdb, hit, stdout);
		if(hit->off + ref_nodev(g->nodes, hit->rid2)->len <= ref_nodev(g->nodes, hit->rid1)->len && hit->rid2 != g->p1 && hit->rid2 != g->p2){
			ref_nodev(g->nodes, (hit->rid1 < hit->rid2)? hit->rid2 : hit->rid1)->closed = 1;
			continue;
		}
		n_ol = hit->n_ol;
		if(hit->dir1){
			rid1 = hit->rid2;
			rid2 = hit->rid1;
			offe = hit->off;
			offb = hit->off + ref_nodev(g->nodes, hit->rid2)->len - ref_nodev(g->nodes, hit->rid1)->len;
		} else {
			rid1 = hit->rid1;
			rid2 = hit->rid2;
			offb = hit->off;
			offe = hit->off + ref_nodev(g->nodes, hit->rid2)->len - ref_nodev(g->nodes, hit->rid1)->len;
		}
		one_bitvec(g->flags, rid1 * g->nodes->size + rid2);
		one_bitvec(g->flags, rid2 * g->nodes->size + rid1);
		n1 = ref_nodev(g->nodes, rid1);
		n2 = ref_nodev(g->nodes, rid2);
		if(n1->closed || n2->closed) continue;
		if(offb == 0 || offe == 0) continue;
		e = next_ref_edgev(n1->edges[0]); e->nid = rid2; e->n_ol = n_ol; e->off = offb; e->hit_idx = i; e->mut_idx = n2->edges[1]->size; e->closed = 0;
		e = next_ref_edgev(n2->edges[1]); e->nid = rid1; e->n_ol = n_ol; e->off = offe; e->hit_idx = i; e->mut_idx = n1->edges[0]->size - 1; e->closed = 0;
	}
}

uint32_t simplify_lgraph(LGraph *g){
	Node *n;
	Edge *e1, *e2;
	uint32_t i, j, k;
	uint32_t ret;
	ret = 0;
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		for(j=n->edges[0]->size;j>1;j--){
			e1 = ref_edgev(n->edges[0], j - 1);
			if(e1->closed) continue;
			for(k=0;k<j-1;k++){
				e2 = ref_edgev(n->edges[0], k);
				if(get_bitvec(g->flags, e1->nid * g->nodes->size + e2->nid)){
					e1->closed = 1;
					ref_edgev(ref_nodev(g->nodes, e1->nid)->edges[1], e1->mut_idx)->closed = 1;
					ret ++;
					break;
				}
			}
		}
	}
	return ret;
}

void print_dot_lgraph(LGraph *g, FILE *out){
	Node *n;
	Edge *e;
	uint32_t i, j;
	fprintf(out, "digraph G%u {\n", g->cns_id);
	fprintf(out, "rankdir=LR\n");
	fprintf(out, "N000 [shape=rect color=red]\n");
	fprintf(out, "N001 [shape=rect color=red]\n");
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->closed) continue;
		for(j=0;j<n->edges[0]->size;j++){
			e = ref_edgev(n->edges[0], j);
			if(e->closed) continue;
			fprintf(out, "N%03u -> N%03u [label=\"%u\"color=blue]\n", i, e->nid, e->off);
		}
	}
	fprintf(out, "}\n");
}

void push_tracker(Tracker *tracker, uint32_t src_nid, uint32_t hit_idx, uint32_t nid, uint32_t n_ol, uint32_t off){
	Trace *t;
	if(count_tracev(tracker->traces) >= tracker->m_trace){
		t = malloc(sizeof(Trace));
		push_tracev(tracker->traces, t);
		tracker->m_trace = count_tracev(tracker->traces);
	} else {
		t = get_tracev(tracker->traces, tracker->traces->size ++);
	}
	t->src_nid = src_nid;
	t->hit_idx = hit_idx;
	t->nid     = nid;
	t->n_ol    = n_ol;
	t->off     = off;
	push_heap(tracker->heap, t);
}

Trace* pop_tracker(Tracker *tracker){ return (Trace*)pop_heap(tracker->heap); }

uint32_t walk_along_lgraph(LGraph *g, layv *lays){
	Layout *lay1, *lay2;
	Node *n1, *n2;
	Edge *e1, *e2;
	uint32_t i, n_next, i_next, n_past, i_past;
	if(lays->size == 0) return 0;
	while(1){
		encap_layv(lays, 1);
		lay1 = ref_layv(lays, lays->size - 1);
		n1 = ref_nodev(g->nodes, lay1->nid);
		n_next = 0; i_next = 0;
		for(i=0;i<n1->edges[lay1->dir]->size;i++){
			e1 = ref_edgev(n1->edges[lay1->dir], i);
			if(e1->closed) continue;
			n_next ++;
			i_next = i;
			if(n_next > 1) break;
		}
		if(n_next != 1) break;
		e1 = ref_edgev(n1->edges[lay1->dir], i_next);
		n2 = ref_nodev(g->nodes, e1->nid);
		if(n2->visit) break;
		n_past = 0; i_past = 0;
		for(i=0;i<n2->edges[!lay1->dir]->size;i++){
			e2 = ref_edgev(n2->edges[!lay1->dir], i);
			if(e2->closed) continue;
			n_past ++;
			i_past = i;
			if(n_past > 1) break;
		}
		if(n_past != 1) break;
		e2 = ref_edgev(n2->edges[!lay1->dir], i_past);
		if(e2->nid != lay1->nid){
			fprintf(stdout, " -- Error in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__); fflush(stdout); abort();
		}
		lay2 = next_ref_layv(lays);
		lay2->nid = e1->nid;
		lay2->dir = lay1->dir;
		lay2->off = lay1->off + e1->off;
		lay2->len = n2->len;
		n2->visit = 1;
	}
	lay1 = ref_layv(lays, lays->size - 1);
	return lay1->off + lay1->len;
}

uint32_t translate_lgraph(LGraph *g){
	SeqNode *sn;
	SeqEdge *se;
	Node *n;
	Edge *e;
	Layout *lay;
	uint32_t i, j;
	clear_snv(g->seqnodes);
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		if(n->visit) continue;
		sn = next_ref_snv(g->seqnodes);
		if(g->seqnodes->size > g->m_seqnode){
			sn->lays = init_layv(32);
			sn->edges[0] = init_sev(4);
			sn->edges[1] = init_sev(4);
			sn->hits = init_seqhitv(1024);
			g->m_seqnode = g->seqnodes->size;
		} else {
			clear_layv(sn->lays);
			clear_sev(sn->edges[0]);
			clear_sev(sn->edges[1]);
			clear_seqhitv(sn->hits);
		}
		lay = next_ref_layv(sn->lays);
		lay->nid = i;
		lay->dir = 1;
		lay->off = 0;
		lay->len = n->len;
		sn->len = walk_along_lgraph(g, sn->lays);
		reverse_layv(sn->lays);
		ref_apply_array(sn->lays->buffer, sn->lays->size, Layout, (a->dir = !a->dir, a->off = sn->len - (a->off + a->len)));
		sn->len = walk_along_lgraph(g, sn->lays);
	}
	for(i=0;i<g->seqnodes->size;i++){
		sn = ref_snv(g->seqnodes, i);
		for(j=0;j<sn->lays->size;j++){
			lay = ref_layv(sn->lays, j);
			n = ref_nodev(g->nodes, lay->nid);
			n->snid = i;
			n->snoff = lay->off;
		}
	}
	for(i=0;i<g->seqnodes->size;i++){
		sn = ref_snv(g->seqnodes, i);
		lay = ref_layv(sn->lays, 0);
		n = ref_nodev(g->nodes, lay->nid);
		for(j=0;j<n->edges[1]->size;j++){
			e = ref_edgev(n->edges[1], j);
			if(e->closed) continue;
			se = next_ref_sev(sn->edges[1]);
			se->snid = ref_nodev(g->nodes, e->nid)->snid;
			se->n_ol = e->n_ol;
			se->mut_idx = 0;
			se->gap  = 0;
			se->closed = 0;
		}
		lay = ref_layv(sn->lays, sn->lays->size - 1);
		n = ref_nodev(g->nodes, lay->nid);
		for(j=0;j<n->edges[0]->size;j++){
			e = ref_edgev(n->edges[0], j);
			if(e->closed) continue;
			se = next_ref_sev(sn->edges[0]);
			se->snid = ref_nodev(g->nodes, e->nid)->snid;
			se->n_ol = e->n_ol;
			se->mut_idx = 0;
			se->gap  = 0;
			se->closed = 0;
		}
	}
	return g->seqnodes->size;
}

void gen_bitseq4seqnode(LGraph *g, uint32_t snid, u64list *bits){
	SeqNode *sn;
	Layout *lay;
	uint32_t i, j;
	uint32_t soff, slen, off, len, c, offset;
	sn = ref_snv(g->seqnodes, snid);
	encap_u64list(bits, (sn->len + 31) / 32);
	soff = slen = 0;
	for(i=0;i<sn->lays->size;i++){
		lay = ref_layv(sn->lays, i);
		if(lay->off + lay->len <= (int)slen) continue;
		off = slen - (soff + lay->off);
		len = lay->len - off;
		offset = sr_rdseq_offset(g->sdb, lay->nid);
		for(j=0;j<len;j++){
			c = bits2bit(g->sdb->rd_seqs->buffer, offset + j + off);
			bit2bits(bits->buffer, soff + j, c);
		}
		soff += lay->off;
		slen = soff + lay->len;
	}
}

uint32_t align_sr2lgraph(LGraph *g, SR_SeqDB *sdb){
	SR_AlnHit *h;
	SeqHit *hit;
	SeqNode *sn;
	uint32_t i, j, ret;
	ret = 0;
	for(i=0;i<g->seqnodes->size;i++){
		sn = ref_snv(g->seqnodes, i);
		gen_bitseq4seqnode(g, i, g->bits);
		clear_sr_hitv(g->aux->hits);
		sr_align_long_sdb(sdb, g->aux, i, g->bits->buffer, 0, sn->len, 3, 0);
		ret += g->aux->hits->size;
		for(j=0;j<g->aux->hits->size;j++){
			h = ref_sr_hitv(g->aux->hits, j);
			hit = next_ref_seqhitv(sn->hits);
			hit->rid = h->rid2;
			hit->dir = h->dir1 ^ h->dir2;
			hit->len = sr_rdseq_length(sdb, hit->rid);
			hit->mm  = h->n_mm;
			if(h->dir1){
				hit->off = sn->len - (h->off + hit->len);
			} else {
				hit->off = h->off;
			}
		}
		//sn->hits is already sorted by hit->rid in dsc order
	}
	return ret;
}

int is_illegal_path(LGraph *g, uint32_t nid, int path_off, uint32_t path_dir){
	int a, b, ret;
	Node *n;
	n = ref_nodev(g->nodes, nid);
	a = n->min_ins;
	b = n->max_ins;
	ret = 1;
	if((n->anchor_type & 0x01U) == path_dir){
		if(path_off < a || path_off > b) ret = 0;
	} else {
		if(((int)g->max_ins) - (path_off + n->len) < a) ret = 0;
		if(((int)g->min_ins) - (path_off + n->len) > b) ret = 0;
	}
	return ret;
}

int layout_lgraph(LGraph *g){
	Node *n;
	Edge *e;
	Trace *t1, *t2;
	Layout *l;
	uint32_t i, bp1, bp2, hid;
	if(count_nodev(g->nodes) < 2) return 0;
	clear_heap(g->trackers[0]->heap);
	clear_tracev(g->trackers[0]->traces);
	clear_heap(g->trackers[1]->heap);
	clear_tracev(g->trackers[1]->traces);
	n = ref_nodev(g->nodes, g->p1);
	n->visit  = 0;
	n->bt     = g->p1;
	n->bt_off = 0;
	n->bt_dir = 0;
	push_tracker(g->trackers[0], g->p1, 0, g->p1, n->len, 0);
	n = ref_nodev(g->nodes, g->p2);
	n->visit  = 0;
	n->bt     = g->p2;
	n->bt_off = 0;
	n->bt_dir = 1;
	push_tracker(g->trackers[1], g->p2, 0, g->p2, n->len, 0);
	bp1 = bp2 = 0;
	hid = 0;
	while((t1 = pop_tracker(g->trackers[0])) && (t2 = pop_tracker(g->trackers[1]))){
		n = ref_nodev(g->nodes, t1->nid);
		if(n->visit == 0){
			n->visit = 1;
			n->bt     = t1->src_nid;
			n->bt_idx = t1->hit_idx;
			n->bt_off = t1->off;
			n->bt_dir = 0;
			for(i=0;i<count_edgev(n->edges[0]);i++){
				e = ref_edgev(n->edges[0], i);
				if(ref_nodev(g->nodes, e->nid)->closed) continue;
				if(e->closed) continue;
				if(g->check_ins && is_illegal_path(g, e->nid, e->off + n->bt_off, 0) == 0) continue;
				push_tracker(g->trackers[0], t1->nid, e->hit_idx, e->nid, e->n_ol, e->off + n->bt_off);
			}
		} else if(n->bt_dir == 1){
			if(n->bt_off + t1->off + n->len >= g->min_ins && n->bt_off + t1->off + n->len <= g->max_ins){
				bp1 = t1->src_nid;
				bp2 = t1->nid;
				hid = t1->hit_idx;
				break;
			}
		}
		n = ref_nodev(g->nodes, t2->nid);
		if(n->visit == 0){
			n->visit = 1;
			n->bt     = t2->src_nid;
			n->bt_idx = t2->hit_idx;
			n->bt_off = t2->off;
			n->bt_dir = 1;
			for(i=0;i<count_edgev(n->edges[1]);i++){
				e = ref_edgev(n->edges[1], i);
				if(e->closed) continue;
				if(ref_nodev(g->nodes, e->nid)->closed) continue;
				if(g->check_ins && is_illegal_path(g, e->nid, e->off + n->bt_off, 1) == 0) continue;
				push_tracker(g->trackers[1], t2->nid, e->hit_idx, e->nid, e->n_ol, e->off + n->bt_off);
			}
		} else if(n->bt_dir == 0){
			if(n->bt_off + t2->off + n->len >= g->min_ins && n->bt_off + t2->off + n->len <= g->max_ins){
				bp1 = t2->nid;
				bp2 = t2->src_nid;
				hid = t2->hit_idx;
				break;
			}
		}
	}
	if(bp1 == bp2) return 0;
	clear_layv(g->lays);
	while(1){
		n = ref_nodev(g->nodes, bp1);
		l = next_ref_layv(g->lays);
		l->nid = bp1;
		l->hit_idx = n->bt_idx;
		if(n->bt == bp1) break;
		bp1 = n->bt;
	}
	reverse_layv(g->lays);
	while(1){
		n = ref_nodev(g->nodes, bp2);
		l = next_ref_layv(g->lays);
		l->nid = bp2;
		l->hit_idx = hid;
		if(n->bt == bp2) break;
		bp2 = n->bt;
		hid = n->bt_idx;
	}
	g->n_layout = g->lays->size;
	return 1;
}

void skeleton_lgraph(LGraph *g){
	Layout *lay;
	SR_AlnHit *hit;
	uint32_t i, slen, soff, off, len;
	clear_string(g->skeleton);
	soff = 0;
	slen = 0;
	for(i=0;i<g->lays->size;i++){
		lay = ref_layv(g->lays, i);
		if(i == 0){
			len = sr_rdseq_length(g->sdb, lay->nid);
			bits2seq(g->chs, ref_u64list(g->sdb->rd_seqs, 0), sr_rdseq_offset(g->sdb, lay->nid), len);
			append_string(g->skeleton, g->chs, len);
			soff = 0;
			slen = len;
			continue;
		}
		hit = ref_sr_hitv(g->aux->hits, lay->hit_idx);
		len = sr_rdseq_length(g->sdb, lay->nid);
		if(soff + hit->off + len <= slen) continue;
		off = slen - (soff + hit->off);
		bits2seq(g->chs, ref_u64list(g->sdb->rd_seqs, 0), sr_rdseq_offset(g->sdb, lay->nid) + off, len - off);
		append_string(g->skeleton, g->chs, len - off);
		slen += len - off;
		soff += hit->off;
	}
	g->cns_len = slen;
	clear_suppv(g->supps);
}

void get_skeleton_lgraph(LGraph *g, String *skeleton){
	append_string(skeleton, g->skeleton->string, g->skeleton->size);
}

void set_skeleton_lgraph(LGraph *g, String *skeleton){
	clear_string(g->skeleton);
	append_string(g->skeleton, skeleton->string, skeleton->size);
	g->cns_len = skeleton->size;
}

void push_aln_lgraph(LGraph *g, Supp *sp1, Supp *sp2){
	Supp *sp;
	sp = next_ref_suppv(g->supps);
	*sp = *sp1;
	sp = next_ref_suppv(g->supps);
	*sp = *sp2;
}

uint32_t align_skeleton_lgraph(LGraph *g, SR_SeqDB *sdb, SR_AlnAux *aux){
	SR_AlnHit *hit, *h;
	Supp *sp;
	uint32_t i, j, last_rid, slen, ret;
	g->sp_sdb = sdb;
	slen = g->cns_len;
	clear_u64list(g->bits);
	encap_u64list(g->bits, (slen + 31) / 32 + 32);
	seq2bits(ref_u64list(g->bits, 0), 0, g->skeleton->string, slen);
	clear_sr_hitv(aux->hits);
	sr_align_long_sdb(sdb, aux, g->cns_id, ref_u64list(g->bits, 0), 0, slen, 3, 2 * 1024);
	if(aux->hits->size < 2) return 0;
	ret = 0;
	last_rid = ref_sr_hitv(aux->hits, 0)->rid2;
	for(i=1;i<aux->hits->size;i++){
		hit = ref_sr_hitv(aux->hits, i);
		if(last_rid & 0x01){
			if(hit->rid2 == (last_rid & (~0x1U))){
				for(j=0;j<2;j++){
					h = ref_sr_hitv(aux->hits, i - j);
					sp = next_ref_suppv(g->supps);
					sp->rid = h->rid2;
					sp->len = sr_rdseq_length(sdb, sp->rid);
					if(h->dir1) sp->off = slen - (h->off + sr_rdseq_length(sdb, h->rid2));
					else sp->off = h->off;
					sp->dir = h->dir1 ^ h->dir2;
				}
				ret ++;
			}
		}
		last_rid = hit->rid2;
	}
	return ret;
}

void consensus_lgraph(LGraph *g){
	comb_base_t *cb;
	Supp *sp, *sp1, *sp2;
	uint64_t seqoff;
	uint32_t i, j, base, b, e, *ary;
	clear_u32list(g->ords);
	for(i=0;i<g->supps->size;i++) push_u32list(g->ords, i);
	ary = ref_u32list(g->ords, 0);
	sort_array(ary, g->ords->size, uint32_t, (((int)ref_suppv(g->supps, a)->off) - ((int)ref_suppv(g->supps, b)->off)));
	clear_u32list(g->bcovs);
	clear_u32list(g->pcovs);
	clear_cbv(g->cbases);
	for(i=0;i<g->cns_len;i++){
		cb = next_ref_cbv(g->cbases);
		cb->bases[0] = 0;
		cb->bases[1] = 0;
		cb->bases[2] = 0;
		cb->bases[3] = 0;
		push_u32list(g->bcovs, 0);
		push_u32list(g->pcovs, 0);
	}
	for(i=0;i<g->supps->size;i++){
		sp = ref_suppv(g->supps, i);
		seqoff = sr_rdseq_offset(g->sp_sdb, sp->rid);
		for(j=0;j<sp->len;j++){
			base = sp->dir? bits2revbit(ref_u64list(g->sp_sdb->rd_seqs, 0), (seqoff + sp->len - 1 - j))
					: bits2bit(ref_u64list(g->sp_sdb->rd_seqs, 0), (seqoff + j));
			ref_cbv(g->cbases, j + sp->off)->bases[base] ++;
			g->bcovs->buffer[j + sp->off] ++;
		}
	}
	for(i=0;i<g->supps->size;i+=2){
		sp1 = ref_suppv(g->supps, i);
		sp2 = ref_suppv(g->supps, i + 1);
		if(sp1->off < sp2->off){ b = sp1->off; e = sp2->off + sp2->len; }
		else { b = sp2->off; e = sp1->off + sp1->len; }
		for(j=b;j<e;j++) g->pcovs->buffer[j] ++;
	}
	clear_string(g->cns_seq);
	append_string(g->cns_seq, g->skeleton->string, g->skeleton->size);
	for(i=0;i<g->cns_len;i++){
		base = 0;
		cb = ref_cbv(g->cbases, i);
		for(j=1;j<4;j++){ if(cb->bases[j] > cb->bases[base]){ base = j; } }
		if(cb->bases[base] == 0) continue;
		g->cns_seq->string[i] = bit_base_table[base];
	}
	b = sr_rdseq_length(g->sdb, 0);
	e = g->cns_len - sr_rdseq_length(g->sdb, 1) + 1;
	if(b + 10 > e){ b = 0; e = g->cns_len; }
	g->min_bcov = 1000;
	g->min_pcov = 1000;
	for(i=b;i<e;i++){
		if(get_u32list(g->bcovs, i) < g->min_bcov) g->min_bcov = get_u32list(g->bcovs, i); 
		if(get_u32list(g->pcovs, i) < g->min_pcov) g->min_pcov = get_u32list(g->pcovs, i); 
	}
}

void output_lgraph(LGraph *g, FILE *seqf, FILE *msaf){
	uint32_t i, j, off1, off2, idx, off;
	uint64_t rdoff, c;
	Supp *sp;
	if(seqf){
		fprintf(seqf, ">ps%u len=%u n_rd=%u,%u,%u min_base_cov=%u min_mate_cov=%u\n", g->cns_id, g->cns_len,
			(unsigned)g->ords->size, g->n_layout, g->n_sr, g->min_bcov, g->min_pcov);
		for(i=0;i<g->cns_len;i++){
			fputc(g->cns_seq->string[i], seqf);
			if(((i + 1) % 100) == 0) fputc('\n', seqf);
		}
		if(i % 100) fputc('\n', seqf);
	}
	if(msaf){
		clear_u64list(g->bits);
		push_u64list(g->bits, (((uint64_t)0) << 32) | g->cns_id);
		push_u64list(g->bits, (((uint64_t)g->cns_len) << 32) | g->ords->size);
		push_u64list(g->bits, (((uint64_t)g->min_bcov) << 32) | g->min_pcov);
		push_u64list(g->bits, (((uint64_t)g->n_sr) << 32) | g->n_layout);
		for(i=0;i<g->ords->size;i++){
			sp = ref_suppv(g->supps, get_u32list(g->ords, i));
			push_u64list(g->bits, *((uint64_t*)sp));
		}
		off1 = 4 + g->ords->size;
		off2 = 0;
		encap_u64list(g->bits, (off2 + g->cns_len + 31) / 32);
		seq2bits(ref_u64list(g->bits, off1), off2, g->cns_seq->string, g->cns_len);
		off2 += g->cns_len;
		for(i=0;i<g->ords->size;i++){
			sp = ref_suppv(g->supps, get_u32list(g->ords, i));
			encap_u64list(g->bits, (off2 + sp->len + 31) / 32);
			rdoff = sr_rdseq_offset(g->sp_sdb, sp->rid);
			if(sp->dir){
				for(j=0;j<sp->len;j++){
					idx = off1 + (off2 >> 5);
					off = off2 & 0x1FU;
					off2 ++;
					c = bits2revbit(ref_u64list(g->sp_sdb->rd_seqs, 0), rdoff + sp->len - j - 1);
					if(off == 0) g->bits->buffer[idx] = 0;
					g->bits->buffer[idx] |= c << (((~off) & 0x1FU) << 1);
				}
			} else {
				for(j=0;j<sp->len;j++){
					idx = off1 + (off2 >> 5);
					off = off2 & 0x1FU;
					off2 ++;
					c = bits2bit(ref_u64list(g->sp_sdb->rd_seqs, 0), rdoff + j);
					if(off == 0) g->bits->buffer[idx] = 0;
					g->bits->buffer[idx] |= c << (((~off) & 0x1FU) << 1);
				}
			}
		}
		idx = off1 + (off2 + 31) / 32;
		g->bits->buffer[0] |= ((uint64_t)idx) << 32;
		fwrite(g->bits->buffer, sizeof(uint64_t), idx, msaf);
	}
}

