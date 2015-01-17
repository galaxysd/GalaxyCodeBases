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

void init2_lgraph(LGraph *g, ATOptions *opt){
	uint32_t i;
	g->opt         = opt;
	g->sdb         = sr_init2_sdb(opt->kmer_size[1], opt->n_seed[1], opt->rd_len, 0, 0.77);
	for(i=0;i<opt->n_seed[1];i++){
		if((opt->seed_skips[1] >> i) & 0x01) sr_mask_seed_sdb(g->sdb, i);
	}
	g->aux         = sr_init_aux(1, MAX_RD_LEN, opt->min_ol[1], opt->max_mm[1], opt->min_sm[1], opt->allow_gap[1], 0, opt->max_hits[1]);
	g->long_aux    = sr_init_aux(3, MAX_RD_LEN, opt->min_ol[2], opt->max_mm[2], opt->min_sm[2], opt->allow_gap[2] & 0x01, 0, opt->max_hits[2]);
	g->min_ol      = opt->min_ol[1];
	g->ol_step     = opt->ol_step;
	g->cur_ol[0]   = 0;
	g->cur_ol[1]   = 0;
	g->nodes       = init_nodev(16);
	g->m_node      = 0;
	g->trackers[0] = malloc(sizeof(Tracker));
	g->trackers[0]->trace_dir = 0;
	g->trackers[0]->traces    = init_tracev(16);
	g->trackers[0]->m_trace   = 0;
	g->trackers[0]->heap      = init_heap(cmp_trace, g);
	g->trackers[1] = malloc(sizeof(Tracker));
	g->trackers[1]->trace_dir = 1;
	g->trackers[1]->traces    = init_tracev(16);
	g->trackers[1]->m_trace   = 0;
	g->trackers[1]->heap      = init_heap(cmp_trace, g);
	g->p1          = 0;
	g->p2          = 1;
	g->n_l         = 0;
	g->n_r         = 0;
	g->lays        = init_layv(16);
	g->skeleton    = init_string(16);
	g->sp_hits     = init_sr_hitv(16);
	g->sp_sdb      = NULL;
	g->cbases      = init_cbv(16);
	g->cigars[0]   = init_cigarm(16);
	g->cigars[1]   = init_cigarm(16);
	g->msa[0]      = init_strm(16);
	g->msa[1]      = init_strm(16);
	g->bcovs       = init_u32list(16);
	g->pcovs       = init_u32list(16);
	g->cns_seq     = init_string(16);
	g->ords[0]     = init_u32list(16);
	g->ords[1]     = init_u32list(16);
	g->bits        = init_u64list(16);
	g->cns_id      = 0;
	g->cns_len     = 0;
	g->check_ins   = 0;
	g->min_ins     = opt->min_ins;
	g->max_ins     = opt->max_ins;
	g->flags       = init_bitvec(64);
}

LGraph* init_lgraph(ATOptions *opt){
	LGraph *g;
	g = malloc(sizeof(LGraph));
	init2_lgraph(g, opt);
	return g;
}

void trunoff_inserts_checking_lgraph(LGraph *g){ g->check_ins = 0; }
void trunon_inserts_checking_lgraph(LGraph *g){ g->check_ins = 1; }
void free2_lgraph(LGraph *g){
	uint32_t i;
	sr_free_sdb(g->sdb);
	sr_free_aux(g->aux);
	sr_free_aux(g->long_aux);
	for(i=0;i<g->m_node;i++){
		free_edgev(ref_nodev(g->nodes, i)->edges[0]);
		free_edgev(ref_nodev(g->nodes, i)->edges[1]);
	}
	free_nodev(g->nodes);
	free_layv(g->lays);
	free_string(g->skeleton);
	free_u32list(g->ords[0]);
	free_u32list(g->ords[1]);
	free_cbv(g->cbases);
	for(i=0;i<g->cigars[0]->size;i++) free_cigarv(get_cigarm(g->cigars[0], i));
	free_cigarm(g->cigars[0]);
	for(i=0;i<g->cigars[1]->size;i++) free_cigarv(get_cigarm(g->cigars[1], i));
	free_cigarm(g->cigars[1]);
	for(i=0;i<g->msa[0]->size;i++) free_b8list(get_strm(g->msa[0], i));
	free_strm(g->msa[0]);
	for(i=0;i<g->msa[1]->size;i++) free_b8list(get_strm(g->msa[1], i));
	free_strm(g->msa[1]);
	free_sr_hitv(g->sp_hits);
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
	sr_clear_aux(g->aux);
	clear_layv(g->lays);
	clear_string(g->cns_seq);
	clear_sr_hitv(g->sp_hits);
	clear_string(g->skeleton);
	clear_string(g->cns_seq);
	g->cns_len = 0;
	g->cns_id = 0;
	g->n_l = 0;
	g->n_r = 0;
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
	if(anchor_type) g->n_r ++;
	else g->n_l ++;
	node->min_ins = min_ins;
	node->max_ins = max_ins;
	node->snid    = 0;
	node->snoff   = 0;
	node->visit   = 0;
	node->last_n_edges[0] = 0;
	node->last_n_edges[1] = 0;
	sr_push_sdb(g->sdb, seq, seqlen);
}

void ready_lgraph(LGraph *g){
	Node *n;
	uint32_t i;
	g->n_sr = g->nodes->size;
	if(count_nodev(g->nodes) < 2) return;
	//if(count_nodev(g->nodes) > 2 * g->opt->limit) return;
	g->cns_id = ref_nodev(g->nodes, 0)->rid >> 1;
	if(g->opt->flags[1] & 0x02) fprintf(stdout, ">alignments%d\n", g->cns_id);
	g->cur_ol[0] = 0;
	g->cur_ol[1] = 0;
	sr_ready_sdb(g->sdb);
	sr_fit_aux2sdb(g->aux, g->sdb);
	sr_clear_aux(g->aux);
	for(i=0;i<g->sdb->n_idx;i+=2) sr_index_sdb(g->sdb, i);
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		n->last_n_edges[0] = 0;
		n->last_n_edges[1] = 0;
	}
}

uint32_t add_hits2edges(LGraph *g, SR_AlnHit *hit, uint32_t hit_idx){
	Node *n1, *n2;
	Edge *e;
	uint32_t offb, offe, rid1, rid2, n_ol;
	if(hit->off + ref_nodev(g->nodes, hit->rid2)->len <= ref_nodev(g->nodes, hit->rid1)->len && hit->rid2 != g->p1 && hit->rid2 != g->p2){
		ref_nodev(g->nodes, (hit->rid1 < hit->rid2)? hit->rid2 : hit->rid1)->closed = 1;
		return 0;
	}
	n_ol = hit->n_ol;
	if(hit->dir1){
		rid1 = hit->rid2;
		rid2 = hit->rid1;
		offe = hit->off;
		if(hit->off + ref_nodev(g->nodes, hit->rid2)->len < ref_nodev(g->nodes, hit->rid1)->len){
			ref_nodev(g->nodes, rid1)->closed = 1; return 0;
		}
		offb = hit->off + ref_nodev(g->nodes, hit->rid2)->len - ref_nodev(g->nodes, hit->rid1)->len;
	} else {
		rid1 = hit->rid1;
		rid2 = hit->rid2;
		offb = hit->off;
		if(hit->off + ref_nodev(g->nodes, hit->rid2)->len < ref_nodev(g->nodes, hit->rid1)->len){
			ref_nodev(g->nodes, rid2)->closed = 1; return 0;
		}
		offe = hit->off + ref_nodev(g->nodes, hit->rid2)->len - ref_nodev(g->nodes, hit->rid1)->len;
	}
	n1 = ref_nodev(g->nodes, rid1);
	n2 = ref_nodev(g->nodes, rid2);
	if(n1->closed || n2->closed) return 0;
	if(offb == 0 || offe == 0) return 0;
	e = next_ref_edgev(n1->edges[0]); e->nid = rid2; e->n_ol = n_ol; e->off = offb; e->hit_idx = hit_idx; e->mut_idx = n2->edges[1]->size; e->closed = 0;
	e = next_ref_edgev(n2->edges[1]); e->nid = rid1; e->n_ol = n_ol; e->off = offe; e->hit_idx = hit_idx; e->mut_idx = n1->edges[0]->size - 1; e->closed = 0;
	return 1;
}

uint32_t next_align_lgraph(LGraph *g){
	SR_AlnHit *hit;
	Node *n;
	uint32_t i, j, beg;
	if(g->cur_ol[0] == 0){ g->cur_ol[0] = g->sdb->max_rd_len; }
	else if(g->cur_ol[1] <= g->min_ol){ return 0; }
	else { g->cur_ol[0] = g->cur_ol[1] - 1; }
	g->cur_ol[1] = (g->cur_ol[0] < g->ol_step + g->min_ol)? g->min_ol : g->cur_ol[0] - g->ol_step;
	sr_set_overlap_aux(g->aux, g->cur_ol[0], g->cur_ol[1]);
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		for(j=0;j<n->edges[0]->size;j++) sr_block_aux(g->aux, ref_edgev(n->edges[0], j)->nid);
		for(j=0;j<n->edges[1]->size;j++) sr_block_aux(g->aux, ref_edgev(n->edges[1], j)->nid);
		beg = g->aux->hits->size;
		sr_align_sdb(g->sdb, i, g->aux);
		for(j=beg;j<count_sr_hitv(g->aux->hits);j++){
			hit = ref_sr_hitv(g->aux->hits, j);
			if(g->opt->flags[1] & 0x02) sr_output_hit(g->sdb, hit, stdout);
			add_hits2edges(g, hit, j);
		}
	}
	return 1;
}

uint32_t simplify_lgraph(LGraph *g){
	Node *n;
	Edge *e1, *e2;
	uint32_t i, j, k;
	uint32_t ret;
	ret = 0;
	while(next_align_lgraph(g));
	encap_bitvec(g->flags, g->nodes->size * g->nodes->size);
	zeros_bitvec(g->flags);
	for(i=0;i<g->nodes->size;i++){
		n = ref_nodev(g->nodes, i);
		sort_array(n->edges[0]->buffer, n->edges[0]->size, Edge, (((int)a.off) - ((int)b.off)));
		for(j=0;j<n->edges[0]->size;j++){
			e1 = ref_edgev(n->edges[0], j);
			if(e1->closed) continue;
			one_bitvec(g->flags, i * g->nodes->size + e1->nid);
		}
		for(j=0;j<n->edges[1]->size;j++){
			e1 = ref_edgev(n->edges[1], j);
			if(e1->closed) continue;
			one_bitvec(g->flags, i * g->nodes->size + e1->nid);
		}
	}
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
	uint32_t i, j, bp1, bp2, hid;
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
	t1 = NULL;
	t2 = NULL;
	while(1){
		if(t1 == NULL) t1 = pop_tracker(g->trackers[0]);
		if(t2 == NULL) t2 = pop_tracker(g->trackers[1]);
		if(t1 == NULL || t2 == NULL){
			if(next_align_lgraph(g) == 0) break;
			else {
				for(i=0;i<g->nodes->size;i++){
					n = ref_nodev(g->nodes, i);
					if(!n->visit) continue;
					for(j=n->last_n_edges[n->bt_dir];j<n->edges[n->bt_dir]->size;j++){
						e = ref_edgev(n->edges[n->bt_dir], j);
						if(ref_nodev(g->nodes, e->nid)->closed) continue;
						if(e->closed) continue;
						if(g->check_ins && is_illegal_path(g, e->nid, e->off + n->bt_off, 0) == 0) continue;
						push_tracker(g->trackers[n->bt_dir], i, e->hit_idx, e->nid, e->n_ol, e->off + n->bt_off);
					}
				}
				continue;
			}
		}
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
			n->last_n_edges[0] = n->edges[0]->size;
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
			n->last_n_edges[1] = n->edges[1]->size;
		} else if(n->bt_dir == 0){
			if(n->bt_off + t2->off + n->len >= g->min_ins && n->bt_off + t2->off + n->len <= g->max_ins){
				bp1 = t2->nid;
				bp2 = t2->src_nid;
				hid = t2->hit_idx;
				break;
			}
		}
		t1 = t2 = NULL;
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
	comb_base_t *cb;
	Layout *lay, *laz;
	cigarv *cv;
	b8list *str;
	SR_AlnHit *hit;
	uint64_t seqoff;
	uint32_t i, j;
	uint8_t base;
	clear_string(g->skeleton);
	if(g->lays->size == 0) return;
	clear_u32list(g->ords[0]);
	for(i=0;i<g->lays->size;i++) push_u32list(g->ords[0], i);
	for(i=g->cigars[0]->size;i<=g->lays->size+1;i++) push_cigarm(g->cigars[0], init_cigarv(1024));
	for(i=g->msa[0]->size;i<=g->lays->size+1;i++) push_strm(g->msa[0], init_b8list(1024));
	laz = ref_layv(g->lays, 0);
	laz->len = ref_nodev(g->nodes, laz->nid)->len;
	laz->off = 0;
	g->cns_len = laz->off + laz->len;
	for(i=1;i<g->lays->size;i++){
		lay = ref_layv(g->lays, i);
		lay->len = ref_nodev(g->nodes, lay->nid)->len;
		hit = ref_sr_hitv(g->aux->hits, lay->hit_idx);
		if(hit->rid1 == laz->nid){
			lay->off = laz->off + hit->off;
		} else {
			lay->off = laz->off + hit->off + ref_nodev(g->nodes, hit->rid2)->len - ref_nodev(g->nodes, hit->rid1)->len;
		}
		laz = lay;
		if((uint32_t)(laz->off + laz->len) > g->cns_len) g->cns_len = laz->off + laz->len;
	}
	clear_cbv(g->cbases);
	for(i=0;i<g->cns_len;i++){
		cb = next_ref_cbv(g->cbases); 
		cb->bases[0] = cb->bases[1] = cb->bases[2] = cb->bases[3] = 0;
	}
	for(i=0;i<g->lays->size;i++){
		lay = ref_layv(g->lays, i);
		cv = get_cigarm(g->cigars[0], i);
		clear_cigarv(cv);
		if(lay->off) cv->size = append_cigars(cv->buffer, cv->size, ALN_CIGAR_TYPE_INS, lay->off);
		cv->size = append_cigars(cv->buffer, cv->size, ALN_CIGAR_TYPE_MAT, lay->len);
		if(lay->off + lay->len < (int)g->cns_len) cv->size = append_cigars(cv->buffer, cv->size, ALN_CIGAR_TYPE_INS, g->cns_len - lay->off - lay->len);
		seqoff = sr_rdseq_offset(g->sdb, lay->nid);
		bits2seq(g->chs[0], g->sdb->rd_seqs->buffer, seqoff, lay->len);
		str = get_strm(g->msa[0], i);
		clear_b8list(str);
		encap_b8list(str, g->cns_len);
		str->size = cigars_seq2aln((char*)str->buffer, cv->buffer, cv->size, 1, g->chs[0]);
		for(j=0;j<lay->len;j++){
			base = bits2bit(ref_u64list(g->sdb->rd_seqs, 0), (seqoff + j));
			ref_cbv(g->cbases, lay->off + j)->bases[base] ++;
		}
	}
	encap_string(g->skeleton, g->cns_len);
	for(i=0;i<g->cns_len;i++){
		cb = ref_cbv(g->cbases, i);
		for(j=1,base=0;j<4;j++){ if(cb->bases[j] > cb->bases[base]){ base = j; } }
		g->skeleton->string[i] = bit_base_table[base];
	}
	g->skeleton->size = i;
	g->skeleton->string[i] = 0;
	for(i=0;i<g->lays->size;i++){
		str = get_strm(g->msa[0], i);
		for(j=0;j<g->cns_len;j++){
			if(str->buffer[j] == '-') continue;
			if(str->buffer[j] == ' '){ str->buffer[j] = '~'; continue; }
			if(str->buffer[j] != g->skeleton->string[j]) str->buffer[j] += 'a' - 'A';
		}
	}
	clear_sr_hitv(g->sp_hits);
}

void get_skeleton_lgraph(LGraph *g, String *skeleton){
	append_string(skeleton, g->skeleton->string, g->skeleton->size);
}

void set_skeleton_lgraph(LGraph *g, String *skeleton){
	clear_string(g->skeleton);
	append_string(g->skeleton, skeleton->string, skeleton->size);
	g->cns_len = skeleton->size;
}

void push_aln_lgraph(LGraph *g, SR_AlnHit *sp1, SR_AlnHit *sp2){
	push_sr_hitv(g->sp_hits, *sp1);
	push_sr_hitv(g->sp_hits, *sp2);
}

// implemented in sr_aln.c
int sr_aln_core_dp(SR_SeqDB *sdb, SR_AlnAux *aux, SR_AlnPtr *ptr);

uint32_t align_skeleton_lgraph(LGraph *g, SR_SeqDB *sdb){
	SR_AlnHit *hit, *h;
	SR_AlnPtr PTR;
	SR_AlnAux *aux;
	uint32_t i, j, last_rid, slen, size, ret;
	ret = 0;
	aux = g->long_aux;
	g->sp_sdb = sdb;
	slen = g->cns_len;
	clear_u64list(g->bits);
	encap_u64list(g->bits, (slen + 31) / 32 + 32);
	seq2bits(ref_u64list(g->bits, 0), 0, g->skeleton->string, slen);
	sr_fit_aux2sdb(aux, sdb);
	sr_clear_aux(aux);
	sr_rd_flags_mode_aux(aux, 0);
	sr_align_long_sdb(sdb, aux, g->cns_id, ref_u64list(g->bits, 0), 0, slen);
	if(aux->hits->size == 0) goto RET;
	for(i=0;i<aux->hits->size;i++){
		h = ref_sr_hitv(aux->hits, i);
		if(h->dir1){
			h->dir1 = 0; h->dir2 = !h->dir2;
			h->off  = g->cns_len - h->off - sr_rdseq_length(sdb, h->rid2);
			reverse_cigars(h->cigars, h->n_cigar);
		}
	}
	sort_array(aux->hits->buffer, aux->hits->size, SR_AlnHit, (((int64_t)b.rid2) - ((int64_t)a.rid2)));
	last_rid = ref_sr_hitv(aux->hits, 0)->rid2;
	for(i=1;i<aux->hits->size;i++){
		hit = ref_sr_hitv(aux->hits, i);
		if(last_rid & 0x01){
			if(hit->rid2 == (last_rid & (~0x1U))){
				for(j=2;j>0;j--){
					h = ref_sr_hitv(aux->hits, i - (j - 1));
					push_sr_hitv(g->sp_hits, *h);
					h->n_ol = 0;
				}
				ret ++;
			}
		}
		last_rid = hit->rid2;
	}
	if((g->opt->allow_gap[2] & 0x02) == 0){
		PTR.dir1 = 0;
		PTR.fix  = 1;
		for(i=0,size=aux->hits->size;i<size;i++){
			encap_sr_hitv(aux->hits, 1);
			hit = ref_sr_hitv(aux->hits, i);
			if(hit->n_ol == 0) continue;
			if(get_bitvec(aux->rd_flags, hit->rid2 ^ 0x01U) == 0) continue;
			PTR.dir2 = !hit->dir2;
			PTR.ptr  = (hit->rid2 ^ 0x01U) - sdb->rd_off;
			if(sr_aln_core_dp(sdb, aux, &PTR) == 0) continue;
			h = ref_sr_hitv(aux->hits, aux->hits->size - 1);
			if((h->rid2 ^ hit->rid2) != 0x01U) continue;
			//sr_output_hit(sdb, h, stdout);
			if(hit->rid2 & 0x01){
				push_sr_hitv(g->sp_hits, *h);
				push_sr_hitv(g->sp_hits, *hit);
			} else {
				push_sr_hitv(g->sp_hits, *hit);
				push_sr_hitv(g->sp_hits, *h);
			}
		}
	}
RET:
	sr_clear_rd_flags_aux(aux);
	sr_rd_flags_mode_aux(aux, 1);
	return ret;
}

void consensus_lgraph(LGraph *g){
	SR_AlnHit *sp, *sp1, *sp2;
	cigarv *cv0, *cv1, *cv2;
	AlnCigarIter *iter;
	b8list *str, *cns;
	uint32_t i, j, base, b, e, cns_len, acgt[6];
	clear_u32list(g->ords[1]);
	for(i=0;i<g->sp_hits->size;i++) push_u32list(g->ords[1], i);
	for(i=g->cigars[1]->size;i<=g->sp_hits->size+1;i++) push_cigarm(g->cigars[1], init_cigarv(1024));
	cv0 = get_cigarm(g->cigars[1], g->sp_hits->size);
	clear_cigarv(cv0); cv0->size = append_cigars(cv0->buffer, cv0->size, ALN_CIGAR_TYPE_MAT, g->cns_len);
	cv2 = get_cigarm(g->cigars[1], g->sp_hits->size + 1);
	iter = init_iter_2cigars(NULL, 0, NULL, 0);
	for(i=0;i<g->sp_hits->size;i++){
		sp = ref_sr_hitv(g->sp_hits, i);
		if(sp->is_gap == 0) continue;
		reset_iter_2cigars(iter, cv0->buffer, cv0->size, sp->cigars, sp->n_cigar);
		clear_cigarv(cv2);
		while(iter_2cigars(iter)){
			if(iter->type1 == iter->type2){
				cv2->size = append_cigars(cv2->buffer, cv2->size, iter->type1, iter->len);
			} else if(iter->type1 == ALN_CIGAR_TYPE_MAT){
				if(iter->type2 == ALN_CIGAR_TYPE_DEL){
					cv2->size = append_cigars(cv2->buffer, cv2->size, ALN_CIGAR_TYPE_DEL, iter->len);
					iter->l1  += iter->len;
				} else {
					cv2->size = append_cigars(cv2->buffer, cv2->size, ALN_CIGAR_TYPE_MAT, iter->len);
				}
			} else if(iter->type1 == ALN_CIGAR_TYPE_DEL){
				cv2->size = append_cigars(cv2->buffer, cv2->size, ALN_CIGAR_TYPE_DEL, iter->len);
				iter->l2 += iter->len;
			} else {
				fprintf(stdout, " -- G%d in %s -- %s:%d --\n", g->cns_id, __FUNCTION__, __FILE__, __LINE__);
				fflush(stdout);
				abort();
			}
		}
		clear_cigarv(cv0); append_cigarv(cv0, cv2);
		//fprintf(stdout, "G%d\tC\t%s\n", g->cns_id, cigars2string(cv0->buffer, cv0->size, g->chs[2])); fflush(stdout);
	}
	cns_len = 0;
	for(i=0;i<cv0->size;i++) cns_len += cv0->buffer[i].len;
	for(i=0;i<g->msa[1]->size;i++){ clear_b8list(get_strm(g->msa[1], i)); encap_b8list(get_strm(g->msa[1], i), cns_len + 1); }
	for(i=g->msa[1]->size;i<=g->sp_hits->size;i++) push_strm(g->msa[1], init_b8list(cns_len + 1));
	cns = get_strm(g->msa[1], g->sp_hits->size);
	cns->size = cigars_seq2aln((char*)cns->buffer, cv0->buffer, cv0->size, 0, g->skeleton->string);
	for(i=0;i<g->sp_hits->size;i++){
		sp = ref_sr_hitv(g->sp_hits, i);
		cv1 = get_cigarm(g->cigars[1], i);
		clear_cigarv(cv1);
		str = get_strm(g->msa[1], i);
		reset_iter_2cigars(iter, cv0->buffer, cv0->size, sp->cigars, sp->n_cigar);
		while(iter_2cigars(iter)){
			if(iter->type1 == ALN_CIGAR_TYPE_MAT) cv1->size = append_cigars(cv1->buffer, cv1->size, iter->type2, iter->len);
			else if(iter->type2 == ALN_CIGAR_TYPE_DEL) cv1->size = append_cigars(cv1->buffer, cv1->size, ALN_CIGAR_TYPE_MAT, iter->len);
			else { cv1->size = append_cigars(cv1->buffer, cv1->size, ALN_CIGAR_TYPE_INS, iter->len); iter->l2 += iter->len; }
		}
		if(sp->dir2) bits2revseq(g->chs[0], ref_u64list(g->sp_sdb->rd_seqs, 0), sr_rdseq_offset(g->sp_sdb, sp->rid2), sr_rdseq_length(g->sp_sdb, sp->rid2));
		else bits2seq(g->chs[0], ref_u64list(g->sp_sdb->rd_seqs, 0), sr_rdseq_offset(g->sp_sdb, sp->rid2), sr_rdseq_length(g->sp_sdb, sp->rid2));
		str->size = cigars_seq2aln((char*)str->buffer, cv1->buffer, cv1->size, 1, g->chs[0]);
	}
	free_iter_2cigars(iter);
	clear_u32list(g->bcovs);
	clear_string(g->cns_seq);
	g->min_bcov = 1000;
	for(i=0;i<cns_len;i++){
		acgt[0] = 0; acgt[1] = 0; acgt[2] = 0; acgt[3] = 0; acgt[4] = 0; acgt[5] = 0;
		for(j=0;j<g->sp_hits->size;j++){
			str = get_strm(g->msa[1], j);
			if(str->buffer[i] == ' ') continue;
			if(str->buffer[i] == '-') acgt[5] ++;
			else {
				acgt[base_bit_table[str->buffer[i]]] ++;
				acgt[4] ++;
			}
		}
		base = 5;
		for(j=0;j<4;j++) if(acgt[j] > acgt[base]){ base = j; }
		push_u32list(g->bcovs, acgt[4]);
		push_u32list(g->pcovs, 0);
		if(acgt[base] > 1) cns->buffer[i] = bit_base_table[base];
		if(base < 4 || acgt[base] == 0) add_char_string(g->cns_seq, cns->buffer[i]);
		if(base != 5 && i >= g->sdb->max_rd_len && i + g->sdb->max_rd_len < cns_len && g->min_bcov > acgt[base]) g->min_bcov = acgt[base];
	}
	if(g->min_bcov == 1000) g->min_bcov = 0;
	g->cns_len = g->cns_seq->size;
	for(i=0;i<g->sp_hits->size;i++){
		str = get_strm(g->msa[1], i);
		for(j=0;j<cns_len;j++){
			if(str->buffer[j] == '-') continue;
			if(str->buffer[j] == ' '){ str->buffer[j] = '-'; continue; }
			if(str->buffer[j] != cns->buffer[j]) str->buffer[j] += 'a' - 'A';
		}
	}
	clear_u32list(g->pcovs);
	for(i=0;i<g->sp_hits->size;i+=2){
		sp1 = ref_sr_hitv(g->sp_hits, i + 0);
		sp2 = ref_sr_hitv(g->sp_hits, i + 1);
		if(sp1->off < sp2->off){
			b = sp1->off; e = sp2->off + sr_rdseq_length(g->sp_sdb, sp2->rid2);
		} else {
			b = sp2->off; e = sp1->off + sr_rdseq_length(g->sp_sdb, sp1->rid2);
		}
		for(j=b;j<e;j++) g->pcovs->buffer[j] ++;
	}
	b = sr_rdseq_length(g->sdb, 0);
	e = cns_len - sr_rdseq_length(g->sdb, 1) + 1;
	if(b + 10 > e){ b = 0; e = cns_len; }
	g->min_pcov = 1000;
	for(i=b;i<e;i++){
		if(get_u32list(g->pcovs, i) < g->min_pcov) g->min_pcov = get_u32list(g->pcovs, i); 
	}
}

void output_lgraph(LGraph *g, FILE *seqf, FILE *msaf){
	uint32_t i;
	SR_AlnHit *sp;
	Layout *lay;
	Node *n;
	cigarv *cv;
	b8list *str, *cns;
	if(seqf){
		fprintf(seqf, ">ps%u len=%u n_anchor=%d,%d n_cns=%d n_realn=%d min_base_cov=%u min_mate_cov=%u\n", g->cns_id, g->cns_len,
				g->n_l, g->n_r, g->n_layout, (unsigned)g->ords[1]->size,  g->min_bcov, g->min_pcov);
		for(i=0;i<g->cns_len;i++){
			fputc(g->cns_seq->string[i], seqf);
			if(((i + 1) % 100) == 0) fputc('\n', seqf);
		}
		if(i % 100) fputc('\n', seqf);
	}
	if(msaf){
		fprintf(msaf, ">ps%u len=%u n_anchor=%d,%d n_cns=%d n_realn=%d min_base_cov=%u min_mate_cov=%u\n", g->cns_id, g->cns_len,
				g->n_l, g->n_r, g->n_layout, (unsigned)g->ords[1]->size,  g->min_bcov, g->min_pcov);
		fprintf(msaf, "skeleton\t+\t%s\n", g->skeleton->string);
		for(i=0;i<g->lays->size;i++){
			lay = ref_layv(g->lays, i);
			n = ref_nodev(g->nodes, lay->nid);
			cv = get_cigarm(g->cigars[0], i);
			str = get_strm(g->msa[0], i);
			fprintf(msaf, "%010u\t%c\t%s\t%s\n", n->rid, "+-"[n->seqdir], str->buffer, cigars2string(cv->buffer, cv->size, g->chs[1]));
		}
#define cns_lg_hit_off(x) ((int)(ref_sr_hitv(g->sp_hits, x)->off))
		sort_array(g->ords[1]->buffer, g->ords[1]->size, uint32_t, ((int)(cns_lg_hit_off(a) - cns_lg_hit_off(b))));
		fprintf(msaf, "final_seq\t+\t%s\n", g->cns_seq->string);
		cns = get_strm(g->msa[1], g->sp_hits->size);
		cv = get_cigarm(g->cigars[1], g->sp_hits->size);
		fprintf(msaf, "full_cns\t+\t%s\t%s\n", cns->buffer, cigars2string(cv->buffer, cv->size, g->chs[0]));
		for(i=0;i<g->ords[1]->size;i++){
			sp = ref_sr_hitv(g->sp_hits, get_u32list(g->ords[1], i));
			cv = get_cigarm(g->cigars[1], i);
			str = get_strm(g->msa[1], get_u32list(g->ords[1], i));
			fprintf(msaf, "%010u\t%c\t%s\t%s\n", sp->rid2, "+-"[sp->dir2], str->buffer, cigars2string(cv->buffer, cv->size, g->chs[1]));
		}
	}
}

